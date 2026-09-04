import numpy as np

# CCSDS 502.0-B-3 Table 4-3 (OMM) / Table 3-3 (OPM): the 21 named
# lower-triangular covariance keywords, each mapped to its (row, col) in
# the 6x6 matrix. Unlike OEM/OCM covariance data lines (positional rows of
# numbers, handled by reconstruct_covariance/make_lower_triangular below),
# the OMM/OPM covariance block is itself KVN -- each matrix cell is its
# own named "KEYWORD = value" line -- so each keyword maps directly to a
# known cell instead of being inferred from row/column position.
_NAMED_COVARIANCE_KEYWORD_CELLS = [
    ("CX_X", 0, 0),
    ("CY_X", 1, 0),
    ("CY_Y", 1, 1),
    ("CZ_X", 2, 0),
    ("CZ_Y", 2, 1),
    ("CZ_Z", 2, 2),
    ("CX_DOT_X", 3, 0),
    ("CX_DOT_Y", 3, 1),
    ("CX_DOT_Z", 3, 2),
    ("CX_DOT_X_DOT", 3, 3),
    ("CY_DOT_X", 4, 0),
    ("CY_DOT_Y", 4, 1),
    ("CY_DOT_Z", 4, 2),
    ("CY_DOT_X_DOT", 4, 3),
    ("CY_DOT_Y_DOT", 4, 4),
    ("CZ_DOT_X", 5, 0),
    ("CZ_DOT_Y", 5, 1),
    ("CZ_DOT_Z", 5, 2),
    ("CZ_DOT_X_DOT", 5, 3),
    ("CZ_DOT_Y_DOT", 5, 4),
    ("CZ_DOT_Z_DOT", 5, 5),
]


def covariance_matrix_from_named_keywords(data: dict) -> np.ndarray | None:
    """
    Builds a 6x6 symmetric covariance matrix from CCSDS's named
    lower-triangular covariance keywords (CX_X, CY_X, CY_Y, ...,
    CZ_DOT_Z_DOT -- see `_NAMED_COVARIANCE_KEYWORD_CELLS`).

    Parameters
    ----------
    data : dict
        Flat keyword -> value mapping (e.g. an already-parsed OMM's data
        section). Only the 21 known covariance keys are read; anything
        else in `data` is ignored. Values may be strings (as parsed from
        KVN) or numbers.

    Returns
    -------
    np.ndarray | None
        A 6x6 symmetric matrix, or None if NONE of the 21 covariance
        keywords are present in `data` (covariance is entirely optional).

    Raises
    ------
    ValueError
        Per CCSDS 502.0-B-3 4.2.4.5, covariance is all-or-nothing: if ANY
        of the 21 named keywords are present, all 21 must be. Raised if
        only some are present.
    """
    present_keys = [key for key, _, _ in _NAMED_COVARIANCE_KEYWORD_CELLS if key in data]
    if not present_keys:
        return None
    if len(present_keys) != len(_NAMED_COVARIANCE_KEYWORD_CELLS):
        missing_keys = [
            key for key, _, _ in _NAMED_COVARIANCE_KEYWORD_CELLS if key not in present_keys
        ]
        raise ValueError(
            "OMM/OPM covariance is all-or-nothing (CCSDS 502.0-B-3 4.2.4.5): "
            f"{len(present_keys)}/{len(_NAMED_COVARIANCE_KEYWORD_CELLS)} named "
            f"covariance keywords are present; missing: {', '.join(missing_keys)}."
        )

    matrix = np.zeros((6, 6))
    for key, row, col in _NAMED_COVARIANCE_KEYWORD_CELLS:
        val = float(data[key])
        matrix[row, col] = val
        matrix[col, row] = val
    return matrix


def covariance_matrix_to_named_keywords(matrix) -> dict[str, float]:
    """
    Inverse of `covariance_matrix_from_named_keywords`: flattens a 6x6
    symmetric covariance matrix into the 21 named lower-triangular
    keywords used by the OMM/OPM covariance block.

    Parameters
    ----------
    matrix : array-like
        A 6x6 (symmetric) covariance matrix.

    Returns
    -------
    dict[str, float]
        `{"CX_X": ..., "CY_X": ..., ..., "CZ_DOT_Z_DOT": ...}`, reading
        values off the lower triangle (as does `make_lower_triangular`).
    """
    matrix = np.asarray(matrix, dtype=float)
    return {key: float(matrix[row, col]) for key, row, col in _NAMED_COVARIANCE_KEYWORD_CELLS}


def reconstruct_covariance(covariance_lines):
    """
    Converts a list of OEM covariance lines (strings) into a full 6x6 symmetric numpy matrix.

    Args:
        covariance_lines (list of str): The lines from the OEM file for a single epoch.
                                        Example: ["1.0", "0.1 2.0", ...]
    Returns:
        np.array: A 6x6 symmetric covariance matrix.
    """
    # Initialize a 6x6 matrix of zeros
    matrix = np.zeros((6, 6))

    # The OEM standard defines 6 rows for the lower triangle
    # We assume the input list contains exactly these rows in order.
    for i, line in enumerate(covariance_lines):
        # Parse the numbers from the string line
        values = [float(x) for x in line.split()]

        # Determine which row index we are filling (0 to 5)
        row_idx = i

        for col_idx, val in enumerate(values):
            # Fill the lower triangle
            matrix[row_idx, col_idx] = val

            # Fill the upper triangle (symmetry)
            # If row_idx == col_idx (the diagonal), this just overwrites itself, which is fine.
            matrix[col_idx, row_idx] = val

    return matrix


def make_lower_triangular(matrix):
    """
    Converts a 6x6 symmetric numpy matrix into a list of OEM lower triangular strings.

    Args:
        matrix (np.array): A 6x6 symmetric covariance matrix.

    Returns:
        list of str: The lower triangular lines formatted for an OEM file.
                     Example: ["1.0000...e+01", "1.0000...e-01 2.0000...e+01", ...]
    """
    lines = []

    # Iterate through the 6 rows defined in OEM standard
    for i in range(6):
        # Slice the row: take elements from index 0 up to and including the diagonal (i)
        # We use :i+1 because Python slices exclude the end index.
        row_values = matrix[i, : i + 1]

        # Convert floats to strings using scientific notation for precision
        # CCSDS/OEM usually requires high precision (e.g., 14 decimal places)
        str_values = [f"{val:.14e}" for val in row_values]

        # Join the values with a space and add to our list
        lines.append(" ".join(str_values))

    return lines
