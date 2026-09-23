import os
import re

import numpy as np
import pandas as pd
from astropy.table import Table

from tudatpy import constants
from tudatpy.astro import time_representation
from tudatpy.data_input.tracking_data.optical_utilities import SPACECRAFT_POSITION_COLUMNS
from tudatpy.data_input.tracking_data.radar_utilities import (
    RADAR_TABLE_META_KEY,
    empty_radar_table,
    radar_data_from_raw,
)

from . import unpackers

PARSED_80COL_COLUMNS = [
    "number",
    "provisional_designation",
    "discovery",
    "epoch",
    "epoch_seconds_UTC",
    "RA",
    "DEC",
    "observatory",
    "magnitude",
    "band",
    "note1",
    "note2",
    "catalog",
    "spacecraft_parallax_type",
    *SPACECRAFT_POSITION_COLUMNS,
]


def _split_80_column_records(lines: list[str]) -> list[str]:
    """Strip line endings and split concatenated records."""
    records = []
    for line in (str(line).rstrip("\r\n") for line in lines):
        if len(line) > 80 and len(line) % 80 == 0:
            records.extend(line[index : index + 80] for index in range(0, len(line), 80))
        else:
            records.append(line)
    return records


def _parse_implicit_decimal(field: str, integer_width: int) -> float:
    """Parse a fixed-width field with an implied decimal point."""
    if not field.strip():
        return np.nan
    return float(f"{field[:integer_width]}.{field[integer_width:]}".replace(" ", "0"))


def _mpc_epochs_utc(year, month, day_fraction, round_to_second=False) -> np.ndarray:
    """Return UTC seconds since J2000 from MPC date fields."""
    day = day_fraction.astype(int)
    epochs = (
        np.array(
            [
                time_representation.date_time_components_to_epoch(
                    int(current_year), int(current_month), int(current_day), 0, 0, 0.0
                )
                for current_year, current_month, current_day in zip(year, month, day)
            ]
        )
        + (day_fraction - day).to_numpy(dtype=float) * constants.JULIAN_DAY
    )
    return np.round(epochs) if round_to_second else epochs


def _unpacked_names(number: pd.Series, provisional_designation: pd.Series) -> pd.Series:
    """Return human-readable object names from packed MPC designations."""
    identification = pd.DataFrame(
        {"number": number, "provisional_designation": provisional_designation}
    ).apply(identify_object, axis=1)
    return identification["unpacked_number"].fillna(identification["unpacked_name"])


def _parse_radar_records(lines: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    """Parse MPC radar record pairs (Note 2 R followed by r)."""
    is_first = lines.str[14] == "R"
    if not is_first.any():
        return empty_radar_table(), is_first
    first = lines[is_first]
    second = lines.shift(-1)[is_first]
    not_paired = second.str[14] != "r"
    if not_paired.any():
        raise ValueError(
            f"Radar Structure Error at Line {not_paired.idxmax() + 1}: "
            "'R' record not followed by an 'r' record."
        )

    def decimal(field: pd.Series, integer_width: int) -> pd.Series:
        return field.map(lambda text: _parse_implicit_decimal(text, integer_width))

    raw = pd.DataFrame(
        {
            "target_body": _unpacked_names(first.str[0:5], first.str[5:12]),
            "epoch_seconds_UTC": _mpc_epochs_utc(
                first.str[15:19].astype(int),
                first.str[20:22].astype(int),
                first.str[23:32].astype(float),
                round_to_second=True,
            ),
            "transmitter": first.str[68:71].str.strip().str.zfill(3),
            "receiver": first.str[77:80].str.strip().str.zfill(3),
            "target_point": second.str[32],
            "transmitter_frequency_hz": decimal(first.str[62:68], 5) * 1.0e6,
            "delay_us": decimal(first.str[32:47], 11),
            "delay_sigma_us": decimal(second.str[33:47], 10),
            "doppler_hz": decimal(first.str[47:62], 11),
            "doppler_sigma_hz": decimal(second.str[47:62], 11),
        },
        index=first.index,
    )
    if second.str[62:68].str.strip().ne("").any():
        raise ValueError(
            "MPC radar records with a transmitter frequency continued on the "
            "'r' record are not supported."
        )
    return radar_data_from_raw(raw, source="MPC"), is_first | is_first.shift(1, fill_value=False)


def _parse_spacecraft_parallax(parallax_lines: pd.Series) -> pd.DataFrame:
    """Return geocentric J2000 spacecraft positions [m] from MPC s records."""
    parallax_type = pd.to_numeric(parallax_lines.str[32], errors="coerce")
    positions = pd.DataFrame(
        {
            column: pd.to_numeric(parallax_lines.str[field].str.replace(" ", ""), errors="coerce")
            for column, field in zip(
                SPACECRAFT_POSITION_COLUMNS,
                [slice(34, 45), slice(46, 57), slice(58, 69)],
            )
        },
        index=parallax_lines.index,
    )
    invalid = ~parallax_type.isin([1, 2]) | positions.isna().any(axis=1)
    if invalid.any():
        raise ValueError(
            f"Satellite Structure Error at Line {invalid.idxmax() + 1}.\n"
            "Parallax line needs type '1' or '2' in column 33 and three numeric components."
        )
    positions = positions.mul(
        np.where(parallax_type == 1, 1.0e3, constants.ASTRONOMICAL_UNIT), axis=0
    )
    positions["spacecraft_parallax_type"] = parallax_type
    return positions


def get_first_failure_reason(row: pd.Series) -> str:
    """
    Returns ONLY the first reason why a row failed validation.

    Parameters
    ----------
    row : pd.Series
        A row from the DataFrame representing a single observation line.

    Returns
    -------
    str
        A descriptive error message indicating the specific validation failure.
    """

    # 1. Check Mandatory Internal Separators (RA/DEC internals)
    sep_checks = {
        "sep_ra_hm": "Space between RA Hour/Min",
        "sep_ra_ms": "Space between RA Min/Sec",
        "sep_dec_dm": "Space between DEC Deg/Min",
        "sep_dec_ms": "Space between DEC Min/Sec",
    }
    for col, name in sep_checks.items():
        val = str(row[col])
        if val != " ":
            return f"Invalid Separator: {name} contains '{val}' (expected space)"

    # 2. Check for Parsing Failures
    numeric_checks = {
        "year": "Year",
        "month": "Month",
        "day_frac": "Day",
        "ra_h": "RA Hours",
        "ra_m": "RA Minutes",
        "ra_s": "RA Seconds",
        "dec_d": "DEC Degrees",
        "dec_m": "DEC Minutes",
        "dec_s": "DEC Seconds",
    }

    for col, name in numeric_checks.items():
        original_str = str(row[col]) if pd.notna(row[col]) else ""
        if pd.isna(row[f"{col}_n"]) and original_str.strip() != "":
            return f"Invalid format in {name} ('{original_str}')"

    # 3. Check Logic/Ranges
    if pd.notna(row["month_n"]) and not (1 <= row["month_n"] <= 12):
        return f"Month '{int(row['month_n'])}' out of range (1-12)"
    if pd.notna(row["ra_h_n"]) and not (0 <= row["ra_h_n"] <= 24):
        return f"RA Hours '{row['ra_h_n']}' out of range (0-24)"
    if pd.notna(row["ra_m_n"]) and row["ra_m_n"] >= 60:
        return f"RA Minutes '{row['ra_m_n']}' >= 60"
    if pd.notna(row["ra_s_n"]) and row["ra_s_n"] >= 60:
        return f"RA Seconds '{row['ra_s_n']}' >= 60"
    if pd.notna(row["dec_m_n"]) and row["dec_m_n"] >= 60:
        return f"DEC Minutes '{row['dec_m_n']}' >= 60"
    if pd.notna(row["dec_s_n"]) and row["dec_s_n"] >= 60:
        return f"DEC Seconds '{row['dec_s_n']}' >= 60"

    # 4. Check Context (is something is wrong and the note IS NOT s, fail)
    note2 = str(row["note2"]) if pd.notna(row["note2"]) else ""
    if note2.lower() != "s":
        return f"Line is invalid and lacks 's' flag (note 2 is '{note2}'), which means this line is not a parallax line."
    return "Unknown validation error"


def parse_80cols_data(lines: list[str]) -> Table:
    """
    Parses MPC observation data in the ASCII 80-column format.

    The function uses vectorized Pandas operations for efficiency. The input
    records may contain optical observations, space-based S/s pairs, and radar
    R/r pairs. Radar observations are stored as a canonical pandas table in
    ``table.meta[RADAR_TABLE_META_KEY]``.

    Parameters
    ----------
    lines : list[str]
        Raw 80-column observation lines.

    Returns
    -------
    Table
        Astropy Table with 'number' column already unpacked to human-readable format.
        (e.g., '0003I' -> '3I', '00433' -> '433')

    Raises
    ------
    ValueError
        If the file/list is empty, if lines are not 80 columns wide, or if
        validation/satellite logic fails.
    """

    # 1. LOAD RAW LINES & CHECK LENGTH
    if not lines:
        raise ValueError("Input list is empty.")

    df = pd.DataFrame({"clean_line": _split_80_column_records(lines)})
    line_length = df["clean_line"].str.len()
    if (line_length != 80).any():
        bad_idx = line_length.ne(80).idxmax()
        raise ValueError(
            f"Line Length Error at Line {bad_idx + 1}.\n"
            f"Expected 80 characters, got {line_length[bad_idx]}.\n"
            f"Content: '{df.at[bad_idx, 'clean_line']}'"
        )

    radar_table, radar_line_mask = _parse_radar_records(df["clean_line"])
    df = df.loc[~radar_line_mask].copy()

    # 2. SLICE COLUMNS
    col_map = {
        "number": slice(0, 5),
        "provisional_designation": slice(5, 12),
        "discovery": slice(12, 13),
        "note1": slice(13, 14),
        "note2": slice(14, 15),
        "year": slice(15, 19),
        "month": slice(20, 22),
        "day_frac": slice(23, 32),
        "ra_h": slice(32, 34),
        "ra_m": slice(35, 37),
        "ra_s": slice(38, 44),
        "dec_sign": slice(44, 45),
        "dec_d": slice(45, 47),
        "dec_m": slice(48, 50),
        "dec_s": slice(51, 56),
        "gap_1": slice(56, 65),
        "magnitude": slice(65, 70),
        "band": slice(70, 71),
        "gap_2": slice(71, 77),
        "observatory": slice(77, 80),
    }

    sep_map = {
        "sep_ra_hm": slice(34, 35),
        "sep_ra_ms": slice(37, 38),
        "sep_dec_dm": slice(47, 48),
        "sep_dec_ms": slice(50, 51),
    }

    for name, sl in {**col_map, **sep_map}.items():
        df[name] = df["clean_line"].str[sl]

    # 3. NUMERIC COERCION
    cols_to_convert = [
        "year",
        "month",
        "day_frac",
        "ra_h",
        "ra_m",
        "ra_s",
        "dec_d",
        "dec_m",
        "dec_s",
        "magnitude",
    ]

    for col in cols_to_convert:
        df[f"{col}_n"] = pd.to_numeric(df[col], errors="coerce")

    # 4. SATELLITE PAIRING VALIDATION
    flag_series = df["note2"].fillna("")
    is_sat_obs = flag_series == "S"
    is_sat_par = flag_series == "s"

    count_obs = is_sat_obs.sum()
    count_par = is_sat_par.sum()

    if count_obs != count_par:
        raise ValueError(
            f"Satellite Structure Error: Mismatch in satellite lines.\n"
            f"Found {count_obs} Observation lines ('S') and {count_par} Parallax lines ('s')."
        )

    if count_obs > 0:
        # MPC spacecraft observations are represented by an observation record
        # ('S') immediately followed by a parallax vector record ('s'). The
        # vector record is not an observation by itself and is joined onto the
        # preceding optical observation row.
        next_is_s = flag_series.shift(-1) == "s"
        valid_pairs = is_sat_obs & next_is_s
        if valid_pairs.sum() != count_obs:
            bad_indices = df.index[is_sat_obs & (~next_is_s)]
            raise ValueError(
                f"Satellite Structure Error at Line {bad_indices[0] + 1}.\n"
                f"Observation 'S' not followed by Parallax 's'."
            )

    spacecraft_positions = _parse_spacecraft_parallax(df.loc[is_sat_par, "clean_line"]).set_axis(
        df.index[is_sat_obs]
    )

    # 5. VALIDATION LOGIC
    is_valid_structure = (
        (df["sep_ra_hm"] == " ")
        & (df["sep_ra_ms"] == " ")
        & (df["sep_dec_dm"] == " ")
        & (df["sep_dec_ms"] == " ")
    )

    is_valid_data = (
        df["month_n"].between(1, 12)
        & df["day_frac_n"].notna()
        & df["ra_h_n"].between(0, 24)
        & (df["ra_m_n"] < 60)
        & (df["ra_s_n"] < 60)
        & (df["dec_m_n"] < 60)
        & (df["dec_s_n"] < 60)
    )

    is_valid_obs = is_valid_structure & is_valid_data
    is_satellite_flag = is_sat_par

    if hasattr(unpackers, "OBS_TYPES_TO_DROP"):
        parser_drop_flags = [flag for flag in unpackers.OBS_TYPES_TO_DROP if flag != "S"]
        is_drop_flag = df["note2"].isin(parser_drop_flags)
    else:
        is_drop_flag = pd.Series(False, index=df.index)

    mask_error = (~is_valid_obs) & (~is_satellite_flag) & (~is_drop_flag)

    if mask_error.any():
        first_idx = df[mask_error].index[0]
        bad_row = df.loc[first_idx]
        error_msg = get_first_failure_reason(bad_row)
        raise ValueError(
            f"Parsing Error at Line {first_idx + 1}.\n"
            f"Reason: {error_msg}\n"
            f"Line Content: '{bad_row['clean_line']}'"
        )

    df_obs = df[is_valid_obs & (~is_drop_flag)].join(spacecraft_positions, how="left")
    if df_obs.empty and radar_table.empty:
        raise ValueError("No valid observation lines found.")
    final_df = (
        _optical_output_frame(df_obs)
        if not df_obs.empty
        else pd.DataFrame(columns=PARSED_80COL_COLUMNS)
    )
    parsed_table = Table.from_pandas(final_df)
    if not radar_table.empty:
        parsed_table.meta[RADAR_TABLE_META_KEY] = radar_table
    return parsed_table


def _optical_output_frame(df_obs: pd.DataFrame) -> pd.DataFrame:
    """Build the parsed optical table from validated fixed-width fields."""
    df_obs = df_obs.copy()
    for column in [
        "number",
        "provisional_designation",
        "discovery",
        "note1",
        "note2",
        "band",
        "observatory",
    ]:
        df_obs[column] = df_obs[column].str.strip().replace({"": None, np.nan: None})
    names = _unpacked_names(df_obs["number"], df_obs["provisional_designation"])
    df_obs["number"] = names.fillna(df_obs["number"])

    epochs_utc = _mpc_epochs_utc(df_obs["year_n"], df_obs["month_n"], df_obs["day_frac_n"])
    ra_rad = np.deg2rad(
        (df_obs["ra_h_n"] + df_obs["ra_m_n"] / 60.0 + df_obs["ra_s_n"] / 3600.0) * 15.0
    )
    dec_sign = np.where(df_obs["dec_sign"] == "-", -1.0, 1.0)
    dec_rad = np.deg2rad(
        (df_obs["dec_d_n"] + df_obs["dec_m_n"] / 60.0 + df_obs["dec_s_n"] / 3600.0) * dec_sign
    )

    return pd.DataFrame(
        {
            "number": df_obs["number"],
            "provisional_designation": df_obs["provisional_designation"],
            "discovery": df_obs["discovery"].eq("*"),
            "epoch": [
                time_representation.seconds_since_epoch_to_julian_day(epoch) for epoch in epochs_utc
            ],
            "epoch_seconds_UTC": epochs_utc,
            "RA": ra_rad,
            "DEC": dec_rad,
            "observatory": df_obs["observatory"],
            "magnitude": df_obs["magnitude_n"],
            "band": df_obs["band"],
            "note1": df_obs["note1"],
            "note2": df_obs["note2"],
            "catalog": None,
            "spacecraft_parallax_type": df_obs["spacecraft_parallax_type"],
            **{column: df_obs[column] for column in SPACECRAFT_POSITION_COLUMNS},
        },
        index=df_obs.index,
    )[PARSED_80COL_COLUMNS]


def parse_80cols_file(filename: str | list[str]) -> Table:
    """Parse MPC 80-column observation files into an astropy table.

    This is a supporting parser used by :func:`read_80_column_data`. In the
    typical Tudat workflow, call :func:`read_80_column_data` instead, so the
    parsed observations are converted to Tudat tracking-data objects.

    Parameters
    ----------
    filename : str | list[str]
        Path to one MPC 80-column file, or paths to multiple files. Each record
        must follow the MPC fixed-width optical, space-based or radar format.

    Returns
    -------
    Table
        Astropy table with standardized optical astrometry columns, including
        object identifiers, observation epochs, right ascension, declination,
        observatory code, magnitude, band, MPC note fields, and catalog field.
        Canonical radar data are stored in
        ``table.meta[RADAR_TABLE_META_KEY]``.
    """
    all_lines = []

    # Handle single string input
    if isinstance(filename, str):
        filepaths = [filename]
    else:
        filepaths = filename

    # Read files
    for path in filepaths:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Observation file not found: {path}")

        with open(path, "r") as f:
            all_lines.extend(f.readlines())

    return parse_80cols_data(all_lines)


def identify_object(row: pd.Series) -> pd.Series:
    """
    Internal helper to apply unpacking logic row-by-row.
    Returns unpacked_number (preferred for asteroids) and unpacked_name (for others).
    """
    # Safely extract strings
    raw_number = row["number"]
    perm_id = str(raw_number).strip() if pd.notna(raw_number) and raw_number else ""

    raw_prov = row["provisional_designation"]
    prov_id = str(raw_prov).strip() if pd.notna(raw_prov) and raw_prov else ""

    result = {"obj_type": "Unknown", "unpacked_name": None, "unpacked_number": None}

    # --- PATH A: PERMANENT ID IS PRESENT ---
    if perm_id:
        # Type 1: Natural Satellite (e.g., J013S)
        if re.match(r"^[JSUND]\d{3}S$", perm_id):
            result["obj_type"] = "Natural Satellite"
            if perm_id[0] in unpackers.PLANET_MAP:
                result["unpacked_name"] = unpackers.unpack_permanent_natural_satellite(perm_id)

        # Type 2: Comets (e.g., 0029P)
        elif re.match(r"^\d{4}[PD]$", perm_id):
            result["obj_type"] = "Comet"
            # Extract number part for potential use
            num_val = int(perm_id[0:4])
            result["unpacked_number"] = str(num_val)
            result["unpacked_name"] = f"{num_val}{perm_id[4]}"

        # Type 3: Interstellar (e.g., 0003I)
        elif re.match(r"^\d{4}I$", perm_id):
            result["obj_type"] = "Interstellar"
            # Interstellars are often referred to by Name/Designation (3I)
            result["unpacked_name"] = f"{int(perm_id[0:4])}I"

        # Type 4: Minor Planets (Standard or Packed)
        else:
            result["obj_type"] = "Minor Planet"
            # We prioritize unpacked_number for asteroids (e.g., "433")
            result["unpacked_number"] = unpackers.unpack_permanent_minor_planet(perm_id)
            result["unpacked_name"] = f"({result['unpacked_number']})"

    # --- PATH B: ONLY PROVISIONAL ID IS PRESENT ---
    elif prov_id:
        try:
            if len(prov_id) == 7 and prov_id[6].isalpha() and prov_id[6] not in ["I", "Z"]:
                result["obj_type"] = "Minor Planet"
                result["unpacked_name"] = unpackers.unpack_provisional_minor_planet(prov_id)
            else:
                result["obj_type"] = "Comet/Satellite"
                result["unpacked_name"] = unpackers.unpack_provisional_comet_or_satellite(prov_id)
        except ValueError:
            result["obj_type"] = "Unknown"
            result["unpacked_name"] = prov_id

    else:
        raise ValueError(
            f"Invalid observation line:\n {row} \n Missing both permanent ID and provisional ID."
        )

    return pd.Series(result)


def enrich_observations(observations: Table) -> Table:
    """
    Takes the raw parsed table and adds 'obj_type', 'unpacked_name', and 'unpacked_number'.

    This acts as the bridge between the raw MPC 80-col format and
    human-readable designations.

    Parameters
    ----------
    observations : Table
        The output from `parse_80cols_file`.

    Returns
    -------
    Table
        The original table with added columns.
    """

    df = observations.to_pandas()
    enrichment = df.apply(identify_object, axis=1, result_type="expand")
    df_enriched = pd.concat([df, enrichment], axis=1)

    return Table.from_pandas(df_enriched)


def parse_packed_permanent_designation(packed_perm_num: str) -> dict[str, str]:
    """
    Parses a packed permanent designation string from the MPC.

    Parameters
    ----------
    packed_perm_num : str
        The packed permanent designation string (e.g., '00433', 'J013S', '0029P').

    Returns
    -------
    dict
        A dictionary containing parsed identification data, such as 'type', 'name', 'number', and 'comettype'.
    """
    ident_data = {}
    packed_perm_num = packed_perm_num.strip()

    # Rule for Natural Satellites (e.g., J013S)
    if re.match(r"^[JSUND]\d{3}S$", packed_perm_num):
        ident_data["type"] = "Natural Satellite"
        ident_data["name"] = unpackers.unpack_permanent_natural_satellite(packed_perm_num)

    # Rule for Comets (e.g., 0029P)
    elif re.match(r"^\d{4}[PD]$", packed_perm_num):
        ident_data["type"] = "Comet"
        ident_data["number"] = str(int(packed_perm_num[0:4]))
        ident_data["comettype"] = packed_perm_num[4]

    # Rule for Interstellar Objects (e.g., 0002I)
    elif re.match(r"^\d{4}I$", packed_perm_num):
        ident_data["type"] = "Interstellar Object"
        number = int(packed_perm_num[0:4])
        ident_data["name"] = f"{number}I"
        ident_data["number"] = str(number)
        ident_data["comettype"] = "I"

    # Rule for Minor Planets (e.g., 00433, A0345, D4341, ~000z)
    elif packed_perm_num:
        ident_data["type"] = "Minor Planet"
        ident_data["number"] = unpackers.unpack_permanent_minor_planet(packed_perm_num)

    return ident_data


def parse_80cols_identification_fields(line: str) -> dict[str, str]:
    """
    Parses the identification part of an 80-column MPC observation line.

    This function extracts and interprets the packed permanent designation,
    provisional designation, and discovery flag from the beginning of an MPC line.

    Parameters
    ----------
    line : str
        A single 80-column string representing an MPC observation.

    Returns
    -------
    dict
        A dictionary containing parsed identification data, such as 'type', 'number', 'desig', and 'discovery'.
    """
    ident_data = {}
    packed_perm_num = line[0:5].strip()
    packed_prov_desig = line[5:12].strip()

    if line[12].strip() == "*":
        ident_data["discovery"] = True

    if packed_perm_num:
        ident_data.update(parse_packed_permanent_designation(packed_perm_num))
        if packed_prov_desig:
            if ident_data.get("type") == "Minor Planet":
                ident_data["desig"] = unpackers.unpack_provisional_minor_planet(packed_prov_desig)
            elif ident_data.get("type") in ["Comet", "Interstellar Object"]:
                ident_data["desig"] = packed_prov_desig.strip() or "NaN"
    elif packed_prov_desig:
        if packed_prov_desig[6].isalpha() and packed_prov_desig[6] not in ["I", "Z"]:
            ident_data["type"] = "Minor Planet"
            ident_data["desig"] = unpackers.unpack_provisional_minor_planet(packed_prov_desig)
        else:
            ident_data["type"] = "Comet"
            ident_data["desig"] = unpackers.unpack_provisional_comet_or_satellite(packed_prov_desig)

    return ident_data
