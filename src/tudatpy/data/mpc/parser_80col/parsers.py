import re
import numpy as np
import pandas as pd
from astropy.table import Table
from tudatpy.astro.time_representation import DateTime
from tudatpy.data.mpc.parser_80col import unpackers
import os
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
        'sep_ra_hm': 'Space between RA Hour/Min',
        'sep_ra_ms': 'Space between RA Min/Sec',
        'sep_dec_dm': 'Space between DEC Deg/Min',
        'sep_dec_ms': 'Space between DEC Min/Sec'
    }
    for col, name in sep_checks.items():
        val = str(row[col])
        if val != ' ':
            return f"Invalid Separator: {name} contains '{val}' (expected space)"

    # 2. Check for Parsing Failures
    numeric_checks = {
        'year': 'Year', 'month': 'Month', 'day_frac': 'Day',
        'ra_h': 'RA Hours', 'ra_m': 'RA Minutes', 'ra_s': 'RA Seconds',
        'dec_d': 'DEC Degrees', 'dec_m': 'DEC Minutes', 'dec_s': 'DEC Seconds'
    }

    for col, name in numeric_checks.items():
        original_str = str(row[col]) if pd.notna(row[col]) else ""
        if pd.isna(row[f'{col}_n']) and original_str.strip() != "":
            return f"Invalid format in {name} ('{original_str}')"

    # 3. Check Logic/Ranges
    if pd.notna(row['month_n']) and not (1 <= row['month_n'] <= 12):
        return f"Month '{int(row['month_n'])}' out of range (1-12)"
    if pd.notna(row['ra_h_n']) and not (0 <= row['ra_h_n'] <= 24):
        return f"RA Hours '{row['ra_h_n']}' out of range (0-24)"
    if pd.notna(row['ra_m_n']) and row['ra_m_n'] >= 60:
        return f"RA Minutes '{row['ra_m_n']}' >= 60"
    if pd.notna(row['ra_s_n']) and row['ra_s_n'] >= 60:
        return f"RA Seconds '{row['ra_s_n']}' >= 60"
    if pd.notna(row['dec_m_n']) and row['dec_m_n'] >= 60:
        return f"DEC Minutes '{row['dec_m_n']}' >= 60"
    if pd.notna(row['dec_s_n']) and row['dec_s_n'] >= 60:
        return f"DEC Seconds '{row['dec_s_n']}' >= 60"

    # 4. Check Context (is something is wrong and the note IS NOT s, fail)
    note2 = str(row['note2']) if pd.notna(row['note2']) else ''
    if note2.lower() != 's':
        return f"Line is invalid and lacks 's' flag (note 2 is '{note2}'), which means this line is not a parallax line."
    return "Unknown validation error"

def parse_80cols_data(lines: list[str]) -> Table:
    """
    Parses MPC observation data (ASCII 80-column format) from a file or a list of strings.

    The function uses vectorized Pandas operations for efficiency.

    Parameters
    ----------
    source : str or list[str]
        The source of the data.
        - If a `str`, it is treated as a filepath to be opened and read.
        - If a `list[str]`, it is treated as the raw lines of data.

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

    df = pd.DataFrame({'raw_line': lines})
    # Remove newlines regardless of input method
    df['clean_line'] = df['raw_line'].astype(str).str.replace('\n', '', regex=False).str.replace('\r', '', regex=False)

    df['len'] = df['clean_line'].str.len()
    if (df['len'] != 80).any():
        bad_idx = df.index[df['len'] != 80][0]
        bad_row = df.iloc[bad_idx]
        raise ValueError(
            f"Line Length Error at Line {bad_idx + 1}.\n"
            f"Expected 80 characters, got {bad_row['len']}.\n"
            f"Content: '{bad_row['clean_line']}'"
        )

    # 2. SLICE COLUMNS
    col_map = {
        'number': slice(0, 5),
        'provisional_designation': slice(5, 12),
        'discovery': slice(12, 13),
        'note1': slice(13, 14),
        'note2': slice(14, 15),
        'year': slice(15, 19),
        'month': slice(20, 22),
        'day_frac': slice(23, 32),
        'ra_h': slice(32, 34),
        'ra_m': slice(35, 37),
        'ra_s': slice(38, 44),
        'dec_sign': slice(44, 45),
        'dec_d': slice(45, 47),
        'dec_m': slice(48, 50),
        'dec_s': slice(51, 56),
        'gap_1': slice(56, 65),
        'magnitude': slice(65, 70),
        'band': slice(70, 71),
        'gap_2': slice(71, 77),
        'observatory': slice(77, 80)
    }

    sep_map = {
        'sep_ra_hm': slice(34, 35),
        'sep_ra_ms': slice(37, 38),
        'sep_dec_dm': slice(47, 48),
        'sep_dec_ms': slice(50, 51),
    }

    for name, sl in {**col_map, **sep_map}.items():
        df[name] = df['clean_line'].str[sl]

    # 3. NUMERIC COERCION
    cols_to_convert = ['year', 'month', 'day_frac', 'ra_h', 'ra_m', 'ra_s',
                       'dec_d', 'dec_m', 'dec_s', 'magnitude']

    for col in cols_to_convert:
        df[f'{col}_n'] = pd.to_numeric(df[col], errors='coerce')

    # 4. SATELLITE PAIRING VALIDATION
    flag_series = df['note2'].fillna('')
    is_sat_obs = flag_series == 'S'
    is_sat_par = flag_series == 's'

    count_obs = is_sat_obs.sum()
    count_par = is_sat_par.sum()

    if count_obs != count_par:
        raise ValueError(
            f"Satellite Structure Error: Mismatch in satellite lines.\n"
            f"Found {count_obs} Observation lines ('S') and {count_par} Parallax lines ('s')."
        )

    if count_obs > 0:
        next_is_s = flag_series.shift(-1) == 's'
        valid_pairs = is_sat_obs & next_is_s
        if valid_pairs.sum() != count_obs:
            bad_indices = df.index[is_sat_obs & (~next_is_s)]
            raise ValueError(
                f"Satellite Structure Error at Line {bad_indices[0] + 1}.\n"
                f"Observation 'S' not followed by Parallax 's'."
            )

    # 5. VALIDATION LOGIC
    is_valid_structure = (
            (df['sep_ra_hm'] == ' ') &
            (df['sep_ra_ms'] == ' ') &
            (df['sep_dec_dm'] == ' ') &
            (df['sep_dec_ms'] == ' ')
    )

    is_valid_data = (
            df['month_n'].between(1, 12) &
            df['day_frac_n'].notna() &
            df['ra_h_n'].between(0, 24) & (df['ra_m_n'] < 60) & (df['ra_s_n'] < 60) &
            (df['dec_m_n'] < 60) & (df['dec_s_n'] < 60)
    )

    is_valid_obs = is_valid_structure & is_valid_data
    is_satellite_flag = is_sat_par

    if hasattr(unpackers, 'OBS_TYPES_TO_DROP'):
        is_drop_flag = df['note2'].isin(unpackers.OBS_TYPES_TO_DROP)
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

    df_obs = df[is_valid_obs & (~is_drop_flag)].copy()
    if df_obs.empty:
        raise ValueError("No valid observation lines found.")
    str_cols = ['number', 'provisional_designation', 'discovery', 'note1', 'note2', 'band', 'observatory']
    for col in str_cols:
        df_obs[col] = df_obs[col].str.strip().replace({"": None, np.nan: None})

    ident_info = df_obs.apply(identify_object, axis=1)

    human_readable_number = ident_info['unpacked_number'].fillna(ident_info['unpacked_name'])
    df_obs['number'] = human_readable_number.fillna(df_obs['number'])

    # -------------------------------------------------------------------------
    # 8. FINAL CALCULATIONS
    # -------------------------------------------------------------------------
    day_int = df_obs['day_frac_n'].astype(int)
    day_remainder = df_obs['day_frac_n'] - day_int

    timestamps = pd.to_datetime({
        'year': df_obs['year_n'],
        'month': df_obs['month_n'],
        'day': day_int
    }, errors='coerce')

    obs_time_utc = timestamps + pd.to_timedelta(day_remainder, unit='D')

    ra_deg = (df_obs['ra_h_n'] + df_obs['ra_m_n']/60.0 + df_obs['ra_s_n']/3600.0) * 15.0
    ra_rad = np.deg2rad(ra_deg)

    dec_sign_mult = np.where(df_obs['dec_sign'] == '-', -1, 1)
    dec_deg = (df_obs['dec_d_n'] + df_obs['dec_m_n']/60.0 + df_obs['dec_s_n']/3600.0) * dec_sign_mult
    dec_rad = np.deg2rad(dec_deg)

    t_obj = [DateTime.from_python_datetime(obs_time) for obs_time in obs_time_utc]

    final_df = pd.DataFrame({
        "number": df_obs['number'], # Now contains human-readable string
        "provisional_designation": df_obs['provisional_designation'],
        "discovery": df_obs['discovery'].eq('*'),
        "epoch": [t.to_julian_day() for t in t_obj],
        "epoch_utc": obs_time_utc,
        "epoch_seconds": [t.to_epoch() for t in t_obj],
        "RA": ra_rad,
        "DEC": dec_rad,
        "observatory": df_obs['observatory'],
        "magnitude": df_obs['magnitude_n'],
        "band": df_obs['band'],
        "note1": df_obs['note1'],
        "note2": df_obs['note2'],
        "catalog": None
    })
    return Table.from_pandas(final_df)


def parse_80cols_file(filename: str | list[str]) -> Table:
    """
    Reads MPC observation data from a file (or list of files) and parses it.

    Parameters
    ----------
    filename : str or list[str]
        A single file path or a list of file paths.

    Returns
    -------
    Table
        Astropy Table with unpacked data.
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

        with open(path, 'r') as f:
            all_lines.extend(f.readlines())

    return parse_80cols_data(all_lines)

# ... [identify_object and enrich_observations remain exactly the same] ...
def identify_object(row: pd.Series) -> pd.Series:
    """
    Internal helper to apply unpacking logic row-by-row.

    Returns unpacked_number (preferred for asteroids) and unpacked_name (for others).

    Parameters
    ----------
    row : pd.Series
        A row from the observations DataFrame containing 'number' and
        'provisional_designation' columns.

    Returns
    -------
    pd.Series
        A Series containing 'obj_type', 'unpacked_name', and 'unpacked_number'.
    """
    # Safely extract strings
    raw_number = row['number']
    perm_id = str(raw_number).strip() if pd.notna(raw_number) and raw_number else ''

    raw_prov = row['provisional_designation']
    prov_id = str(raw_prov).strip() if pd.notna(raw_prov) and raw_prov else ''

    result = {
        'obj_type': 'Unknown',
        'unpacked_name': None,
        'unpacked_number': None
    }

    # --- PATH A: PERMANENT ID IS PRESENT ---
    if perm_id:
        if re.match(r"^[JSUND]\d{3}S$", perm_id):
            result['obj_type'] = 'Natural Satellite'
            if perm_id[0] in unpackers.PLANET_MAP:
                result['unpacked_name'] = unpackers.unpack_permanent_natural_satellite(perm_id)

        elif re.match(r"^\d{4}[PD]$", perm_id):
            result['obj_type'] = 'Comet'
            num_val = int(perm_id[0:4])
            result['unpacked_number'] = str(num_val)
            result['unpacked_name'] = f"{num_val}{perm_id[4]}"

        elif re.match(r"^\d{4}I$", perm_id):
            result['obj_type'] = 'Interstellar'
            result['unpacked_name'] = f"{int(perm_id[0:4])}I"

        else:
            result['obj_type'] = 'Minor Planet'
            result['unpacked_number'] = unpackers.unpack_permanent_minor_planet(perm_id)
            result['unpacked_name'] = f"({result['unpacked_number']})"

    # --- PATH B: ONLY PROVISIONAL ID IS PRESENT ---
    elif prov_id:
        if len(prov_id) == 7 and prov_id[6].isalpha() and prov_id[6] not in ['I', 'Z']:
            result['obj_type'] = 'Minor Planet'
            result['unpacked_name'] = unpackers.unpack_provisional_minor_planet(prov_id)
        else:
            result['obj_type'] = 'Comet/Satellite'
            result['unpacked_name'] = unpackers.unpack_provisional_comet_or_satellite(prov_id)

    else:
        raise ValueError('Observation line does not have permanent nor provisional ID.')

    return pd.Series(result)

# -----------------------------------------------------------------------------
# 2. IDENTIFICATION & AVAILABLE ENRICHMENT (Post-Processing)
# -----------------------------------------------------------------------------

def identify_object(row: pd.Series) -> pd.Series:
    """
    Internal helper to apply unpacking logic row-by-row.
    Returns unpacked_number (preferred for asteroids) and unpacked_name (for others).
    """
    # Safely extract strings
    raw_number = row['number']
    perm_id = str(raw_number).strip() if pd.notna(raw_number) and raw_number else ''

    raw_prov = row['provisional_designation']
    prov_id = str(raw_prov).strip() if pd.notna(raw_prov) and raw_prov else ''

    result = {
        'obj_type': 'Unknown',
        'unpacked_name': None,
        'unpacked_number': None
    }

    # --- PATH A: PERMANENT ID IS PRESENT ---
    if perm_id:
        # Type 1: Natural Satellite (e.g., J013S)
        if re.match(r"^[JSUND]\d{3}S$", perm_id):
            result['obj_type'] = 'Natural Satellite'
            if perm_id[0] in unpackers.PLANET_MAP:
                result['unpacked_name'] = unpackers.unpack_permanent_natural_satellite(perm_id)

        # Type 2: Comets (e.g., 0029P)
        elif re.match(r"^\d{4}[PD]$", perm_id):
            result['obj_type'] = 'Comet'
            # Extract number part for potential use
            num_val = int(perm_id[0:4])
            result['unpacked_number'] = str(num_val)
            result['unpacked_name'] = f"{num_val}{perm_id[4]}"

        # Type 3: Interstellar (e.g., 0003I)
        elif re.match(r"^\d{4}I$", perm_id):
            result['obj_type'] = 'Interstellar'
            # Interstellars are often referred to by Name/Designation (3I)
            result['unpacked_name'] = f"{int(perm_id[0:4])}I"

        # Type 4: Minor Planets (Standard or Packed)
        else:
            result['obj_type'] = 'Minor Planet'
            # We prioritize unpacked_number for asteroids (e.g., "433")
            result['unpacked_number'] = unpackers.unpack_permanent_minor_planet(perm_id)
            result['unpacked_name'] = f"({result['unpacked_number']})"

    # --- PATH B: ONLY PROVISIONAL ID IS PRESENT ---
    elif prov_id:
        if len(prov_id) == 7 and prov_id[6].isalpha() and prov_id[6] not in ['I', 'Z']:
            result['obj_type'] = 'Minor Planet'
            result['unpacked_name'] = unpackers.unpack_provisional_minor_planet(prov_id)
        else:
            result['obj_type'] = 'Comet/Satellite'
            result['unpacked_name'] = unpackers.unpack_provisional_comet_or_satellite(prov_id)

    else:
        raise ValueError(f'Invalid observation line:\n {row} \n Missing both permanent ID and provisional ID.')

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
    enrichment = df.apply(identify_object, axis=1, result_type='expand')
    df_enriched = pd.concat([df, enrichment], axis=1)

    return Table.from_pandas(df_enriched)