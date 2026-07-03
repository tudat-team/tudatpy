import re
import numpy as np
import pandas as pd
from astropy.table import Table
from tudatpy.astro import time_representation
from tudatpy.astro.time_representation import DateTime
from tudatpy.data.mpc.parser_80col import unpackers
import os

ASTRONOMICAL_UNIT_METERS = 1.49597870700e11
SPEED_OF_LIGHT_METERS_PER_SECOND = 299792458.0


def _parse_implicit_decimal(field: str, integer_width: int) -> float:
    """Parse fixed-width MPC fields with an implicit decimal point."""
    if field.strip() == "":
        return np.nan

    padded_field = field.ljust(integer_width)
    integer_part = padded_field[:integer_width].replace(" ", "0")
    fractional_part = padded_field[integer_width:].replace(" ", "0")
    sign = ""
    if integer_part and integer_part[0] in ["+", "-"]:
        sign = integer_part[0]
        integer_part = integer_part[1:]

    return float(f"{sign}{integer_part}.{fractional_part}")


def _parse_mpc_datetime(year, month, day_frac):
    day_integer = int(float(day_frac))
    day_remainder = float(day_frac) - day_integer
    timestamp = pd.to_datetime(
        {
            "year": [int(year)],
            "month": [int(month)],
            "day": [day_integer],
        },
        errors="coerce",
    )[0]
    return timestamp + pd.to_timedelta(day_remainder, unit="D")


def _utc_datetime_to_tdb_epoch(timestamp):
    date_time = DateTime.from_iso_string(timestamp.isoformat(sep=" "))
    epoch_utc = date_time.to_epoch()
    epoch_tdb = time_representation.default_time_scale_converter().convert_time(
        input_scale=time_representation.utc_scale,
        output_scale=time_representation.tdb_scale,
        input_value=epoch_utc,
    )
    return date_time.to_julian_day(), epoch_utc, epoch_tdb


def _parse_radar_observation_pairs(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """Parse MPC two-line radar observations and return center-of-mass rows."""
    radar_rows = []
    radar_line_mask = pd.Series(False, index=df.index)
    time_scale_converter = time_representation.default_time_scale_converter()

    index = 0
    while index < len(df):
        first_record = df.iloc[index]["clean_line"]
        if len(first_record) < 15 or first_record[14] not in ["R", "V"]:
            index += 1
            continue

        if index + 1 >= len(df):
            raise ValueError(f"Radar Structure Error at Line {index + 1}: missing second record.")

        record_type = first_record[14]
        expected_second_record_type = record_type.lower()
        second_record = df.iloc[index + 1]["clean_line"]
        if len(second_record) < 33 or second_record[14] != expected_second_record_type:
            raise ValueError(
                f"Radar Structure Error at Line {index + 1}: radar '{record_type}' record "
                f"is not followed by a radar '{expected_second_record_type}' record."
            )

        radar_line_mask.iloc[index] = True
        radar_line_mask.iloc[index + 1] = True

        target_point = second_record[32]
        if target_point != "C":
            index += 2
            continue

        measurement = _parse_implicit_decimal(first_record[32:47], integer_width=11)
        measurement_sigma = _parse_implicit_decimal(second_record[33:47], integer_width=10)
        if np.isnan(measurement) or np.isnan(measurement_sigma):
            index += 2
            continue

        ident_info = identify_object(
            pd.Series(
                {
                    "number": first_record[0:5],
                    "provisional_designation": first_record[5:12],
                }
            )
        )
        number = ident_info["unpacked_number"] or ident_info["unpacked_name"]
        timestamp = _parse_mpc_datetime(
            first_record[15:19],
            first_record[20:22],
            first_record[23:32],
        )
        date_time = DateTime.from_iso_string(timestamp.isoformat(sep=" "))
        epoch_utc = date_time.to_epoch()
        epoch_tdb = time_scale_converter.convert_time(
            input_scale=time_representation.utc_scale,
            output_scale=time_representation.tdb_scale,
            input_value=epoch_utc,
        )

        transmitter = first_record[68:71].strip().zfill(3)
        receiver = first_record[77:80].strip().zfill(3)
        radar_frequency_mhz = _parse_implicit_decimal(
            first_record[62:68],
            integer_width=5,
        )
        radar_row = {
            "number": number,
            "provisional_designation": first_record[5:12].strip() or None,
            "discovery": first_record[12] == "*",
            "epoch": date_time.to_julian_day(),
            "epoch_seconds_UTC": epoch_utc,
            "epoch_seconds_TDB": epoch_tdb,
            "RA": np.nan,
            "DEC": np.nan,
            "observatory": receiver,
            "magnitude": np.nan,
            "band": None,
            "note1": first_record[13].strip() or None,
            "note2": record_type,
            "catalog": None,
            "spacecraft_parallax_type": np.nan,
            "spacecraft_position_x": np.nan,
            "spacecraft_position_y": np.nan,
            "spacecraft_position_z": np.nan,
            "radar_target_point": target_point,
            "radar_delay_us": np.nan,
            "radar_delay_sigma_us": np.nan,
            "radar_range": np.nan,
            "radar_range_sigma": np.nan,
            "radar_doppler_shift": np.nan,
            "radar_doppler_frequency": np.nan,
            "radar_doppler_frequency_sigma": np.nan,
            "radar_transmitter": transmitter,
            "radar_receiver": receiver,
            "radar_frequency_mhz": radar_frequency_mhz,
        }
        if record_type == "R":
            radar_row["radar_delay_us"] = measurement
            radar_row["radar_delay_sigma_us"] = measurement_sigma
            radar_row["radar_range"] = SPEED_OF_LIGHT_METERS_PER_SECOND * measurement * 1.0e-6
            radar_row["radar_range_sigma"] = (
                SPEED_OF_LIGHT_METERS_PER_SECOND * measurement_sigma * 1.0e-6
            )
        else:
            radar_row["radar_doppler_shift"] = measurement
            radar_row["radar_doppler_frequency"] = radar_frequency_mhz * 1.0e6 + measurement
            radar_row["radar_doppler_frequency_sigma"] = measurement_sigma
        radar_rows.append(radar_row)
        index += 2

    return pd.DataFrame(radar_rows), radar_line_mask


def _split_80_column_records(lines: list[str]) -> list[str]:
    """Split raw MPC records, including Astroquery's concatenated satellite pairs."""
    records = []
    for line in lines:
        clean_line = str(line).replace("\n", "").replace("\r", "")
        if len(clean_line) > 80 and len(clean_line) % 80 == 0:
            records.extend(clean_line[i : i + 80] for i in range(0, len(clean_line), 80))
        else:
            records.append(clean_line)
    return records


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

    df = pd.DataFrame({"raw_line": _split_80_column_records(lines)})
    # Remove newlines regardless of input method
    df["clean_line"] = (
        df["raw_line"]
        .astype(str)
        .str.replace("\n", "", regex=False)
        .str.replace("\r", "", regex=False)
    )

    df["len"] = df["clean_line"].str.len()
    if (df["len"] != 80).any():
        bad_idx = df.index[df["len"] != 80][0]
        bad_row = df.iloc[bad_idx]
        raise ValueError(
            f"Line Length Error at Line {bad_idx + 1}.\n"
            f"Expected 80 characters, got {bad_row['len']}.\n"
            f"Content: '{bad_row['clean_line']}'"
        )

    radar_table, radar_line_mask = _parse_radar_observation_pairs(df)
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

    satellite_col_map = {
        "satellite_parallax_type": slice(32, 33),
        "satellite_x": slice(34, 45),
        "satellite_y": slice(46, 57),
        "satellite_z": slice(58, 69),
    }

    sep_map = {
        "sep_ra_hm": slice(34, 35),
        "sep_ra_ms": slice(37, 38),
        "sep_dec_dm": slice(47, 48),
        "sep_dec_ms": slice(50, 51),
    }

    for name, sl in {**col_map, **sep_map, **satellite_col_map}.items():
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
        next_is_s = flag_series.shift(-1) == "s"
        valid_pairs = is_sat_obs & next_is_s
        if valid_pairs.sum() != count_obs:
            bad_indices = df.index[is_sat_obs & (~next_is_s)]
            raise ValueError(
                f"Satellite Structure Error at Line {bad_indices[0] + 1}.\n"
                f"Observation 'S' not followed by Parallax 's'."
            )

        parallax_rows = df.loc[is_sat_par].copy()
        parallax_rows["satellite_parallax_type_n"] = pd.to_numeric(
            parallax_rows["satellite_parallax_type"], errors="coerce"
        )
        if (~parallax_rows["satellite_parallax_type_n"].isin([1, 2])).any():
            bad_idx = parallax_rows.index[~parallax_rows["satellite_parallax_type_n"].isin([1, 2])][
                0
            ]
            raise ValueError(
                f"Satellite Structure Error at Line {bad_idx + 1}.\n"
                "Parallax line must specify parallax type '1' or '2' in column 33."
            )

        for component in ["satellite_x", "satellite_y", "satellite_z"]:
            parallax_rows[f"{component}_n"] = pd.to_numeric(
                parallax_rows[component].astype(str).str.replace(" ", "", regex=False),
                errors="coerce",
            )

        invalid_component = (
            parallax_rows[["satellite_x_n", "satellite_y_n", "satellite_z_n"]].isna().any(axis=1)
        )
        if invalid_component.any():
            bad_idx = parallax_rows.index[invalid_component][0]
            raise ValueError(
                f"Satellite Structure Error at Line {bad_idx + 1}.\n"
                "Could not parse one or more satellite parallax vector components."
            )

        type_1 = parallax_rows["satellite_parallax_type_n"] == 1
        scale = np.where(type_1, 1000.0, ASTRONOMICAL_UNIT_METERS)
        for component in ["satellite_x", "satellite_y", "satellite_z"]:
            parallax_rows[f"{component}_m"] = parallax_rows[f"{component}_n"] * scale

        satellite_state_columns = [
            "satellite_parallax_type_n",
            "satellite_x_m",
            "satellite_y_m",
            "satellite_z_m",
        ]
        satellite_parallax_data = parallax_rows.loc[:, satellite_state_columns].copy()
        satellite_parallax_data.index = satellite_parallax_data.index - 1
    else:
        satellite_parallax_data = pd.DataFrame(
            columns=[
                "satellite_parallax_type_n",
                "satellite_x_m",
                "satellite_y_m",
                "satellite_z_m",
            ]
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
        is_drop_flag = df["note2"].isin(unpackers.OBS_TYPES_TO_DROP)
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
    df_obs = df_obs.join(satellite_parallax_data, how="left")
    if df_obs.empty and radar_table.empty:
        raise ValueError("No valid observation lines found.")
    if df_obs.empty:
        return Table.from_pandas(radar_table)
    str_cols = [
        "number",
        "provisional_designation",
        "discovery",
        "note1",
        "note2",
        "band",
        "observatory",
    ]
    for col in str_cols:
        df_obs[col] = df_obs[col].str.strip().replace({"": None, np.nan: None})

    ident_info = df_obs.apply(identify_object, axis=1)

    human_readable_number = ident_info["unpacked_number"].fillna(ident_info["unpacked_name"])
    df_obs["number"] = human_readable_number.fillna(df_obs["number"])

    # -------------------------------------------------------------------------
    # 8. FINAL CALCULATIONS
    # -------------------------------------------------------------------------
    day_int = df_obs["day_frac_n"].astype(int)
    day_remainder = df_obs["day_frac_n"] - day_int

    timestamps = pd.to_datetime(
        {"year": df_obs["year_n"], "month": df_obs["month_n"], "day": day_int}, errors="coerce"
    )

    obs_times_utc_datetime = timestamps + pd.to_timedelta(day_remainder, unit="D")

    ra_deg = (df_obs["ra_h_n"] + df_obs["ra_m_n"] / 60.0 + df_obs["ra_s_n"] / 3600.0) * 15.0
    ra_rad = np.deg2rad(ra_deg)

    dec_sign_mult = np.where(df_obs["dec_sign"] == "-", -1, 1)
    dec_deg = (
        df_obs["dec_d_n"] + df_obs["dec_m_n"] / 60.0 + df_obs["dec_s_n"] / 3600.0
    ) * dec_sign_mult
    dec_rad = np.deg2rad(dec_deg)

    obs_times_utc_iso_string = [t.isoformat(sep=" ") for t in obs_times_utc_datetime]
    dt_objects = [DateTime.from_iso_string(t_iso) for t_iso in obs_times_utc_iso_string]
    float_epochs_utc = [dt.to_epoch() for dt in dt_objects]
    time_scale_converter = time_representation.default_time_scale_converter()
    float_epochs_tdb = [
        time_scale_converter.convert_time(
            input_scale=time_representation.utc_scale,
            output_scale=time_representation.tdb_scale,
            input_value=float_epoch_utc,
        )
        for float_epoch_utc in float_epochs_utc
    ]

    final_df = pd.DataFrame()
    if not df_obs.empty:
        final_df = pd.DataFrame(
            {
                "number": df_obs["number"],  # Now contains human-readable string
                "provisional_designation": df_obs["provisional_designation"],
                "discovery": df_obs["discovery"].eq("*"),
                "epoch": [t.to_julian_day() for t in dt_objects],
                "epoch_seconds_UTC": float_epochs_utc,
                "epoch_seconds_TDB": float_epochs_tdb,
                "RA": ra_rad,
                "DEC": dec_rad,
                "observatory": df_obs["observatory"],
                "magnitude": df_obs["magnitude_n"],
                "band": df_obs["band"],
                "note1": df_obs["note1"],
                "note2": df_obs["note2"],
                "catalog": None,
                "spacecraft_parallax_type": df_obs["satellite_parallax_type_n"],
                "spacecraft_position_x": df_obs["satellite_x_m"],
                "spacecraft_position_y": df_obs["satellite_y_m"],
                "spacecraft_position_z": df_obs["satellite_z_m"],
            }
        )

    if not radar_table.empty:
        final_df = pd.concat([final_df, radar_table], ignore_index=True, sort=False)
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

        with open(path, "r") as f:
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
    raw_number = row["number"]
    perm_id = str(raw_number).strip() if pd.notna(raw_number) and raw_number else ""

    raw_prov = row["provisional_designation"]
    prov_id = str(raw_prov).strip() if pd.notna(raw_prov) and raw_prov else ""

    result = {"obj_type": "Unknown", "unpacked_name": None, "unpacked_number": None}

    # --- PATH A: PERMANENT ID IS PRESENT ---
    if perm_id:
        if re.match(r"^[JSUND]\d{3}S$", perm_id):
            result["obj_type"] = "Natural Satellite"
            if perm_id[0] in unpackers.PLANET_MAP:
                result["unpacked_name"] = unpackers.unpack_permanent_natural_satellite(perm_id)

        elif re.match(r"^\d{4}[PD]$", perm_id):
            result["obj_type"] = "Comet"
            num_val = int(perm_id[0:4])
            result["unpacked_number"] = str(num_val)
            result["unpacked_name"] = f"{num_val}{perm_id[4]}"

        elif re.match(r"^\d{4}I$", perm_id):
            result["obj_type"] = "Interstellar"
            result["unpacked_name"] = f"{int(perm_id[0:4])}I"

        else:
            result["obj_type"] = "Minor Planet"
            result["unpacked_number"] = unpackers.unpack_permanent_minor_planet(perm_id)
            result["unpacked_name"] = f"({result['unpacked_number']})"

    # --- PATH B: ONLY PROVISIONAL ID IS PRESENT ---
    elif prov_id:
        if prov_id[0:3] in unpackers.SURVEY_MAP or (
            len(prov_id) == 7 and prov_id[6].isalpha() and prov_id[6] not in ["I", "Z"]
        ):
            result["obj_type"] = "Minor Planet"
            result["unpacked_name"] = unpackers.unpack_provisional_minor_planet(prov_id)
        else:
            result["obj_type"] = "Comet/Satellite"
            result["unpacked_name"] = unpackers.unpack_provisional_comet_or_satellite(prov_id)

    else:
        raise ValueError("Observation line does not have permanent nor provisional ID.")

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
        if prov_id[0:3] in unpackers.SURVEY_MAP or (
            len(prov_id) == 7 and prov_id[6].isalpha() and prov_id[6] not in ["I", "Z"]
        ):
            result["obj_type"] = "Minor Planet"
            result["unpacked_name"] = unpackers.unpack_provisional_minor_planet(prov_id)
        else:
            result["obj_type"] = "Comet/Satellite"
            result["unpacked_name"] = unpackers.unpack_provisional_comet_or_satellite(prov_id)

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


import re
import datetime
import numpy as np
from astropy.table import Table
from tudatpy.astro.time_representation import DateTime
from tudatpy.astro import time_representation

# Import the refactored unpacker functions and constants
from . import unpackers


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
