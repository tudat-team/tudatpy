import re
import numpy as np
from astropy.table import Table
from datetime import timedelta, datetime
from tudatpy.astro.time_representation import DateTime

# Import the refactored unpacker functions and constants
from . import unpackers

def parse_packed_permanent_designation(packed_perm_num: str) -> dict[str,str]:
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

    if not packed_perm_num:
        raise ValueError("Permanent designation string is empty.")

    # Rule for Natural Satellites (e.g., J013S)
    if re.match(r"^[JSUND]\d{3}S$", packed_perm_num):
        ident_data['type'] = 'Natural Satellite'
        ident_data['name'] = unpackers.unpack_permanent_natural_satellite(packed_perm_num)

    # Rule for Comets (e.g., 0029P)
    elif re.match(r"^\d{4}[PD]$", packed_perm_num):
        ident_data['type'] = 'Comet'
        ident_data['number'] = str(int(packed_perm_num[0:4]))
        ident_data['comettype'] = packed_perm_num[4]

    # Rule for Interstellar Objects
    elif re.match(r"^\d{4}I$", packed_perm_num):
        ident_data['type'] = 'Interstellar Object'
        number = int(packed_perm_num[0:4])
        ident_data['name'] = f"{number}I"
        ident_data['number'] = str(number)
        ident_data['comettype'] = 'I'

    # Rule for Minor Planets (e.g., 00433, A0345, D4341, ~000z)
    elif packed_perm_num:
        ident_data['type'] = 'Minor Planet'
        ident_data['number'] = unpackers.unpack_permanent_minor_planet(packed_perm_num)

    else:
        print('Type not recognized')

    return ident_data


def parse_80cols_identification_fields(line: str) -> dict[str,str]:
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

    if not packed_perm_num and not packed_prov_desig:
        raise ValueError(f"Line: {line} contains no Permanent or Provisional identification. "
                         f"Please make sure each line of your file follows the MPC format: "
                         f"https://minorplanetcenter.net/iau/info/OpticalObs.html")

    # Discovery flag
    if line[12] == '*':
        ident_data['discovery'] = True

    if packed_perm_num:
        ident_data.update(parse_packed_permanent_designation(packed_perm_num))
        if packed_prov_desig:
            if ident_data.get('type') == 'Minor Planet':
                ident_data['desig'] = unpackers.unpack_provisional_minor_planet(packed_prov_desig)
            elif ident_data.get('type') in ['Comet', 'Interstellar Object']:
                ident_data['desig'] = packed_prov_desig.strip() or 'NaN'
    elif packed_prov_desig:
        if len(packed_prov_desig) == 7 and packed_prov_desig[6].isalpha() and packed_prov_desig[6] not in ['I', 'Z']:
            ident_data['type'] = 'Minor Planet'
            ident_data['desig'] = unpackers.unpack_provisional_minor_planet(packed_prov_desig)
        else:
            ident_data['type'] = 'Comet'
            ident_data['desig'] = unpackers.unpack_provisional_comet_or_satellite(packed_prov_desig)

    return ident_data


def parse_80cols_file(filename: str) -> Table:
    """
    Reads an MPC observation file (ASCII 80-column format) and parses each line,
    returning an Astropy Table with parsed observations. Lines that are shorter
    than 80 characters or that cause parsing errors are skipped.

    Parameters
    ----------
    filename : str
        The path to the 80-column MPC observation file.

    Returns
    -------
    Table
        An Astropy Table containing the parsed observations. Each row represents
        an observation with columns like 'number', 'provisional_designation',
        'epoch', 'RA', 'DEC', 'magnitude', and 'observatory.
    """
    parsed_observations = []
    skipped_lines = {'length': 0, 'obs_type': 0, 'date': 0, 'coordinate': 0, 'other': 0}

    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, start=1):
            line = line.rstrip('\n')
            if not line:
                continue

            # Validate line length
            if len(line) != 80:
                print(f'Warning: Line {line_num}: Wrong length ({len(line)} != 80) - skipping')
                skipped_lines['length'] += 1
                continue

            try:
                # Parse identification fields
                ident_data = parse_80cols_identification_fields(line)

                # Check observation type filter
                note_2 = line[14].strip()
                if note_2 in unpackers.OBS_TYPES_TO_DROP:
                    reason = unpackers.map_reason_to_drop.get(note_2, "Flag listed in drop list")
                    print(f"Line {line_num}: Skipped - {reason} (Flag: '{note_2}')")
                    skipped_lines['obs_type'] += 1
                    continue

                # Parse time fields
                try:
                    year = int(line[15:19])
                    month = int(line[20:22])
                    day_str = line[23:32]
                    full_day_float = float(day_str)
                    day_int = int(full_day_float)
                    day_frac = full_day_float - day_int

                    base_dt = datetime(year, month, day_int)
                    obs_time_utc = base_dt + timedelta(days=day_frac)
                    dt_object = DateTime.from_python_datetime(obs_time_utc)

                except (ValueError, OverflowError) as e:
                    print(f"Line {line_num}: Invalid date '{line[15:32]}' - {e} - skipping")
                    skipped_lines['date'] += 1
                    continue

                # Parse coordinates and magnitude
                try:
                    ra_h_str, ra_m_str, ra_s_str = line[32:34], line[35:37], line[38:44]
                    dec_sign, dec_d_str, dec_m_str, dec_s_str = line[44], line[45:47], line[48:50], line[51:56]
                    mag_str = line[65:70].strip()

                    ra_deg = (int(ra_h_str) + int(ra_m_str)/60.0 + float(ra_s_str)/3600.0) * 15.0
                    ra_rad = np.deg2rad(ra_deg)

                    dec_deg_abs = int(dec_d_str) + int(dec_m_str)/60.0 + float(dec_s_str)/3600.0
                    dec_deg = -dec_deg_abs if dec_sign == '-' else dec_deg_abs
                    dec_rad = np.deg2rad(dec_deg)

                    mag_val = float(mag_str) if mag_str else None

                except (ValueError, IndexError) as e:
                    print(f"Line {line_num}: Invalid coordinate/magnitude data - {e} - skipping")
                    skipped_lines['coordinate'] += 1
                    continue

                # Build final data dictionary
                final_data = {
                    "number": ident_data.get('number'),
                    "provisional_designation": ident_data.get('desig'),
                    "discovery": ident_data.get('discovery', False),
                    "epoch": dt_object.to_julian_day(),
                    "epoch_utc": obs_time_utc,
                    "epoch_seconds": dt_object.to_epoch(),
                    "RA": ra_rad,
                    "DEC": dec_rad,
                    "observatory": line[77:80].strip(),
                    "magnitude": mag_val,
                    "band": line[70].strip() or None,
                    "note1": line[13].strip() or None,
                    "note2": line[14].strip() or None,
                    "catalog": None
                }

                parsed_observations.append(final_data)

            except Exception as e:
                # Catch-all for unexpected errors
                print(f"Line {line_num}: Unexpected error - {type(e).__name__}: {e} - skipping")
                skipped_lines['other'] += 1
                continue

    # Summary if successful
    total_skipped = sum(skipped_lines.values())
    if total_skipped > 0:
        print(f"\nParsing summary: {len(parsed_observations)} lines successfully parsed, {total_skipped} lines skipped")
        print(f"  - Length errors: {skipped_lines['length']}")
        print(f"  - Observation type filtered: {skipped_lines['obs_type']}")
        print(f"  - Date errors: {skipped_lines['date']}")
        print(f"  - Coordinate errors: {skipped_lines['coordinate']}")
        print(f"  - Other errors: {skipped_lines['other']}")

    return Table(parsed_observations)