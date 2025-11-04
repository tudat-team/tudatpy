import re
import datetime
import numpy as np
from astropy.table import Table
from tudatpy.astro.time_representation import DateTime

# Import the refactored unpacker functions and constants
from . import unpackers

def parse_packed_permanent_designation(packed_perm_num: str) -> dict:
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
        ident_data['type'] = 'Natural Satellite'
        ident_data['name'] = unpackers.unpack_permanent_natural_satellite(packed_perm_num)

    # Rule for Comets (e.g., 0029P)
    elif re.match(r"^\d{4}[PD]$", packed_perm_num):
        ident_data['type'] = 'Comet'
        ident_data['number'] = str(int(packed_perm_num[0:4]))
        ident_data['comettype'] = packed_perm_num[4]

    # Rule for Interstellar Objects (e.g., 0002I)
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

    return ident_data


def parse_80cols_identification_fields(line: str) -> dict:
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

    if line[12].strip() == '*':
        ident_data['discovery'] = True

    if packed_perm_num:
        ident_data.update(parse_packed_permanent_designation(packed_perm_num))
        if packed_prov_desig:
            if ident_data.get('type') == 'Minor Planet':
                ident_data['desig'] = unpackers.unpack_provisional_minor_planet(packed_prov_desig)
            elif ident_data.get('type') in ['Comet', 'Interstellar Object']:
                ident_data['desig'] = packed_prov_desig.strip() or 'NaN'
    elif packed_prov_desig:
        if packed_prov_desig[6].isalpha() and packed_prov_desig[6] not in ['I', 'Z']:
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
        'epoch', 'RA', 'DEC', 'magnitude', and 'observatory'.
    """
    parsed_observations = []
    with open(filename, 'r') as f:
        for line in f:
            if len(line) < 80:
                continue

            try:
                ident_data = parse_80cols_identification_fields(line)
            except (ValueError, IndexError):
                continue

            try:
                # Filter based on the imported constant
                if line[14].strip() in unpackers.OBS_TYPES_TO_DROP:
                    continue

                # --- Time Parsing ---
                year = int(line[15:19])
                month = int(line[20:22])
                day_frac = float(line[23:32])
                day = int(day_frac)
                remainder = day_frac - day
                total_seconds = remainder * 86400
                hours = int(total_seconds / 3600)
                minutes = int((total_seconds % 3600) / 60)
                seconds = total_seconds % 60
                microseconds = int((seconds - int(seconds)) * 1_000_000)

                obs_time_utc = datetime.datetime(year, month, day, hours, minutes, int(seconds), microseconds)
                dt_object = DateTime.from_python_datetime(obs_time_utc)

                # --- RA/Dec Parsing ---
                ra_parts = line[32:44].strip().split()
                ra_hr, ra_min, ra_sec = int(ra_parts[0]), int(ra_parts[1]), float(ra_parts[2])
                ra_deg = (ra_hr * 15 + ra_min * 0.25 + ra_sec * (15 / 3600))
                ra_rad = np.radians(np.degrees((np.radians(ra_deg) + np.pi) % (2 * np.pi) - np.pi))

                dec_parts = line[44:56].strip().split()
                dec_sign = -1 if dec_parts[0].startswith('-') else 1
                dec_deg_val = int(dec_parts[0].replace('-', '').replace('+', ''))
                dec_min, dec_sec = int(dec_parts[1]), float(dec_parts[2])
                dec_deg = dec_sign * (dec_deg_val + dec_min / 60 + dec_sec / 3600)
                dec_rad = np.radians(dec_deg)

                mag_str = line[65:70].strip()

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
                    "magnitude": float(mag_str) if mag_str else None,
                    "band": line[70].strip() or None,
                    "note1": line[13].strip() or None,
                    "note2": line[14].strip() or None,
                    "catalog": None
                }
                parsed_observations.append(final_data)
            except (ValueError, IndexError):
                continue

    return Table(parsed_observations)