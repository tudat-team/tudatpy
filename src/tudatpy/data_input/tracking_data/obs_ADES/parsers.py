import pandas as pd
from astropy.table import Table
import astropy.units as u
from astroquery.mpc import MPC
import numpy as np
import re

# check these imports - maybe useful for radar??
from tudatpy.data_input.tracking_data.optical_utilities import SPACECRAFT_POSITION_COLUMNS
from tudatpy.data_input.tracking_data.radar_utilities import (
    RADAR_TABLE_META_KEY,
    empty_radar_table,
    radar_data_from_raw,
)
import xml.etree.ElementTree as ET

from ades.psvtoxml import psvtoxml

# conversion between note2 and mode according to https://github.com/IAU-ADES/ADES-Master/blob/master/Python/ades/packUtil.py
MODE_DICT = {
    "PHo": "P",
    "PHO": " ",
    "ENC": "e",
    "CCD": "C",
    "CMO": "B",
    "MER": "T",
    "MIC": "M",
    "ccd": "c",
    "OCC": "E",
    "OFF": "O",
    "PMT": "H",
    "NOR": "N",
    "VID": "n",
}

CATALOG_DICT = {
    "UNK": " ",
    "USNOA1": "a",
    "USNOSA1": "b",
    "USNOA2": "c",
    "USNOSA2": "d",
    "UCAC1": "e",
    "Tyc1": "f",
    "Tyc2": "g",
    "GSC1.0": "h",
    "GSC1.1": "i",
    "GSC1.2": "j",
    "GSC2.2": "k",
    "ACT": "l",
    "GSCACT": "m",
    "SDSS8": "n",
    "USNOB1": "o",
    "PPM": "p",
    "UCAC4": "q",
    "UCAC2": "r",
    "USNOB2": "s",
    "PPMXL": "t",
    "UCAC3": "u",
    "NOMAD": "v",
    "CMC14": "w",
    "Hip2": "x",
    "Hip1": "y",
    "GSC": "z",
    "AC": "A",
    "SAO1984": "B",
    "SAO": "C",
    "AGK3": "D",
    "FK4": "E",
    "ACRS": "F",
    "LickGas": "G",
    "Ida93": "H",
    "Perth70": "I",
    "COSMOS": "J",
    "Yale": "K",
    "2MASS": "L",
    "GSC2.3": "M",
    "SDSS7": "N",
    "SSTRC1": "O",
    "MPOSC3": "P",
    "CMC15": "Q",
    "SSTRC4": "R",
    "URAT1": "S",
    "URAT2": "T",
    "Gaia1": "U",
    "Gaia2": "V",
    "Gaia3": "W",
    "Gaia3E": "X",
    "UCAC5": "Y",
    "ATLAS2": "Z",
    "IHW": "0",
    "PS1_DR1": "1",
    "PS1_DR2": "2",
    "Gaia_Int": "3",
    "GZ": "4",
    "UBSC": "9",
    "Gaia_2016": "6",
    "ZZCAT": "7",
    "APASS": "8",
    "AKARI": "!",
    "AG": "@",
    "WISE": "#",
    "LSST2502": "$",
    "IRASPSC": "^",
}


def _tag_name(element):
    """Return a tag name without an optional XML namespace."""
    return element.tag.rsplit("}", 1)[-1]


def _first_present_column(df, *names):
    if not names:
        raise ValueError("At least one column name must be provided")

    for name in names:
        if name in df.columns:
            return df[name]
    return None


def _roving_obs_columns(df, result_data):
    # add the corresponding columns if a roving observatory is present
    # the positions and velocities of the spacecraft are returned in [m] and [m/s]
    spacecraft_columns = {
        "spacecraft_position_x": "pos1",
        "spacecraft_position_y": "pos2",
        "spacecraft_position_z": "pos3",
        "spacecraft_velocity_vx": "vel1",
        "spacecraft_velocity_vy": "vel2",
        "spacecraft_velocity_vz": "vel3",
        "ctr": "ctr",
        "sys": "sys",
        "posCov11": "posCov11",
        "posCov12": "posCov12",
        "posCov13": "posCov13",
        "posCov22": "posCov22",
        "posCov23": "posCov23",
        "posCov33": "posCov33",
    }

    roving_columns = {
        "roving_position_1": "pos1",
        "roving_position_2": "pos2",
        "roving_position_3": "pos3",
        "ctr": "ctr",
        "sys": "sys",
    }

    if result_data is None:
        raise ValueError("The input result table is empty")

    if len(df) != len(result_data):
        raise ValueError("df and result_data must have the same number of rows")

    output_columns = ["observatoryTipe", *spacecraft_columns, *roving_columns]

    # if columns not already presenty in the result_data dataframe add them
    for column in output_columns:
        if column not in result_data.columns:
            result_data[column] = None

    found_roving = False

    # loop over each row of the df dataframe
    for position in range(len(df)):
        row = df.iloc[position]
        sys = row.get("sys")

        if sys in ("ICRF_AU", "ICRF_KM"):
            observatory_type = "Spacecraft"
        elif sys in ("WGS84", "ITRF", "IAU"):
            observatory_type = "Earth Roving"
        else:
            continue

        found_roving = True
        result_data.iat[position, result_data.columns.get_loc("observatoryTipe")] = observatory_type

        # insert in the corresponding row all avalable info about the roving observatory;
        # if info is unavailable the value is None
        if observatory_type == "Spacecraft":
            for output_column, input_column in spacecraft_columns.items():
                value = row.get(input_column)

                if sys == "ICRF_AU":
                    if input_column in ("pos1", "pos2", "pos3"):
                        value = (float(value) * u.au).to_value(u.m)
                    elif input_column in ("vel1", "vel2", "vel3") and pd.notna(value):
                        value = (float(value) * u.au / u.day).to_value(u.m / u.s)
                    elif input_column in (
                        "posCov11",
                        "posCov12",
                        "posCov13",
                        "posCov22",
                        "posCov23",
                        "posCov33",
                    ) and pd.notna(value):
                        value = (float(value) * u.au**2).to_value(u.m**2)

                elif sys == "ICRF_KM":
                    if input_column in ("pos1", "pos2", "pos3"):
                        value = (float(value) * u.km).to_value(u.m)
                    elif input_column in ("vel1", "vel2", "vel3") and pd.notna(value):
                        value = (float(value) * u.km / u.s).to_value(u.m / u.s)
                    elif input_column in (
                        "posCov11",
                        "posCov12",
                        "posCov13",
                        "posCov22",
                        "posCov23",
                        "posCov33",
                    ) and pd.notna(value):
                        value = (float(value) * u.km**2).to_value(u.m**2)

                else:
                    raise ValueError("Unsupported reference frame definition")

                # Store missing pandas values as Python None.
                if pd.isna(value):
                    value = None

                result_data.iat[position, result_data.columns.get_loc(output_column)] = value

        elif observatory_type == "Earth Roving":
            for output_column, input_column in roving_columns.items():
                value = row.get(input_column)

                # Store missing pandas values as Python None.
                if pd.isna(value):
                    value = None

                result_data.iat[position, result_data.columns.get_loc(output_column)] = value

    if not found_roving:
        print("No roving observatories in this data file")

    return result_data


def _check_possible_failures(df, obs_kind):
    # check that all required fields are present and of the expected type
    # does not deal with non-required fields

    if obs_kind == "optical":
        """
        Required fields include:
        - mode
        - station
        - obsTime
        - astCat
        - one between permID, provID, artSat, trkSub
        - ra
        - dec

        If the location group is present then the required elements include:
        - sys
        - ctr
        - pos1
        - pos2
        - pos3
        """
        # check that all required columns are present
        for column in ("mode", "stn", "obsTime", "astCat", "ra", "dec"):
            if column not in df.columns:
                raise ValueError(f"Missing required column: {column}")
            if df[column].isna().any():
                raise ValueError(f"Column '{column}' contains missing values")

        # check that at least one identifier is present
        id_columns = ["permID", "provID", "artSat", "trkSub"]
        if not df.reindex(columns=id_columns).notna().any(axis=1).all():
            raise ValueError("Each row must have at least one non-missing identifier")

        # check that if the location group is present that the required fields are there
        columns = ["sys", "ctr", "pos1", "pos2", "pos3"]
        present = df.reindex(columns=columns).notna()
        if (present.any(axis=1) & ~present.all(axis=1)).any():
            raise ValueError(
                "If any of sys, ctr, pos1, pos2, pos3 is provided, all must be provided"
            )

        # check that if present pos1, pos2 and pos3 are convertible to floats
        position_columns = ["pos1", "pos2", "pos3"]
        if (
            all(col in df.columns for col in position_columns)
            and not df.loc[present.all(axis=1), position_columns]
            .apply(pd.to_numeric, errors="coerce")
            .notna()
            .all()
            .all()
        ):
            raise ValueError("pos1, pos2, and pos3 must be convertible to floats")

        # check that ra and dec are convertible to floats
        for column in ("ra", "dec"):
            try:
                df[column] = pd.to_numeric(df[column]).astype(float)
            except:
                raise ValueError(f"Column '{column}' must contain values convertible to floats")

        # check that if presnt, rmsRA, rmsDec, rmsCorr and rmsTime are also convertible to floats
        rms_columns = ["rmsRA", "rmsDec", "rmsCorr", "rmsTime"]
        for column in rms_columns:
            if column in df.columns:
                values = df[column].dropna()
                if not pd.to_numeric(values, errors="coerce").notna().all():
                    raise ValueError(f"Column '{column}' must be convertible to floats")

        # check that obsTime is in the right format
        pattern = r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{1,4})?Z"

        df["obsTime"] = df["obsTime"].astype("string").str.strip()
        if not df["obsTime"].str.fullmatch(pattern).fillna(False).all():
            raise ValueError("Column 'obsTime' must match yyyy-mm-ddThh:mm:ss[.ssss]Z")

        # check that the station is part of the MPC list (a recognized station code)
        observatories_table = MPC.get_observatory_codes().to_pandas()
        stn_codes = observatories_table["Code"].astype("string").tolist()
        # extract the stations from the input file to check them
        stn = df["stn"].astype("string").str.strip()
        invalid = stn.notna() & ~stn.isin(stn_codes)
        if invalid.any():
            raise ValueError(f"Unknown MPC observatory codes: {stn[invalid].unique().tolist()}")

        # if all tests are passed then return df
        return df

    # elif obs_kind == "radar":

    # elif obs_kind == "offset":

    # elif obs_kind == "occultation":

    else:
        raise ValueError("The input observation kind is invalid")


def _epochs_UTC_to_seconds_UTC(df):
    # transform the given obsTime in seconds since J2000 UTC to be able to
    # provide this to the tracking data object
    epochs = pd.to_datetime(
        df["obsTime"],
        format="ISO8601",
        utc=True,
    )

    j2000 = pd.Timestamp("2000-01-01T12:00:00Z")

    df["epoch_seconds_UTC"] = (epochs - j2000).dt.total_seconds()

    return df


def _mode_to_note2(df):
    # converts the mode entry of ADES to the note2 format of the 80m columns format
    # this allows compatibility with the algorithm that assigns weights according to VFCC17
    # which uses the value of note2 as reference
    df["note2"] = df["mode"].map({key: value for key, value in MODE_DICT.items()})
    return df


def _astCat_to_catalog(df):
    # converts the astCat entry of ADES to the catalog format of the 80m columns format
    # this allows compatibility with the algorithm that assigns weights according to VFCC17
    # which uses the value of catalog as reference
    df["catalog"] = df["astCat"].map({key: value for key, value in CATALOG_DICT.items()})
    return df


def parse_ades_file(file_path: str):  # -> Table:
    """
    Parses MPC observation data in the ADES format.

    The function uses vectorized Pandas operations for efficiency. The input
    records may contain optical observations, space-based S/s pairs, and radar
    R/r pairs. Radar observations are stored as a canonical pandas table in
    ``table.meta[RADAR_TABLE_META_KEY]``. A blank radar bounce-point field is
    interpreted as centre of mass (``C``), as found in published historical
    records whose corresponding JPL entries use that value.

    Parameters
    ----------
    file_path : str
        The path to the ADES files

    format : str
        Either PSV or XML, indicates the format of the ADES file used as input

    Returns
    -------
    Table
        Astropy Table

    Raises
    ------
    ValueError
        If the file is empty, if the format is not the expected psv or xml format
    """

    # 1. LOAD RAW LINES & CHECK LENGTH
    if not file_path:
        raise ValueError("There is no input file")

    if file_path.split(".")[1] == "psv":
        # transform the PSV file on XML format and proceed with the XML parser
        file_path_xml = file_path.split(".")[0] + ".xml"
        # transfrom the psv file to an xml file using the MPC ades functionalities
        psvtoxml(file_path, file_path_xml)
        file_path = file_path_xml

    # read the xml elements and define an astropy table
    tree = ET.parse(f"{file_path}")
    root = tree.getroot()

    # separate the observations according to their type: optical, radar, offset, occultation
    # this creates a dictionary with the key = obs_kind and the item = list of rows corresponding to that observation kind
    rows_by_kind = {}
    for obs_data in root.findall(".//obsData"):
        for obs_type in obs_data:
            obs_kind = _tag_name(obs_type)

            if obs_kind not in {"optical", "radar", "offset", "occultation"}:
                raise ValueError(f"Unsupported observation type: {obs_kind}")

            # Include direct child fields, e.g. permID, obsTime, ra, dec.
            row = {_tag_name(child): (child.text or "").strip() for child in obs_type}

            rows_by_kind.setdefault(obs_kind, []).append(row)

    for obs_kind, rows in rows_by_kind.items():

        df = pd.DataFrame(rows)

        # in creating the previsous dataframe, if one key is missing from the optical element, its value is assigned
        #  to Nan, None or pd.NA - change all of these and all empty strings to None for consistency
        df = df.replace(r"^\s*$", None, regex=True)
        df = df.astype(object).where(pd.notna(df), None)

        # validate  the dataframe, making sure that all required elements are present (according to the observation kind)
        df = _check_possible_failures(df, obs_kind)

        # create the columns containing seconds since J2000 UTC
        df = _epochs_UTC_to_seconds_UTC(df)

        # create the 'note2' and 'catalog' columns to allow application of the VCFF17 weighing scheme
        df = _mode_to_note2(df)
        df = _astCat_to_catalog(df)

        # optical observations
        if obs_kind == "optical":
            result_data_optical = pd.DataFrame(
                {
                    "number": _first_present_column(df, "permID", "provID", "artSat", "trkSub"),
                    "provisional_designation": df.get(
                        "provID"
                    ),  # this is not a required field so it could be None
                    "epoch": df["obsTime"],
                    "epoch_seconds_UTC": df["epoch_seconds_UTC"],
                    "RA": ((pd.to_numeric(df["ra"]).to_numpy() * u.deg).to(u.rad).value + np.pi)
                    % (2 * np.pi)
                    - np.pi,
                    "DEC": (pd.to_numeric(df["dec"]).to_numpy() * u.deg).to(u.rad).value,
                    "observatory": df["stn"],
                    "magnitude": pd.to_numeric(
                        df["mag"], errors="coerce"
                    ),  # this is not a required field so it could be None
                    "band": df.get("band"),  # this is not a required field so it could be None
                    "astCat": df["astCat"],
                    "mode": df["mode"],
                    "note2": df["note2"],
                    "catalog": df["catalog"],
                }
            )

            # add the rmsRA, rmsDec, rmsTime if available
            # for the rows where the info is not available it puts Nan
            for column in ("rmsRA", "rmsDec", "rmsCorr", "rmsTime"):
                if column in df.columns:
                    result_data_optical[column] = pd.to_numeric(df[column])

            # add the info about roving observatories if needed - new corresponding columns
            # are created, if the value is missing for some rows then insert None
            result_data_optical = _roving_obs_columns(df, result_data_optical)

            # transform all missing values or empty values to None in the final table
            result_data_optical = result_data_optical.replace(r"^\s*$", None, regex=True)
            result_data_optical = result_data_optical.astype(object).where(
                pd.notna(result_data_optical), None
            )

            # transform into an astropy Table
            optical_table = Table.from_pandas(result_data_optical)

            # this astropy table should be readable by the optical_utilites/read_astropy_optical_data()
            # to be able to create a tracking data object

            return optical_table
