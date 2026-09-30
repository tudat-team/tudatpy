import pandas as pd
from astropy.table import Table
import astropy.units as u
from astroquery.mpc import MPC
import numpy as np

# check these imports - maybe useful for radar??
from tudatpy.data_input.tracking_data.optical_utilities import SPACECRAFT_POSITION_COLUMNS
from tudatpy.data_input.tracking_data.radar_utilities import (
    RADAR_TABLE_META_KEY,
    empty_radar_table,
    radar_data_from_raw,
)
import xml.etree.ElementTree as ET


def _tag_name(element):
    """Return a tag name without an optional XML namespace."""
    return element.tag.rsplit("}", 1)[-1]


def _first_present_column(df, *names):
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
        "spaecraft_velocity_vx": "vel1",
        "spaecraft_velocity_vy": "vel2",
        "spaecraft_velocity_vz": "vel3",
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

    output_columns = ["observatoryTipe", *roving_columns]

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
                        value = (float(value) * u.au).to(u.m)
                    elif input_column in ("vel1", "vel2", "vel3"):
                        value = (float(value) * u.au / u.day).to(u.m / u.s)
                    elif input_column in (
                        "posCov11",
                        "posCov12",
                        "posCov13",
                        "posCov22",
                        "posCov23",
                        "posCov33",
                    ):
                        value = (float(value) * u.au**2).to(u.m**2)

                elif sys == "ICRF_KM":
                    if input_column in ("pos1", "pos2", "pos3"):
                        value = (float(value) * u.km).to(u.m)
                    elif input_column in ("vel1", "vel2", "vel3"):
                        value = (float(value) * u.km / u.s).to(u.m / u.s)
                    elif input_column in (
                        "posCov11",
                        "posCov12",
                        "posCov13",
                        "posCov22",
                        "posCov23",
                        "posCov33",
                    ):
                        value = (float(value) * u.km**2).to(u.m**2)

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
        for column in ("mode", "station", "obsTime", "astCat", "ra", "dec"):
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

        # check that ra, dec, pos1, pos2, pos3 are convertible to floats
        for column in ("ra", "dec", "pos1", "pos2", "pos3"):
            try:
                df[column] = pd.to_numeric(df[column], errors="raise").astype(float)
            except:
                raise ValueError(f"Column '{column}' must contain values convertible to floats")

        # check format of obsTime (yyyy-mm-ddThh:mm:ss.ssZ)
        pattern = r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{2}Z"
        # eliminate trailing or leading balnk spaces
        df["obsTime"] = df["obsTime"].astype("string").str.strip()
        if not df["obsTime"].astype("string").str.fullmatch(pattern).fillna(False).all():
            raise ValueError("Column 'obsTime' must match yyyy-mm-ddThh:mm:ss.ssZ")

        # check that the station is part of the MPC list
        # load the observatories table from the MPC
        observatories_table = MPC.get_observatory_codes().to_pandas()
        stn_codes = observatories_table["code"].astype("string").tolist()
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
        format="%Y-%m-%dT%H:%M:%S.%fZ",
        utc=True,
    )

    j2000 = pd.Timestamp("2000-01-01T12:00:00Z")

    df["epoch_seconds_UTC"] = (epochs - j2000).dt.total_seconds()

    return df


def parse_ades_file(file_path: str, format: str):  # -> Table:
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

    if format == "XML":
        # read the xml elements and define an astropy table
        tree = ET.parse(f"{file_path}")
        root = tree.getroot()

        # separate the observations according to their type: optical, radar, offset, occultation
        # this creates a dictionary with the key = obs_kind and the item = list of rows corresponding to that observation kind
        rows_by_kind = {}
        for obs_type in root.iter():
            obs_kind = _tag_name(obs_type)

            if obs_kind not in ["optical", "radar", "offset", "occultation"]:
                raise ValueError("The observation type in not supported")

            row = {_tag_name(child): (child.text or "").strip() for child in obs_type}

            if obs_kind not in rows_by_kind:
                rows_by_kind[obs_kind] = []

            rows_by_kind[obs_kind].append(row)

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
                        "RA": (((pd.to_numeric(df["ra"])) * u.deg).to(u.rad).value + np.pi)
                        % (2 * np.pi)
                        - np.pi,
                        "DEC": ((pd.to_numeric(df["dec"])) * u.deg).to(u.rad).value,
                        "observatory": df["stn"],
                        "magnitude": df.get(
                            "mag"
                        ),  # this is not a required field so it could be None
                        "band": df.get("band"),  # this is not a required field so it could be None
                        "catalog": df["astCat"],
                        "mode": df["mode"],
                    }
                )

                # add the rmsRA, rmsDec, rmsTime if available
                # for the rows where the info is not available it puts Nan
                for column in ("rmsRA", "rmsDec", "rmsTime"):
                    if column in df.columns:
                        result_data_optical[column] = pd.to_numeric(df[column])

                # add the info about roving observatories if needed - new corresponding columns
                # are created, if the value is missing for some rows then insert None
                result_data_optical = _roving_obs_columns(df, result_data_optical)

                # transform all missing values or empty values to None in the final table
                result_data_optical = result_data_optical.replace(r"^\s*$", None, regex=True)
                result_data_optical = result_data_optical.astype(object).where(pd.notna(df), None)

                # transform into an astropy Table
                optical_table = Table.from_pandas(result_data_optical)

                # this astropy table should be readable by the optical_utilites/read_astropy_optical_data()
                # to be able to create a tracking data object

                return optical_table

    # elif format == 'PSV':

    else:
        raise ValueError("The given format in invalid; specify either PSV or XML")
