import pandas as pd
from astropy.table import Table
import astropy.units as u
from astroquery.mpc import MPC
import numpy as np
from pathlib import Path
import tempfile
import warnings

from tudatpy.data_input.tracking_data.radar_utilities import (
    RADAR_TABLE_META_KEY,
    radar_data_from_raw,
)
import xml.etree.ElementTree as ET

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

# ADES fields with the uncertainties of an observation as reported by the observer (observation fields), and the corresponding
# uncertainties adopted by an orbit computing center (residual fields), which also provide the adopted biases
OBSERVATION_UNCERTAINTY_FIELDS = ("rmsRA", "rmsDec", "rmsCorr", "rmsTime")
RESIDUAL_UNCERTAINTY_FIELDS = ("sigRA", "sigDec", "sigCorr", "sigTime")
RESIDUAL_BIAS_FIELDS = ("biasRA", "biasDec", "biasTime")
RESIDUAL_COLUMNS = [*RESIDUAL_UNCERTAINTY_FIELDS, *RESIDUAL_BIAS_FIELDS]

# output columns that are kept as floats (with NaN for missing values) in the parsed table
NUMERIC_OUTPUT_COLUMNS = [
    "epoch_seconds_UTC",
    "RA",
    "DEC",
    "magnitude",
    "rmsRA",
    "rmsDec",
    "rmsCorr",
    "rmsTime",
    "uncTime",
    *RESIDUAL_COLUMNS,
    "spacecraft_position_x",
    "spacecraft_position_y",
    "spacecraft_position_z",
    "spacecraft_velocity_vx",
    "spacecraft_velocity_vy",
    "spacecraft_velocity_vz",
    "posCov11",
    "posCov12",
    "posCov13",
    "posCov22",
    "posCov23",
    "posCov33",
    "roving_position_1",
    "roving_position_2",
    "roving_position_3",
]

SPACECRAFT_UNITS = {
    "ICRF_AU": (u.au, u.au / u.day, u.au**2),
    "ICRF_KM": (u.km, u.km / u.s, u.km**2),
}

SI_UNITS = (u.m, u.m / u.s, u.m**2)

SPACECRAFT_COLUMNS = {
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

ROVING_COLUMNS = {
    "roving_position_1": "pos1",
    "roving_position_2": "pos2",
    "roving_position_3": "pos3",
    "ctr": "ctr",
    "sys": "sys",
}

POSITION_COLS = ("pos1", "pos2", "pos3")
VELOCITY_COLS = ("vel1", "vel2", "vel3")
COVARIANCE_COLS = ("posCov11", "posCov12", "posCov13", "posCov22", "posCov23", "posCov33")


def _tag_name(element):
    """Return a tag name without an optional XML namespace."""
    return element.tag.rsplit("}", 1)[-1]


def _first_present_column(df, *names):
    if not names:
        raise ValueError("At least one column name must be provided")

    result = None
    for name in names:
        if name in df.columns:
            values = df[name].replace(r"^\s*$", None, regex=True)
            result = values if result is None else result.fillna(values)
    return result


def _unit_for_column(input_column, sys):
    """Return (input_unit, output_unit) for a column, or None if no conversion applies."""
    pos_u, vel_u, cov_u = SPACECRAFT_UNITS[sys]
    si_pos, si_vel, si_cov = SI_UNITS

    if input_column in POSITION_COLS:
        return pos_u, si_pos
    elif input_column in VELOCITY_COLS:
        return vel_u, si_vel
    elif input_column in COVARIANCE_COLS:
        return cov_u, si_cov
    else:
        return None


def _roving_obs_columns(df, result_data):
    # add the corresponding columns if a roving observatory is present
    # the positions and velocities of the spacecraft are returned in [m] and [m/s]
    if result_data is None:
        raise ValueError("The input result table is empty")

    if len(df) != len(result_data):
        raise ValueError("df and result_data must have the same number of rows")

    output_columns = ["observatoryType", *SPACECRAFT_COLUMNS, *ROVING_COLUMNS]

    # if columns not already present in the result_data dataframe add them
    for column in output_columns:
        if column not in result_data.columns:
            result_data[column] = None

    # loop over each row of the df dataframe
    for position in range(len(df)):
        row = df.iloc[position]
        sys = row.get("sys")

        if sys in ("ICRF_AU", "ICRF_KM"):
            observatory_type = "Spacecraft"
        elif sys in ("WGS84", "ITRF", "IAU"):
            observatory_type = "Earth Roving"
        elif pd.isna(sys):
            continue
        else:
            raise ValueError(f"Unsupported sys value for the observer location: {sys}")

        result_data.iat[position, result_data.columns.get_loc("observatoryType")] = observatory_type

        if observatory_type == "Spacecraft":
            for output_column, input_column in SPACECRAFT_COLUMNS.items():
                value = row.get(input_column)
                if pd.notna(value):
                    units = _unit_for_column(input_column, sys)
                    if units is not None:
                        unit_in, unit_out = units
                        value = (float(value) * unit_in).to_value(unit_out)
                elif pd.isna(value):
                    value = None

                result_data.iat[position, result_data.columns.get_loc(output_column)] = value

        elif observatory_type == "Earth Roving":
            for output_column, input_column in ROVING_COLUMNS.items():
                value = row.get(input_column)

                # Store missing pandas values as Python None.
                if pd.isna(value):
                    value = None

                result_data.iat[position, result_data.columns.get_loc(output_column)] = value

    return result_data


def _check_required_columns(df, column_list):
    for column in column_list:
        if column not in df.columns:
            raise ValueError(f"Missing required column: {column}")
        if df[column].isna().any():
            raise ValueError(f"Column '{column}' contains missing values")


def _check_identifier(df, identifier_list):
    if not df.reindex(columns=identifier_list).notna().any(axis=1).all():
        raise ValueError("Each row must have at least one non-missing identifier")


def _check_obstime_format(df):
    pattern = r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{1,6})?Z"
    df["obsTime"] = df["obsTime"].astype("string").str.strip()
    if not df["obsTime"].str.fullmatch(pattern).fillna(False).all():
        raise ValueError("Column 'obsTime' must match yyyy-mm-ddThh:mm:ss[.ssssss]Z")


def _check_stn_validity(df, column_list):
    observatories_table = MPC.get_observatory_codes().to_pandas()
    stn_codes = observatories_table["Code"].astype("string").tolist()
    # extract the stations from the input file to check them
    for el in column_list:
        stn = df[el].astype("string").str.strip()
        invalid = stn.notna() & ~stn.isin(stn_codes)
        if invalid.any():
            raise ValueError(f"Unknown MPC observatory codes: {stn[invalid].unique().tolist()}")


def _check_float_conversion(df, column_list):
    # this only checks float conversion, it doesn't check if the column is present or not
    for column in column_list:
        if column in df.columns:
            values = df[column].dropna()
            if not pd.to_numeric(values, errors="coerce").notna().all():
                raise ValueError(f"Column '{column}' must be convertible to floats")


def _check_location_group(df, column_list, position_columns):
    present = df.reindex(columns=column_list).notna()
    if (present.any(axis=1) & ~present.all(axis=1)).any():
        raise ValueError(
            "If any of the identifier of the location group is provided, all must be provided"
        )

    # check that if present pos1, pos2 and pos3 are convertible to floats
    if (
        all(col in df.columns for col in position_columns)
        and not df.loc[present.all(axis=1), position_columns]
        .apply(pd.to_numeric, errors="coerce")
        .notna()
        .all()
        .all()
    ):
        raise ValueError("The position values must be convertible to floats")


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
        _check_required_columns(df, ["mode", "stn", "obsTime", "astCat", "ra", "dec"])

        # check that at least one identifier is present
        id_columns = ["permID", "provID", "artSat", "trkSub"]
        _check_identifier(df, id_columns)

        # check that if the location group is present that the required fields are there
        columns = ["sys", "ctr", "pos1", "pos2", "pos3"]
        position_columns = ["pos1", "pos2", "pos3"]
        _check_location_group(df, columns, position_columns)

        # check that ra and dec are convertible to floats
        _check_float_conversion(
            df, ["ra", "dec", "rmsRA", "rmsDec", "rmsCorr", "rmsTime", "uncTime", *RESIDUAL_COLUMNS]
        )

        # check that obsTime is in the right format
        _check_obstime_format(df)

        # check that the station is part of the MPC list (a recognized station code)
        _check_stn_validity(df, ["stn"])

        # if all tests are passed then return df
        return df

    elif obs_kind == "radar":
        """
        Required fields include:
        - one between permID, provID or ArtSat
        - trx
        - rcv
        - obsTime
        - doppler OR delay
        - rmsDoppler or rmsDelay
        - frq
        """
        # check that all required columns are present
        _check_required_columns(df, ["trx", "rcv", "obsTime", "frq"])

        # check that at least one identifier is present
        id_columns = ["permID", "provID", "artSat"]
        _check_identifier(df, id_columns)

        # check that at least one between delay and doppler is present in each row
        id_columns = ["doppler", "delay"]
        if not df.reindex(columns=id_columns).notna().any(axis=1).all():
            raise ValueError("Each row must have at least one between doppler and delay")

        # if delay is present, rmsDelay must be present too (same for doppler / rmsDoppler)
        pairs = {"delay": "rmsDelay", "doppler": "rmsDoppler"}
        for value_col, rms_col in pairs.items():
            has_value = df.reindex(columns=[value_col])[value_col].notna()
            has_rms = df.reindex(columns=[rms_col])[rms_col].notna()
            if (has_value & ~has_rms).any():
                raise ValueError(f"Each row with {value_col} must also have {rms_col}")

        # check that if present, delay, rmsDelay, doppler, rmsDoppler, frq are convertible to floats
        _check_float_conversion(df, ["delay", "rmsDelay", "doppler", "rmsDoppler", "frq"])

        # check that obsTime is in the right format
        _check_obstime_format(df)

        # check that both transmitter and receiver stations are part of the MPC list (a recognized station code)
        _check_stn_validity(df, ["rcv", "trx"])

        # if all tests are passed then return df
        return df

    elif obs_kind == "occultation":
        """
        Required fields include:
        - one between permID, provID, artSat, trkSub
        - mode
        - stn
        - obsTime
        - raStar
        - decStar
        - deltaRa, deltaDec OR dist, pa
        - astCat


        If the location group is present then the required elements include:
        - sys
        - ctr
        - pos1
        - pos2
        - pos3
        """
        # check that all required columns are present
        _check_required_columns(df, ["mode", "stn", "obsTime", "astCat", "raStar", "decStar"])

        # check that at least one identifier is present
        id_columns = ["permID", "provID", "artSat", "trkSub"]
        _check_identifier(df, id_columns)

        # check that at least one between the combinations deltaRa, deltaDec and dist, pa is present in each row
        o = (
            df.reindex(columns=["deltaRA", "deltaDec", "dist", "pa"])
            .apply(pd.to_numeric, errors="coerce")
            .notna()
        )
        bad = ~((o["deltaRA"] & o["deltaDec"]) | (o["dist"] & o["pa"]))
        if bad.any():
            raise ValueError(
                "Each row must have at least one between the combinations deltaRa/deltaDec or dist/pa"
            )

        # check that if the location group is present that the required fields are there
        columns = ["sys", "ctr", "pos1", "pos2", "pos3"]
        position_columns = ["pos1", "pos2", "pos3"]
        _check_location_group(df, columns, position_columns)

        # check that raStar and decStar are convertible to floats
        _check_float_conversion(
            df,
            [
                "raStar",
                "decStar",
                "deltaRA",
                "deltaDec",
                "dist",
                "pa",
                "rmsRA",
                "rmsDec",
                "rmsDist",
                "rmsPA",
                "rmsCorr",
                "rmsTime",
                "uncTime",
                *RESIDUAL_COLUMNS,
            ],
        )

        # check that obsTime is in the right format
        _check_obstime_format(df)

        # check that the station is part of the MPC list (a recognized station code)
        _check_stn_validity(df, ["stn"])

        # if all tests are passed then return df
        return df

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
    # converts the mode entry of ADES to the note2 format of the 80-column format
    # this allows compatibility with the algorithm that assigns weights according to VFCC17
    # which uses the value of note2 as reference
    df["note2"] = df["mode"].map({key: value for key, value in MODE_DICT.items()})
    return df


def _astCat_to_catalog(df):
    # converts the astCat entry of ADES to the catalog format of the 80-column format
    # this allows compatibility with the algorithm that assigns weights according to VFCC17
    # which uses the value of catalog as reference
    df["catalog"] = df["astCat"].map({key: value for key, value in CATALOG_DICT.items()})
    return df


def _numeric(df, column, scale=1.0):
    """Numeric column times `scale`; all-NaN if the column is absent."""
    if column not in df.columns:
        return pd.Series(np.nan, index=df.index)
    return pd.to_numeric(df[column], errors="coerce") * scale


def _transform_delta_to_radec(df):
    # from the coordinates of the star and the provided offset define
    # the coordinates in Ra/Dec of the small body
    # return ra and Dec in degrees and save them in the dataframe
    df["ra"] = (
        ((pd.to_numeric(df["raStar"]).to_numpy() * u.deg))
        + (pd.to_numeric(df["deltaRA"]).to_numpy() * u.arcsec).to(u.deg)
        / np.cos((pd.to_numeric(df["decStar"]).to_numpy() * u.deg).to(u.rad))
    ).value
    df["dec"] = (
        ((pd.to_numeric(df["decStar"]).to_numpy() * u.deg))
        + (pd.to_numeric(df["deltaDec"]).to_numpy() * u.arcsec).to(u.deg)
    ).value

    # in this case the uncertainties are already correctly reported in Ra and Dec so there is no need for a transformation

    return df


def _transform_distpa_to_radec(df):
    dist = (
        ((pd.to_numeric(df["dist"]).to_numpy() * u.arcsec)).to(u.deg).value
    )  # transform in degrees
    pa = np.deg2rad(pd.to_numeric(df["pa"]).to_numpy())  # transform in radians
    d_ra = dist * np.sin(pa)
    d_dec = dist * np.cos(pa)
    # return ra and dec in degrees and save them in the dataframe
    df["ra"] = pd.to_numeric(df["raStar"]).to_numpy() + d_ra / np.cos(
        (pd.to_numeric(df["decStar"]).to_numpy() * u.deg).to(u.rad).value
    )
    df["dec"] = pd.to_numeric(df["decStar"]).to_numpy() + d_dec

    # transform the dist and pa uncertainties in uncertainties in RA and Dec
    # identify if the 3 uncertainties columns are present (correlation included)
    if not all(col in df.columns for col in ["rmsDist", "rmsPA"]):
        return df  # no uncertainties reported

    dist = pd.to_numeric(
        df["dist"]
    ).to_numpy()  # here in arcsec to be consistent with units of rmsDist
    rms_dist = pd.to_numeric(df["rmsDist"], errors="coerce").to_numpy()  # in arcsec
    pa = np.deg2rad(pd.to_numeric(df["pa"]).to_numpy())  # transform in radians
    rms_pa = np.deg2rad(pd.to_numeric(df["rmsPA"]).to_numpy())

    if "rmsCorr" in df.columns:
        rho = np.nan_to_num(pd.to_numeric(df["rmsCorr"], errors="coerce").to_numpy(), nan=0.0)
    else:
        rho = np.zeros(len(df))

    n = len(df)
    # build covariance matrix of dist, pa taking correlation into account when present
    C = np.zeros((n, 2, 2))
    C[:, 0, 0] = rms_dist**2
    C[:, 1, 1] = rms_pa**2
    C[:, 0, 1] = C[:, 1, 0] = rho * rms_dist * rms_pa

    # build the jacobian matrix to perform the transformation to the covariance in RA and Dec
    # pa has been previously already transformed in radians
    J = np.zeros((n, 2, 2))
    J[:, 0, 0] = np.sin(pa)
    J[:, 0, 1] = dist * np.cos(pa)
    J[:, 1, 0] = np.cos(pa)
    J[:, 1, 1] = -dist * np.sin(pa)

    # define the covariance matrix in RA and Dec and extract the uncertainties (in arcsec)
    C_radec = J @ C @ np.swapaxes(J, 1, 2)
    var_ra = C_radec[:, 0, 0]
    var_dec = C_radec[:, 1, 1]
    cov_radec = C_radec[:, 0, 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov_radec / np.sqrt(var_ra * var_dec)
    df["rmsCorr"] = np.clip(
        corr, -1.0, 1.0
    )  # define the correlation (in this case is almost always different than 0 because the rotation introduces correlation)
    df["rmsRA"] = np.sqrt(np.clip(var_ra, 0, None))  # on-sky (RA*cos(dec))
    df["rmsDec"] = np.sqrt(np.clip(var_dec, 0, None))

    return df


def _create_optical_dataframe(df):
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
                df.get("mag"), errors="coerce"
            ),  # this is not a required field so it could be None
            "band": df.get("band"),  # this is not a required field so it could be None
            "astCat": df["astCat"],
            "mode": df["mode"],
            "note2": df["note2"],
            "catalog": df["catalog"],
        }
    )

    # add the uncertainties and biases (observation and residual fields) if available
    # for the rows where the info is not available it puts Nan
    for column in ("rmsRA", "rmsDec", "rmsCorr", "rmsTime", "uncTime", *RESIDUAL_COLUMNS):
        if column in df.columns:
            result_data_optical[column] = pd.to_numeric(df[column])

    return result_data_optical


def parse_ades_file(file_path: str):  # -> Table:
    """
    Parses MPC observation data in the ADES format.

    The function uses vectorized Pandas operations for efficiency. The input
    records may contain optical observations, space-based S/s pairs, and radar
    R/r pairs. Radar observations are stored as a canonical pandas table in
    ``table.meta[RADAR_TABLE_META_KEY]``. A blank radar bounce-point field is
    interpreted as centre of mass (``C``), as found in published historical
    records whose corresponding JPL entries use that value. Offset observations,
    and radar observations that are not w.r.t. the centre of mass, are not yet
    supported, and are skipped with a warning.

    Parameters
    ----------
    file_path : str
        The path to the ADES file (in PSV or XML format)

    Returns
    -------
    Table
        Astropy Table

    Raises
    ------
    ValueError
        If the file is empty, if the format is not the expected psv or xml format
    """

    # check that the input file is given as input
    if not file_path:
        raise ValueError("There is no input file")

    # read the xml elements and define an astropy table
    if Path(file_path).suffix == ".psv":
        # transform the PSV file to a temporary XML file using the MPC ades functionalities
        from ades.psvtoxml import psvtoxml

        with tempfile.TemporaryDirectory() as temporary_directory:
            file_path_xml = str(Path(temporary_directory) / f"{Path(file_path).stem}.xml")
            psvtoxml(file_path, file_path_xml)
            root = ET.parse(file_path_xml).getroot()
    else:
        root = ET.parse(f"{file_path}").getroot()

    # Collect observations, including those outside obsData blocks.
    rows_by_kind = {}
    for obs_type in root.iter():
        obs_kind = _tag_name(obs_type)

        if obs_kind not in {"optical", "radar", "offset", "occultation"}:
            continue

        # Include direct child fields, e.g. permID, obsTime, ra, dec.
        row = {_tag_name(child): (child.text or "").strip() for child in obs_type}

        rows_by_kind.setdefault(obs_kind, []).append(row)

    result_data_optical = pd.DataFrame()
    result_data_radar = pd.DataFrame()
    result_data_occultation = pd.DataFrame()

    for obs_kind, rows in rows_by_kind.items():
        if obs_kind not in {"optical", "radar", "occultation"}:
            warnings.warn(f"Skipping unsupported observation kind: {obs_kind}")
            continue

        df = pd.DataFrame(rows)

        # in creating the previous dataframe, if one key is missing from the optical element, its value is assigned
        #  to Nan, None or pd.NA - change all of these and all empty strings to None for consistency
        df = df.replace(r"^\s*$", None, regex=True)
        df = df.astype(object).where(pd.notna(df), None)

        # validate the dataframe, making sure that all required elements are present (according to the observation kind)
        df = _check_possible_failures(df, obs_kind)

        # create the columns containing seconds since J2000 UTC
        df = _epochs_UTC_to_seconds_UTC(df)

        # optical observations
        if obs_kind == "optical":

            # create the 'note2' and 'catalog' columns to allow application of the VFCC17 weighing scheme
            df = _mode_to_note2(df)
            df = _astCat_to_catalog(df)

            result_data_optical = _create_optical_dataframe(df)

            # add the info about roving observatories if needed - new corresponding columns
            # are created, if the value is missing for some rows then insert None
            result_data_optical = _roving_obs_columns(df, result_data_optical)

        elif obs_kind == "radar":
            # only radar observations w.r.t. the centre of mass (com = 1) are currently supported;
            # observations w.r.t. the peak power (com = 0) are skipped
            if "com" in df.columns:
                not_com = df["com"].notna() & (df["com"].astype(str).str.strip() != "1")
                if not_com.any():
                    warnings.warn(
                        f"Skipping {int(not_com.sum())} radar observation(s) that are not "
                        "w.r.t. the centre of mass (com = 0): this is not presently supported"
                    )
                    df = df.loc[~not_com].copy()

            # this is an outdated entry only present in 80 cols file, added here for consistency
            df["target_point"] = "C"

            result_data_radar = pd.DataFrame(
                {
                    "target_body": _first_present_column(df, "permID", "provID", "artSat"),
                    "epoch": df["obsTime"],
                    "epoch_seconds_UTC": df["epoch_seconds_UTC"],
                    "transmitter": df["trx"],
                    "receiver": df["rcv"],
                    "target_point": df["target_point"],
                    "transmitter_frequency_hz": _numeric(df, "frq", 1e6),
                    # in the obs either the delay or the doppler is reported so depending on the case some of these columns might be = None
                    "delay_us": _numeric(df, "delay", 1e6),
                    "delay_sigma_us": _numeric(df, "rmsDelay"),
                    "doppler_hz": _numeric(df, "doppler"),
                    "doppler_sigma_hz": _numeric(df, "rmsDoppler"),
                }
            )

            result_data_radar = radar_data_from_raw(result_data_radar, source="MPC")

        elif obs_kind == "occultation":
            # create the 'note2' and 'catalog' columns to allow application of the VFCC17 weighing scheme
            df = _mode_to_note2(df)
            df = _astCat_to_catalog(df)

            # transform the entries in the occultation table to standard RA, Dec and rmsRA, rmsDec;
            # each row provides either deltaRA/deltaDec or dist/pa, so both sets of rows are
            # transformed separately
            has_delta = df.reindex(columns=["deltaRA", "deltaDec"]).notna().all(axis=1)
            transformed = []
            if has_delta.any():
                transformed.append(_transform_delta_to_radec(df.loc[has_delta].copy()))
            if (~has_delta).any():
                transformed.append(_transform_distpa_to_radec(df.loc[~has_delta].copy()))
            df = pd.concat(transformed).sort_index()

            result_data_occultation = _create_optical_dataframe(df)

            # add the info about roving observatories if needed - new corresponding columns
            # are created, if the value is missing for some rows then insert None
            result_data_occultation = _roving_obs_columns(df, result_data_occultation)

    if result_data_optical.empty and result_data_radar.empty and result_data_occultation.empty:
        # return an empty table
        return Table()

    if not result_data_occultation.empty:
        # append the occultation table to the optical table (they have the same format and are all treated as optical observations)
        result_data_optical = pd.concat(
            [result_data_optical, result_data_occultation], ignore_index=True
        )

    # keep the numeric columns as floats (missing values are NaN), and transform all missing
    # or empty values in the remaining (text) columns to None in the final table
    numeric_columns = result_data_optical.columns.intersection(NUMERIC_OUTPUT_COLUMNS)
    text_columns = result_data_optical.columns.difference(NUMERIC_OUTPUT_COLUMNS)
    result_data_optical[numeric_columns] = result_data_optical[numeric_columns].apply(pd.to_numeric)
    text_data = result_data_optical[text_columns].replace(r"^\s*$", None, regex=True)
    result_data_optical[text_columns] = text_data.astype(object).where(pd.notna(text_data), None)

    optical_table = Table.from_pandas(result_data_optical)

    # if available add radar data as metadata of the optical astropy table
    if not result_data_radar.empty:
        optical_table.meta[RADAR_TABLE_META_KEY] = result_data_radar

    return optical_table
