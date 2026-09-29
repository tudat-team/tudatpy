import os
import re

import numpy as np
import pandas as pd
from astropy.table import Table
import astropy.units as u

from tudatpy import constants
from tudatpy.astro import time_representation
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
    for column in df.columns:
        if df.get("sys") is not None:
            if df["sys"].isin(["ICRF_AU", "ICRF_KM"]):
                result_data["observatoryTipe"] = "Spacecraft"
                result_data["x"] = df["pos1"]
                result_data["y"] = df["pos2"]
                result_data["z"] = df["pos3"]
                result_data["vx"] = df["vel1"]
                result_data["vy"] = df["vel2"]
                result_data["vz"] = df["vel3"]
                result_data["ctr"] = df["ctr"]
                result_data["sys"] = df["sys"]
                # the covariance data is optional, if unavailable returns None
                result_data["posCov11"] = df.get("posCov11")
                result_data["posCov12"] = df.get("posCov12")
                result_data["posCov13"] = df.get("posCov13")
                result_data["posCov22"] = df.get("posCov22")
                result_data["posCov23"] = df.get("posCov23")
                result_data["posCov33"] = df.get("posCov33")

            elif df["sys"].isin(["WGS84", "ITRF", "IAU"]):
                result_data["observatoryTipe"] = "Earth roving observatory"
                result_data["x"] = df["pos1"]
                result_data["y"] = df["pos2"]
                result_data["z"] = df["pos3"]
                result_data["vx"] = df["vel1"]
                result_data["vy"] = df["vel2"]
                result_data["vz"] = df["vel3"]
                result_data["ctr"] = df["ctr"]
                result_data["sys"] = df["sys"]
                # the covariance data is optional, if unavailable returns None
                result_data["posCov11"] = df.get("posCov11")
                result_data["posCov12"] = df.get("posCov12")
                result_data["posCov13"] = df.get("posCov13")
                result_data["posCov22"] = df.get("posCov22")
                result_data["posCov23"] = df.get("posCov23")
                result_data["posCov33"] = df.get("posCov33")
        else:
            print("No roving observatories in this data file")

    return result_data


# def _check_possible_failures(df):


def parse_ades_data(file_path: str, format: str):  # -> Table:
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
        Astropy Table with 'number' column already unpacked to human-readable format.
        (e.g., '0003I' -> '3I', '00433' -> '433')

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
        rows = []

        for obs_type in root.iter():
            # optocal observations
            if _tag_name(obs_type) == "optical":
                row = {_tag_name(child): (child.text or "").strip() for child in obs_type}
                rows.append(row)

        df = pd.DataFrame(rows)

        result_data = pd.DataFrame(
            {
                "number": _first_present_column(df, "permID", "provID", "artSat", "trkSub"),
                "provisional_designation": df.get("provID"),
                "epoch": df["obsTime"],
                "epoch_seconds_UTC": None,  # transform obsTime in secs since J2000 UTC. obs_time.astype("int64") / 1_000_000_000,
                "RA": ((df["ra"]).float() * u.deg).to(u.rad).value,
                "DEC": ((df["dec"]).float() * u.deg).to(u.rad).value,
                "observatory": df["stn"],
                "magnitude": df["mag"],
                "band": df["band"],
                "catalog": df["astCat"],
                "mode": df["mode"],
            }
        )

        # add the rmsRA, rmsDec, rmsTime if availabke
        for column in ("rmsRA", "rmsDec", "rmsTime"):
            if column in df.columns:
                result_data[column] = pd.to_numeric(df[column])

        # add the info about roving observatories if needed
        result_data = _roving_obs_columns(df, result_data)

    # elif format == 'PSV':

    else:
        raise ValueError("The given format in invalid; specify either PSV or XML")
