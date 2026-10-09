"""Utilities for converting optical astrometry tables to tracking data.

The public readers in this module accept already-loaded pandas or astropy
tables. They share the same optical-data conversion path used by the MPC and
MPC 80-column readers.
"""

import pandas as pd
import numpy as np
import astropy
from astropy_healpix import HEALPix
from astropy.units import Quantity
import astropy.units as u
import os
import re
from tudatpy.astro import time_representation
from tudatpy.data_input.tracking_data import (
    TrackingData,
    TrackingSupplementaryData,
    TranslationalStateSupplementaryData,
)

BIAS_LOWRES_FILE = os.path.join(
    os.path.expanduser("~"),
    ".tudat",
    "resource",
    "star_catalog_biases",
    "debias_2018",
    "bias.dat",
)

# Described here:
# https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
DEFAULT_CATALOG_FLAGS = [
    "a",
    "b",
    "c",
    "d",
    "e",
    "g",
    "i",
    "j",
    "l",
    "m",
    "n",
    "o",
    "p",
    "q",
    "r",
    "t",
    "u",
    "v",
    "w",
    "L",
    "N",
    "Q",
    "R",
    "S",
    "U",
    "Y",
]

REQUIRED_OPTICAL_COLUMNS = ["number", "epoch", "RA", "DEC", "observatory"]
SPACECRAFT_POSITION_COLUMNS = [
    "spacecraft_position_x",
    "spacecraft_position_y",
    "spacecraft_position_z",
]
SPACECRAFT_VELOCITY_COLUMNS = [
    "spacecraft_velocity_vx",
    "spacecraft_velocity_vy",
    "spacecraft_velocity_vz",
]
ROVING_POSITION_COLUMNS = ["roving_position_1", "roving_position_2", "roving_position_3"]
ANCILLARY_STRING_COLUMNS = [
    "band",
    "phottype",
    "catalog",
    "note2",
    "custom_name",
    "mag",
    "discovery",
]

NAIF_IDS = {
    0: "Ssb",
    1: "Mercury Barycenter",
    2: "Venus Barycenter",
    3: "Earth Barycenter",
    4: "Mars Barycenter",
    5: "Jupiter Barycenter",
    6: "Saturn Barycenter",
    7: "Uranus Barycenter",
    8: "Neptune Barycenter",
    9: "Pluto Barycenter",
    10: "Sun",
    199: "Mercury",
    299: "Venus",
    399: "Earth",
    301: "Moon",
    499: "Mars",
    401: "Phobos",
    402: "Deimos",
    599: "Jupiter",
    501: "Io",
    502: "Europa",
    503: "Ganymede",
    504: "Callisto",
    505: "Amalthea",
    506: "Himalia",
    507: "Elara",
    508: "Pasiphae",
    509: "Sinope",
    510: "Lysithea",
    511: "Carme",
    512: "Ananke",
    513: "Leda",
    514: "Thebe",
    515: "Adrastea",
    516: "Metis",
    517: "Callirrhoe",
    518: "Themisto",
    519: "Megaclite",
    520: "Taygete",
    521: "Chaldene",
    522: "Harpalyke",
    523: "Kalyke",
    524: "Iocaste",
    525: "Erinome",
    526: "Isonoe",
    527: "Praxidike",
    528: "Autonoe",
    529: "Thyone",
    530: "Hermippe",
    531: "Aitne",
    532: "Eurydome",
    533: "Euanthe",
    534: "Euporie",
    535: "Orthosie",
    536: "Sponde",
    537: "Kale",
    538: "P asithee",
    539: "Hegemone",
    540: "Mneme",
    541: "Aoede",
    542: "Thelxinoe",
    543: "Arche",
    544: "Kallichore",
    545: "Helike",
    546: "Carpo",
    547: "Eukelade",
    548: "Cyllene",
    549: "Kore",
    550: "Herse",
    553: "Dia",
    699: "Saturn",
    601: "Mimas",
    602: "Enceladus",
    603: "Tethys",
    604: "Dione",
    605: "Rhea",
    606: "Titan",
    607: "Hyperion",
    608: "Iapetus",
    609: "Phoebe",
    610: "Janus",
    611: "Epimetheus",
    612: "Helene",
    613: "Telesto",
    614: "Calypso",
    615: "Atlas",
    616: "Prometheus",
    617: "Pandora",
    618: "Pan",
    619: "Ymir",
    620: "Paaliaq",
    621: "Tarvos",
    622: "Ijiraq",
    623: "Suttungr",
    624: "Kiviuq",
    625: "Mundilfari",
    626: "Albiorix",
    627: "Skathi",
    628: "Erriapus",
    629: "Siarnaq",
    630: "Thrymr",
    631: "Narvi",
    632: "Methone",
    633: "Pallene",
    634: "Polydeuces",
    635: "Daphnis",
    636: "Aegir",
    637: "Bebhionn",
    638: "Bergelmir",
    639: "Bestla",
    640: "Farbauti",
    641: "Fenrir",
    642: "Fornjot",
    643: "Hati",
    644: "Hyrrokkin",
    645: "Kari",
    646: "Loge",
    647: "Skoll",
    648: "Surtur",
    649: "Anthe",
    650: "Jarnsaxa",
    651: "Greip",
    652: "Tarqeq",
    653: "Aegaeon",
    799: "Uranus",
    701: "Ariel",
    702: "Umbriel",
    703: "Titania",
    704: "Oberon",
    705: "Miranda",
    706: "Cordelia",
    707: "Ophelia",
    708: "Bianca",
    709: "Cressida",
    710: "Desdemona",
    711: "Juliet",
    712: "Portia",
    713: "Rosalind",
    714: "Belinda",
    715: "Puck",
    716: "Caliban",
    717: "Sycorax",
    718: "Prospero",
    719: "Setebos",
    720: "Stephano",
    721: "Trinculo",
    722: "Francisco",
    723: "Margaret",
    724: "Ferdinand",
    725: "Perdita",
    726: "Mab",
    727: "Cupid",
    899: "Neptune",
    801: "Triton",
    802: "Nereid",
    803: "Naiad",
    804: "Thalassa",
    805: "Despina",
    806: "Galatea",
    807: "Larissa",
    808: "Proteus",
    809: "Halimede",
    810: "Psamathe",
    811: "Sao",
    812: "Laomedeia",
    813: "Neso",
    999: "Pluto",
    901: "Charon",
    902: "Nix",
    903: "Hydra",
    904: "Kerberos",
    905: "Styx",
}


def standardize_optical_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure IDs are strings and observatory codes are zero-padded.

    This is a supporting utility for inspecting or preparing optical tables.
    In the typical Tudat workflow, call :func:`read_optical_data` or a
    source-specific optical reader instead.

    Parameters
    ----------
    df : pandas.DataFrame
        Optical astrometry table to standardize.

    Returns
    -------
    pandas.DataFrame
        Standardized copy of the input table.
    """
    df = df.copy()

    if "observatory" in df.columns:
        df["observatory"] = df["observatory"].astype(str).str.strip().str.zfill(3)

    if "number" in df.columns:
        df["number"] = df["number"].astype(str).str.strip()

    return df


def validate_optical_table(
    table: astropy.table.QTable | astropy.table.Table | pd.DataFrame, frame: str
) -> None:
    """Validate the common optical astrometry table input.

    This is a supporting utility used by the optical readers. In the typical
    Tudat workflow, call :func:`read_optical_data` or a source-specific optical
    reader instead.

    Parameters
    ----------
    table : astropy.table.QTable | astropy.table.Table | pandas.DataFrame
        Optical astrometry table to validate.
    frame : str
        Reference frame of the input observations.

    Returns
    -------
    None
        The function returns ``None`` if validation succeeds.
    """
    if frame != "J2000":
        raise NotImplementedError("Only observations in J2000 are supported currently")

    if isinstance(table, (astropy.table.QTable, astropy.table.Table)):
        colnames = table.colnames
        nrows = len(table)
    elif isinstance(table, pd.DataFrame):
        colnames = table.columns
        nrows = len(table)
    else:
        raise TypeError(f"Unsupported table type: {type(table).__name__}")

    if not set(REQUIRED_OPTICAL_COLUMNS).issubset(set(colnames)):
        raise ValueError(
            f"Table must include a set of mandatory columns: {REQUIRED_OPTICAL_COLUMNS}"
        )

    if nrows == 0:
        raise ValueError("Table contains zero rows: no valid observations were parsed.")


def create_augmented_optical_table(
    table: astropy.table.QTable | astropy.table.Table | pd.DataFrame,
    in_degrees: bool = True,
    frame: str = "J2000",
    custom_name: str | None = None,
) -> pd.DataFrame:
    """Create the shared augmented optical table used by all optical sources.

    This is a supporting utility for inspecting or preparing optical tables.
    In the typical Tudat workflow, call :func:`read_optical_data` or a
    source-specific optical reader instead.

    Parameters
    ----------
    table : astropy.table.QTable | astropy.table.Table | pandas.DataFrame
        Optical astrometry table to augment.
    in_degrees : bool, default True
        Whether right ascension and declination are provided in degrees.
    frame : str, default "J2000"
        Reference frame of the input observations.
    custom_name : str | None, default None
        Optional target name assigned to all observations.

    Returns
    -------
    pandas.DataFrame
        Augmented optical astrometry table.
    """
    if isinstance(table, (astropy.table.QTable, astropy.table.Table)):
        augmented_table = table.to_pandas()
    else:
        augmented_table = table.copy()

    augmented_table = standardize_optical_dataframe(augmented_table)
    validate_optical_table(augmented_table, frame)

    if custom_name is not None or "custom_name" not in augmented_table.columns:
        augmented_table["custom_name"] = [custom_name] * len(augmented_table)

    if in_degrees:
        augmented_table = augmented_table.assign(
            RA=lambda x: (np.radians(x.RA) + np.pi) % (2 * np.pi) - np.pi,
            DEC=lambda x: np.radians(x.DEC),
        )

    if "mag" not in augmented_table.columns and "magnitude" in augmented_table.columns:
        augmented_table["mag"] = augmented_table["magnitude"]

    for column in ["band", "phottype", "catalog", "note2", "mag", "discovery"]:
        if column not in augmented_table.columns:
            augmented_table[column] = None

    if "epoch_seconds_UTC" not in augmented_table.columns:
        # the entry is already computed in both the ADES and 80 columns parser
        augmented_table["epoch_seconds_UTC"] = [
            time_representation.julian_day_to_seconds_since_epoch(jd)
            for jd in list(augmented_table["epoch"])
        ]

    return augmented_table


def _spacecraft_observation_mask(table: pd.DataFrame) -> pd.Series:
    """Identify optical rows that include a complete spacecraft position."""
    if not set(SPACECRAFT_POSITION_COLUMNS).issubset(table.columns):
        return pd.Series(False, index=table.index)
    return table[SPACECRAFT_POSITION_COLUMNS].notna().all(axis=1)


def _roving_observation_mask(table: pd.DataFrame) -> pd.Series:
    """Identify optical rows that include a complete roving observatory position."""
    if not set(ROVING_POSITION_COLUMNS).issubset(table.columns):
        return pd.Series(False, index=table.index)
    return table[ROVING_POSITION_COLUMNS].notna().all(axis=1)


def _build_spacecraft_supplementary_data(table: pd.DataFrame) -> list[TrackingSupplementaryData]:
    """Receiver-state supplementary data for space-based observations.

    In 80-columns, MPC parallax records give geocentric J2000 positions at the UTC observation
    epochs; Tudat derives velocities by finite differences.
    In ADES often the velocity is also defined.
    """
    supplementary_data = []
    # sort the observations per observatory
    for observatory, group in table.loc[_spacecraft_observation_mask(table)].groupby(
        "observatory", sort=False
    ):
        # if in ADES format check if the ctr value is present; if multiple ctr are
        # defined for the same observatory raise and error. Standard practice is to have 399 as ctr,
        # exceptions are only allowed with the approval of the MPC
        # (this could be improved - allowing all spacecraft obs to be processed singularly without raising the error)
        ctr_values = group["ctr"].dropna().unique() if "ctr" in group.columns else []
        if len(ctr_values) > 1:
            raise ValueError(f"Multiple ctr values for observatory {observatory}: {ctr_values}")

        # define the frame origin of the spacecraft's state
        if len(ctr_values) == 1:
            frame_origin = NAIF_IDS[int(ctr_values[0])]
        else:
            frame_origin = "Earth"

        receiver_data = TrackingSupplementaryData(str(observatory), "")

        # check if info about velocity is availble; if so use it to build the state history
        has_velocity_columns = all(
            column in group.columns for column in SPACECRAFT_VELOCITY_COLUMNS
        )
        columns_spacecraft = SPACECRAFT_POSITION_COLUMNS.copy()
        if has_velocity_columns:
            columns_spacecraft += SPACECRAFT_VELOCITY_COLUMNS
            state_data = group.groupby("epoch_seconds_UTC")[columns_spacecraft].mean()
            receiver_data.translational_state_supplementary_data = (
                TranslationalStateSupplementaryData(
                    state_history={
                        float(epoch): state
                        for epoch, state in zip(state_data.index, state_data.to_numpy())
                    },
                    frame_origin=frame_origin,
                    is_velocity_defined=True,
                    time_scale="UTC",
                    frame_orientation="J2000",
                )
            )
            supplementary_data.append(receiver_data)

        else:
            positions = group.groupby("epoch_seconds_UTC")[columns_spacecraft].mean()
            receiver_data.translational_state_supplementary_data = (
                TranslationalStateSupplementaryData(
                    state_history={
                        float(epoch): np.concatenate((position, np.zeros(3)))
                        for epoch, position in zip(positions.index, positions.to_numpy())
                    },
                    frame_origin=frame_origin,
                    is_velocity_defined=False,
                    time_scale="UTC",
                    frame_orientation="J2000",
                )
            )
            supplementary_data.append(receiver_data)

    return supplementary_data

    # def _build_roving_obs_supplementary_data(table: pd.DataFrame) -> list[TrackingSupplementaryData]:
    """Receiver-state supplementary data for Earth-based roving observatories observations.

    The position can be given in WGS84, ITRF, IAU frames and the type of given position
    will also depend on the reference frame.
    """

    # WGS84: pos1 = east longitude (deg), pos2 = latitude (deg(, pos3 = altitude (m) - already geodetic (as Tudat
    # accepts them) using Earth shape model WGS84

    # ITRF: pos1 = east longitude (degrees), pos2 = Rxy (km), pos3 = Rz (km) - convert to cartesian coordinates
    # to define the ground station on Earth, already body-fixed (ITRF is a body fixed reference frame)
    # x = Rxy cos L
    # y = Rxy sin L
    # z = Rz

    # IAU: pos1 = longitude (deg), pos2 = latitude (m), pos3 = altitude (m) - used to define a roving observatory
    # on a body different than Earth - currently not supported

    # tudat in the definition of a basic station allows for cartesian, spherical or geodetic body-fixed definitions


def datetime_to_utc_seconds(epoch) -> float:
    """Convert an epoch to UTC seconds since J2000.

    Parameters
    ----------
    epoch : float | DateTime | Time | datetime.datetime
        Epoch as UTC seconds since J2000, a Tudat ``DateTime`` or ``Time``, or a
        Python datetime interpreted as UTC.

    Returns
    -------
    float
        UTC seconds since J2000.
    """
    if hasattr(epoch, "to_epoch"):
        return float(epoch.to_epoch())
    if hasattr(epoch, "to_float"):
        return float(epoch.to_float())
    if hasattr(epoch, "year") and hasattr(epoch, "month") and hasattr(epoch, "day"):
        return float(time_representation.DateTime.from_python_datetime(epoch).to_epoch())
    return float(epoch)


def filter_augmented_optical_table(
    table: pd.DataFrame,
    epoch_start=None,
    epoch_end=None,
    observatories: list[str | int] | None = None,
    observatories_exclude: list[str | int] | None = None,
) -> pd.DataFrame:
    """Filter an augmented optical astrometry table.

    This is a supporting utility for inspecting or manipulating already
    augmented optical data. It is not part of the typical Tudat workflow for
    loading and processing tracking data.

    Epoch filters are interpreted as UTC seconds since J2000, or as Python
    datetimes that are converted to UTC seconds since J2000.

    Parameters
    ----------
    table : pandas.DataFrame
        Augmented optical astrometry table to filter.
    epoch_start : float | datetime.datetime | None, default None
        Optional lower epoch bound.
    epoch_end : float | datetime.datetime | None, default None
        Optional upper epoch bound.
    observatories : list[str | int] | None, default None
        Optional list of observatory identifiers to keep.
    observatories_exclude : list[str | int] | None, default None
        Optional list of observatory identifiers to remove.

    Returns
    -------
    pandas.DataFrame
        Filtered augmented optical astrometry table.
    """

    if "epoch_seconds_UTC" not in table.columns:
        raise ValueError("Expected an augmented optical table with an epoch_seconds_UTC column.")

    filtered = table.copy()
    if epoch_start is not None:
        filtered = filtered.loc[
            filtered["epoch_seconds_UTC"] >= datetime_to_utc_seconds(epoch_start)
        ]
    if epoch_end is not None:
        filtered = filtered.loc[filtered["epoch_seconds_UTC"] <= datetime_to_utc_seconds(epoch_end)]

    if observatories is not None:
        included = {str(observatory).strip().zfill(3) for observatory in observatories}
        filtered = filtered.loc[filtered["observatory"].isin(included)]
    if observatories_exclude is not None:
        excluded = {str(observatory).strip().zfill(3) for observatory in observatories_exclude}
        filtered = filtered.loc[~filtered["observatory"].isin(excluded)]

    return filtered.reset_index(drop=True)


def _resolve_optical_target_names(table):
    # Resolve names per MPC identifier, across all observatories. Missing names
    # inherit the identifier's supplied name; conflicting names are ambiguous.
    target_names = {}
    for target, group in table.groupby("number", sort=False):
        names = group.get("custom_name", pd.Series(dtype=str)).dropna().astype(str).str.strip()
        names = names[names != ""].unique()
        if len(names) > 1:
            raise ValueError(f"Conflicting custom names for MPC target {target}: {names.tolist()}")
        target_names[target] = names[0] if len(names) else str(target)

    return target_names


def _apply_ades_weights(group, tracking_data_object):
    required_columns = ["rmsRA", "rmsDec"]
    if group[required_columns].isna().any().any():
        raise ValueError("ADES weighing requires non-missing rmsRA and rmsDec values")
    weightsRA = (
        (
            (pd.to_numeric(group["rmsRA"]).to_numpy(dtype=float) * u.arcsec)
            / np.cos(
                pd.to_numeric(group["DEC"]).to_numpy(dtype=float)
            )  # the value of DEC is already in radians
        ).to_value(u.rad)
    ) ** (
        -2
    )  # scale the uncertainty by cos(Dec) and transform in radians
    weightsDec = (pd.to_numeric(group["rmsDec"]).to_numpy(dtype=float) * u.arcsec).to_value(
        u.rad
    ) ** (-2)

    return weightsRA, weightsDec


def _apply_ades_weights_per_night(group, tracking_data_object):
    # apply ADES weights and de-weight by factor sqrt(N) if N observations have been taken during the same night
    required_columns = ["rmsRA", "rmsDec"]
    if group[required_columns].isna().any().any():
        raise ValueError("ADES weighing requires non-missing rmsRA and rmsDec values")

    # Count observations for each calendar date, once for every row.
    dates = group["epoch"].str[:10]
    observations_per_night = dates.map(dates.value_counts()).tolist()

    weightsRA = []
    weightsDec = []
    for count, rms_ra, rms_dec, dec in zip(
        observations_per_night, group["rmsRA"], group["rmsDec"], group["DEC"]
    ):
        rms_weights_ra = (
            (float(rms_ra) * u.arcsec)
            / np.cos(float(dec))  # the value of DEC is already in radians
        ).to_value(u.rad)

        rms_weights_dec = (float(rms_dec) * u.arcsec).to_value(u.rad)

        weightsRA.append(1 / (rms_weights_ra**2 * np.sqrt(count)))
        weightsDec.append(1 / (rms_weights_dec**2 * np.sqrt(count)))

    return weightsRA, weightsDec


def optical_table_to_tracking_data(
    table: pd.DataFrame,
    weighing_scheme: str | None = "",  # add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
):
    """Convert an augmented optical table to TrackingData and supplementary data lists.

    This is a supporting conversion utility used by the optical readers. In the
    typical Tudat workflow, call :func:`read_optical_data` or a source-specific
    optical reader instead.

    Parameters
    ----------
    table : pandas.DataFrame
        Augmented optical astrometry table to convert.
    add_weights : bool, default False
        Whether to assign the default optical weighing scheme of :cite:p:`veres2017`.
    add_star_catalog_corrections : bool, default False
        Whether to attach star-catalog bias corrections following :cite:p:`eggl2020`.
    add_ancillary_data : bool, default False
        Whether to attach available optical ancillary metadata.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.
    """
    table = create_augmented_optical_table(table, in_degrees=False)
    # weighing_scheme = "VFCC17" if add_weights else ""

    if add_star_catalog_corrections:
        RA_corr, DEC_corr = get_biases_EFCC18(mpc_table=table)
        table = table.assign(_RA_corr=RA_corr, _DEC_corr=DEC_corr)

    metadata_columns = set(ANCILLARY_STRING_COLUMNS) if add_ancillary_data else set()
    # add data required to adopt the VFCC17 weighing scheme
    if weighing_scheme == "VFCC17":
        metadata_columns.update(["note2", "catalog"])
    metadata_columns.add("number")

    target_names = _resolve_optical_target_names(table)

    spacecraft_mask = _spacecraft_observation_mask(table)
    table = table.assign(_is_spacecraft_observation=spacecraft_mask.to_numpy(dtype=bool))
    supplementary_data = _build_spacecraft_supplementary_data(table)

    """
    To develop - currently no tracking data structure to store ground station information
    Needed to be able to work with Earth-based roving observatories 
    """
    # roving_mask = _roving_observation_mask(table)
    # table = table.assign(_is_roving_observation=roving_mask.to_numpy(dtype=bool))
    # add supplementary data for the roving observatory location
    # TO DO! WHERE TO STORE THEM????

    tracking_data_objects = []
    for (target, observatory, is_spacecraft), group in table.groupby(
        ["number", "observatory", "_is_spacecraft_observation"]
    ):
        observable_type, reference_link_end_type = "AngularPosition", "receiver"

        # Ground astrometry uses an Earth ground-station receiver. Space-based
        # astrometry uses the observatory code as the receiver body and attaches
        # its translational state through TrackingSupplementaryData.
        receiver_link_end = (
            (str(observatory), "") if bool(is_spacecraft) else ("Earth", str(observatory))
        )
        link_ends = [
            ((target_names[target], ""), "transmitter"),
            (receiver_link_end, reference_link_end_type),
        ]
        observations = [np.array([ra, dec]) for ra, dec in zip(group["RA"], group["DEC"])]

        tracking_data_object = TrackingData(
            observable_type=observable_type,
            link_ends=link_ends,
            observations=observations,
            epochs=group["epoch_seconds_UTC"],
            reference_link_end=reference_link_end_type,
            time_scale="UTC",
            weighing_scheme=weighing_scheme,
        )

        # currently only applicabile if all observations in the dataset have available rmaRA and rmsDec uncertainties
        if weighing_scheme == "ADES":
            weightsRA, weightsDec = _apply_ades_weights(group, tracking_data_object)
            tracking_data_object.set_observation_weights(
                [
                    np.array([ra_weight, dec_weight], dtype=np.float64)
                    for ra_weight, dec_weight in zip(weightsRA, weightsDec)
                ]
            )

        if weighing_scheme == "ADES_NIGHT":
            weightsRA, weightsDec = _apply_ades_weights_per_night(group, tracking_data_object)
            tracking_data_object.set_observation_weights(
                [
                    np.array([ra_weight, dec_weight], dtype=np.float64)
                    for ra_weight, dec_weight in zip(weightsRA, weightsDec)
                ]
            )

        if add_star_catalog_corrections:
            corrections_list = [
                -np.array([ra_c, dec_c])
                for ra_c, dec_c in zip(group["_RA_corr"], group["_DEC_corr"])
            ]
            tracking_data_object.set_observation_corrections(corrections_list)

        for column in sorted(metadata_columns):
            tracking_data_object.add_observation_metadata(
                column, group[column].fillna("").astype(str).tolist()
            )

        tracking_data_objects.append(tracking_data_object)

    return tracking_data_objects, supplementary_data


def read_optical_data(
    table: pd.DataFrame | astropy.table.QTable | astropy.table.Table,
    in_degrees: bool = True,
    frame: str = "J2000",
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
):
    """Read optical astrometry from a table into TrackingData objects.

    The table must contain the standard optical astrometry columns used by the
    Tudat optical-data pipeline. pandas and astropy tables are both accepted;
    source-specific file parsing is handled by the dedicated MPC and 80-column
    readers.
    Depending on the input type, this function dispatches to
    :func:`read_pandas_optical_data` or :func:`read_astropy_optical_data`.
    Those functions create the augmented optical table with
    :func:`create_augmented_optical_table` and convert it to Tudat optical
    tracking data with :func:`optical_table_to_tracking_data`.

    Parameters
    ----------
    table : pandas.DataFrame | astropy.table.QTable | astropy.table.Table
        Optical astrometry table to read.
    in_degrees : bool, default True
        Whether right ascension and declination are provided in degrees.
    frame : str, default "J2000"
        Reference frame of the input observations.
    custom_name : str | None, default None
        Optional target name assigned to all observations.
    add_weights : bool, default False
        Whether to assign the default optical weighing scheme of :cite:p:`veres2017`.
    add_star_catalog_corrections : bool, default False
        Whether to attach star-catalog bias corrections following :cite:p:`eggl2020`.
    add_ancillary_data : bool, default False
        Whether to attach available optical ancillary metadata.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.
    """
    if isinstance(table, pd.DataFrame):
        return read_pandas_optical_data(
            table,
            in_degrees,
            frame,
            custom_name,
            add_weights,
            add_star_catalog_corrections,
            add_ancillary_data,
        )

    if isinstance(table, (astropy.table.QTable, astropy.table.Table)):
        return read_astropy_optical_data(
            table,
            in_degrees,
            frame,
            custom_name,
            add_weights,
            add_star_catalog_corrections,
            add_ancillary_data,
        )

    raise TypeError("read_optical_data expects a pandas DataFrame or astropy Table/QTable.")


def read_pandas_optical_data(
    table: pd.DataFrame,
    in_degrees: bool = True,
    frame: str = "J2000",
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
):
    """Read optical astrometry from a pandas table into TrackingData objects.

    This is a supporting table-specific reader. In the typical Tudat workflow
    for already-loaded optical tables, call :func:`read_optical_data` instead.

    Parameters
    ----------
    table : pandas.DataFrame
        Optical astrometry table to read.
    in_degrees : bool, default True
        Whether right ascension and declination are provided in degrees.
    frame : str, default "J2000"
        Reference frame of the input observations.
    custom_name : str | None, default None
        Optional target name assigned to all observations.
    add_weights : bool, default False
        Whether to assign the default optical weighing scheme of :cite:p:`veres2017`.
    add_star_catalog_corrections : bool, default False
        Whether to attach star-catalog bias corrections following :cite:p:`eggl2020`.
    add_ancillary_data : bool, default False
        Whether to attach available optical ancillary metadata.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.
    """
    augmented_table = create_augmented_optical_table(table, in_degrees, frame, custom_name)
    return optical_table_to_tracking_data(
        augmented_table,
        add_weights,
        add_star_catalog_corrections,
        add_ancillary_data,
    )


def read_astropy_optical_data(
    table: astropy.table.QTable | astropy.table.Table,
    in_degrees: bool = True,
    frame: str = "J2000",
    custom_name: str | None = None,
    weighing_scheme: str | None = "",
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
):
    """Read optical astrometry from an astropy table into TrackingData objects.

    This is a supporting table-specific reader. In the typical Tudat workflow
    for already-loaded optical tables, call :func:`read_optical_data` instead.

    Parameters
    ----------
    table : astropy.table.QTable | astropy.table.Table
        Optical astrometry table to read.
    in_degrees : bool, default True
        Whether right ascension and declination are provided in degrees.
    frame : str, default "J2000"
        Reference frame of the input observations.
    custom_name : str | None, default None
        Optional target name assigned to all observations.
    add_weights : bool, default False
        Whether to assign the default optical weighing scheme of :cite:p:`veres2017`.
    add_star_catalog_corrections : bool, default False
        Whether to attach star-catalog bias corrections following :cite:p:`eggl2020`.
    add_ancillary_data : bool, default False
        Whether to attach available optical ancillary metadata.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.
    """
    augmented_table = create_augmented_optical_table(table, in_degrees, frame, custom_name)
    return optical_table_to_tracking_data(
        augmented_table,
        weighing_scheme,
        add_star_catalog_corrections,
        add_ancillary_data,
    )


def load_bias_file(
    filepath: str,
    Nside: int | None = None,
    catalog_flags: list | None = None,
) -> tuple[pd.DataFrame, int]:
    """Load a HEALPix star-catalog debias file.

    This is a supporting utility for inspecting or applying optical
    star-catalog bias corrections. It is not part of the typical Tudat workflow
    for loading and processing tracking data.

    The file is processed into a dataframe, and the NSIDE parameter is
    retrieved automatically when possible.
    The debias files are used for star-catalog bias corrections following
    :cite:p:`eggl2020`.

    Parameters
    ----------
    filepath : str
        Filepath of debias file.
    Nside : int | None, optional
        NSIDE value. This can usually be left as None, in which case it is
        retrieved automatically by the function.
    catalog_flags : list | None, optional
        List of catalog flags. This can usually be left at its default value.

    Returns
    -------
    tuple[pd.DataFrame, int]
        Dataframe with biases in multi-index format ((Npix x Ncat) x Nvals)
        and the NSIDE value.

    Raises
    ------
    ValueError
        If NSIDE cannot be retrieved automatically.
    """
    # auto retrieve NSIDE
    if Nside is None:
        counter = 0
        with open(filepath, "r") as file:
            while counter < 10:
                line = file.readline()
                pattern = r"! NSIDE=\s*(\d+)"
                match = re.search(pattern, line)
                if match:
                    Nside = int(match.group(1))
                    break
                counter += 1
        if Nside is None:
            raise ValueError(
                "Could not automatically retrieve NSIDE, please provide it as a parameter"
            )

    if catalog_flags is None:
        catalog_flags = DEFAULT_CATALOG_FLAGS
    catalog_flags = catalog_flags + ["unknown"]

    values = ["RA", "DEC", "PMRA", "PMDEC"]

    # create a multi_index, this effectively creates a df with 3 dimensions. [row, catalog, value]
    m_index = pd.MultiIndex.from_product(
        [catalog_flags, values],
        names=["catalog", "value"],
    )

    bias_dataframe = pd.read_csv(
        filepath,
        sep=" ",
        skiprows=23,
        skipinitialspace=True,
        index_col=None,
        header=None,
    ).iloc[:, :-1]

    # we add a set of 'unknown' columns to speed up assignment later
    len_df = bias_dataframe.shape[0]
    unknown_columns = np.zeros(shape=(len_df, 4))
    bias_dataframe[["aa", "bb", "cc", "dd"]] = unknown_columns

    # apply the multi_index
    bias_dataframe.columns = m_index
    # stack it so it goes from a Npix x (Ncat x Nvals) to (Npix x Ncat) x Nvals shape
    bias_dataframe = bias_dataframe.stack(level=0, future_stack=True)

    return bias_dataframe, Nside


def get_biases_EFCC18(
    mpc_table: pd.DataFrame,
    bias_file: str | None = BIAS_LOWRES_FILE,
    Nside: int | None = None,
    catalog_flags: list[str] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate and return star catalog bias values as described by
    :cite:t:`eggl2020`.

    This is a supporting utility for inspecting or applying optical
    star-catalog bias corrections. It is not part of the typical Tudat workflow
    for loading and processing tracking data.

    Uses the regular bias set by default. A high-resolution version of the bias
    map can be retrieved from the paper and selected using the ``bias_file``
    parameter.

    Parameters
    ----------
    mpc_table : pd.DataFrame
        Table retrieved by calling the mpc.BatchMPC.table property. Must contain
        'RA', 'DEC', 'epoch_seconds_UTC' and 'catalog' columns.
    bias_file : str | None, optional
        Optional bias file location, used to load alternative debias
        coefficients. By default, coefficients are retrieved from Tudat
        resources.
    Nside : int | None, optional
        Optional NSIDE value. This can usually be left as None.
    catalog_flags : list[str] | None, optional
        List of catalog values to use. This can usually be left at its default
        value.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Right ascension corrections and declination corrections.
    """

    if bias_file is None:
        bias_file = BIAS_LOWRES_FILE

    RA = mpc_table["RA"].to_numpy()
    DEC = np.asarray(mpc_table["DEC"].to_numpy(), dtype=float)  # mpc_table["DEC"].to_numpy()
    if "epoch_seconds_UTC" in mpc_table.columns:
        epochs = mpc_table["epoch_seconds_UTC"].to_numpy()
    else:
        epochs = np.array(
            [
                time_representation.julian_day_to_seconds_since_epoch(jd)
                for jd in mpc_table["epoch"].to_numpy()
            ]
        )
    catalog = mpc_table["catalog"].to_numpy()

    # load bias file
    # index matches the pixels
    # this is effectively a 3d table with axes: (pixel, star catalog), value) using pandas multiindex
    bias_df, nside = load_bias_file(filepath=bias_file, Nside=Nside, catalog_flags=catalog_flags)

    # find nearest tile using HEALPix Algorithm and get indices
    hp_obj = HEALPix(nside=nside)

    pixels = hp_obj.lonlat_to_healpix(Quantity(RA, unit=u.rad), Quantity(DEC, unit=u.rad))

    # retrieve bias values from bias file using indices
    # result is N x 4 biases for the correct star catalog
    all_catalog_ids = bias_df.index.levels[1].to_list()
    # this changes all ids not present in the bias file to unknown, resulting in zero bias
    catalog = ["unknown" if (cat not in all_catalog_ids) else cat for cat in catalog]

    # create combinations of pixel id and catalog then retrieve biases
    targets = [(pix, cat) for pix, cat in zip(pixels, catalog)]
    biases = bias_df.loc[targets, ["RA", "DEC", "PMRA", "PMDEC"]].to_numpy()

    # same as find_orb -> bias.cpp
    # https://github.com/Bill-Gray/find_orb/blob/master/bias.cpp#L213
    epochs_years = [
        time_representation.seconds_since_epoch_to_julian_years_since_epoch(epoch)
        for epoch in epochs
    ]

    # from the bias file readme.txt:
    RA_correction = biases[:, 0] + (epochs_years * (biases[:, 2] / 1000))
    RA_correction = RA_correction / np.cos(DEC)  # DEC here in radians because of cosine
    DEC_correction = biases[:, 1] + (epochs_years * (biases[:, 3] / 1000))

    # convert from arcsec to radians
    RA_correction = Quantity(RA_correction, unit=u.arcsec).to(u.rad).value
    DEC_correction = Quantity(DEC_correction, unit=u.arcsec).to(u.rad).value

    return RA_correction, DEC_correction
