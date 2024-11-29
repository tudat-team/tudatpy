from tudatpy.numerical_simulation import environment_setup  # type:ignore
from tudatpy.numerical_simulation import estimation, environment  # type:ignore
from tudatpy.numerical_simulation.estimation_setup import observation  # type:ignore
from tudatpy.numerical_simulation.environment_setup import add_gravity_field_model
from tudatpy.numerical_simulation.environment_setup.gravity_field import central_sbdb

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.dates as mdates
import datetime

from astroquery.mpc import MPC
from astropy.time import Time
from astropy_healpix import HEALPix
from astropy.units import Quantity
from astropy.time import Time
import astropy.units as u
import astropy

from typing import Union, Tuple, List, Dict

import copy

import os
import re


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


def load_bias_file(
    filepath: str,
    Nside: Union[None, int] = None,
    catalog_flags: list = DEFAULT_CATALOG_FLAGS,
) -> Tuple[pd.DataFrame, int]:
    """Loads a healpix star catalog debias file and processes it into a dataframe. Automatically retrieves NSIDE parameter.

    Parameters
    ----------
    filepath : str
        Filepath of debias file.
    Nside : Union[None, int], optional
        NSIDE value, to be left None in most cases as this is retrieved automatically by the function, by default None
    catalog_flags : Union[None, list], optional
        list of catalog flags, should be left default in most cases, by default None

    Returns
    -------
    Tuple[pd.DataFrame, int]
        Dataframe with biases in multiindex format ((Npix x Ncat) x Nvals), the numpix value

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
    bias_dataframe = bias_dataframe.stack(level=0)

    return bias_dataframe, Nside


def get_biases_EFCC18(
    RA: Union[float, np.ndarray, list],
    DEC: Union[float, np.ndarray, list],
    epochJ2000secondsTDB: Union[float, np.ndarray, list],
    catalog: Union[str, np.ndarray, list],
    bias_file: Union[str, None] = BIAS_LOWRES_FILE,
    Nside: Union[int, None] = None,
    catalog_flags: List[str] = DEFAULT_CATALOG_FLAGS,
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate and return star catalog bias values as described in:
    "Star catalog position and proper motion corrections in asteroid astrometry II: The Gaia era" by Eggl et al. (2018).
    Uses the regular bias set by default. A high res version of the bias map can be retrieved from the paper.
    This can then be selected using the bias_file paramater.

    Parameters
    ----------
    RA : Union[float, np.ndarray, list]
        Right Ascension value in radians
    DEC : Union[float, np.ndarray, list]
        Declination value in radians
    epochJ2000secondsTDB : Union[float, np.ndarray, list]
        Time in seconds since J2000 TDB.
    catalog : Union[str, np.ndarray, list]
        Star Catalog code as described by MPC: https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
    bias_file : Tuple[str, None], optional
        Optional bias file location, used to load in alternative debias coefficients. By default coefficients are retrieved from Tudat resources, by default None
    Nside : Tuple[int, None], optional
        Optional Nside value, should be left None in most cases, by default None
    catalog_flags : Tuple[List[str], None], optional
        List of catalog values to use, should be left None in most cases, by default None

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Right Ascencion Corrections, Declination corrections

    Raises
    ------
    ValueError
        If all mandatory inputs are not matching in size.
    """

    if bias_file is None:
        bias_file = BIAS_LOWRES_FILE

    # transform input to numpy arrays
    if not isinstance(RA, np.ndarray):
        RA = np.array([RA]).flatten()
    if not isinstance(DEC, np.ndarray):
        DEC = np.array([DEC]).flatten()
    if not isinstance(epochJ2000secondsTDB, np.ndarray):
        epochJ2000secondsTDB = np.array([epochJ2000secondsTDB]).flatten()
    if not isinstance(catalog, np.ndarray):
        catalog = np.array([catalog]).flatten()

    if not (len(RA) == len(DEC) == len(epochJ2000secondsTDB) == len(catalog)):
        raise ValueError("All inputs must have same size")

    # load bias file
    # index matches the pixels
    # this is effectively a 3d table with axes: (pixel, star catalog), value) using pandas multiindex
    bias_df, nside = load_bias_file(
        filepath=bias_file, Nside=Nside, catalog_flags=catalog_flags
    )

    # find nearest tile using HEALPix Algorithm and get indices
    # ideally nside should be retrieved from the load_bias_file function
    hp_obj = HEALPix(nside=nside)

    pixels = hp_obj.lonlat_to_healpix(
        Quantity(RA, unit=u.rad), Quantity(DEC, unit=u.rad)
    )

    # retrieve bias values from bias file using indices
    # result is N x 4 biases for the correct star catalog
    all_catalog_ids = bias_df.index.levels[1].to_list()
    # this changes all ids not present in the bias file to unknown, resulting in zero bias
    catalog = ["unknown" if (cat not in all_catalog_ids) else cat for cat in catalog]

    # create combinations of pixel id and catalog then retrieve biases
    targets = [(pix, cat) for pix, cat in zip(pixels, catalog)]
    biases = bias_df.loc[targets, ["RA", "DEC", "PMRA", "PMDEC"]].to_numpy()

    # transform the values based on the paper, time to JD TT
    J2000_jd = 2451545.0
    epochs = (epochJ2000secondsTDB / 86400) + J2000_jd
    epochs = Time(epochs, format="jd", scale="tdb")

    # this was taken from find_orb -> bias.cpp
    # https://github.com/Bill-Gray/find_orb/blob/master/bias.cpp#L213
    epochs_years = (epochs.tt.value - J2000_jd) / 365.25

    # from the bias file readme.txt:
    RA_correction = biases[:, 0] + (epochs_years * (biases[:, 2] / 1000))
    RA_correction = RA_correction / np.cos(DEC)  # DEC here in radians because of cosine
    DEC_correction = biases[:, 1] + (epochs_years * (biases[:, 3] / 1000))

    # convert from arcsec to radians
    RA_correction = Quantity(RA_correction, unit=u.arcsec).to(u.rad).value
    DEC_correction = Quantity(DEC_correction, unit=u.arcsec).to(u.rad).value

    return RA_correction, DEC_correction


def get_weights_VFCC17(
    MPC_codes: Union[pd.Series, list, np.ndarray, None] = None,
    epochUTC: Union[pd.Series, list, np.ndarray, None] = None,
    observation_type: Union[pd.Series, list, np.ndarray, None] = None,
    observatory: Union[pd.Series, list, np.ndarray, None] = None,
    star_catalog: Union[pd.Series, list, np.ndarray, None] = None,
    mpc_table: Union[pd.DataFrame, None] = None,
    return_full_table=False,
) -> Union[np.ndarray, pd.DataFrame]:
    """Retrieves observation weights using the weighting scheme presented in
    "Statistical analysis of astrometric errors for the most productive
    asteroid surveys" by Veres et al. (2017). Input may be provided using
    either a full MPC table (e.g. from BatchMPC) or using the individual
    variables.

    Observation types: "x", "X", "V", "v", "W", "w", "R", "r", "Q", "q", "O",
    are not described by the paper and receive a placeholder weight of 1/100
    if provided.

    Parameters
    ----------
    MPC_codes : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with the MPC target codes, e.g. 433 for Eros. Size must match
        other iterables, by default None
    epochUTC : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with UTC times. Size must match other iterables, by default None
    observation_type : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with the observation types in MPC format.
        See the NOTE2 section of the MPC format description for the exact encoding:
        https://minorplanetcenter.net/iau/info/OpticalObs.html.
        Size must match other iterables, by default None
    observatory : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with the MPC target codes, e.g. 433 for Eros.
        Size must match other iterables, by default None
    star_catalog : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with the star catalog codes.
        See the MPC catalog codes page for the exact encoding:
        https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html.
        Size must match other iterables, by default None
    mpc_table : Union[pd.DataFrame, None], optional
        Table retrieved by calling the mpc.BatchMPC.table property.
        Set None when using iterable input.
        Set others None when using table, by default None
    return_full_table : bool, optional
        Return the table with all intermediate calculations if True,
        return a numpy array if False, by default False

    Returns
    -------
    np.ndarray
        If `return_full_table` is False, numpy array with weights with same size as input.
    pd.DataFrame
        If `return_full_table` is True, pandas table with all intermediate calculations.

    Raises
    ------
    ValueError
        MPC_codes, epochUTC, observation_type, observatory and star_catalog must all
        be not None and the same size. mpc_table must be None.
        If table input is used, the remaining input parameters must be done.
    """

    # Input handling
    if (
        (mpc_table is None)
        and (epochUTC is not None)
        and (observation_type is not None)
        and (observatory is not None)
        and (star_catalog is not None)
    ):
        if not (
            len(epochUTC)
            == len(observation_type)
            == len(observatory)
            == len(star_catalog)
        ):
            raise ValueError("All inputs must have same size")

        table_dict = {
            "number": MPC_codes,
            "epochUTC": epochUTC,
            "note2": observation_type,
            "observatory": observatory,
            "catalog": star_catalog,
        }
        table = pd.DataFrame.from_dict(table_dict)
    elif (
        (mpc_table is not None)
        and (epochUTC is None)
        and (observation_type is None)
        and (observatory is None)
        and (star_catalog is None)
    ):
        table = mpc_table.copy()
    else:
        raise ValueError(
            "Must provide either parameters: `epochUTC`, `observation_type`, `observatory` and `star_catalog` OR `mpc_table`."
        )

    # create col for Julian Day and the inverted weight
    table = table.assign(epochJD=lambda x: Time(x.epochUTC).jd1 + Time(x.epochUTC).jd2)
    # NOTE 1000 is a placeholder. The following observation types are not processed and receive the placeholder value:
    # first_discoveries = ["x", "X"], roaming = ["V", "v", "W", "w"], radar = ["R", "r", "Q", "q"], offset = ["O"]
    table = table.assign(inv_w=lambda _: 1000)

    # get an approximate timezone based on the observatory code's longitude
    observatories_table = MPC.get_observatory_codes().to_pandas()
    observatories_table = (
        observatories_table.assign(
            lon_wrapping=lambda x: (x.Longitude + 180) % 360 - 180
        )
        .assign(approx_tz=lambda x: ((x.lon_wrapping / 180) * 12))
        .assign(jd_tz=lambda x: (x.lon_wrapping / 360).fillna(0))
        .loc[:, ["Code", "approx_tz", "jd_tz"]]
    )

    # the reset_index + set_index is a cheeky way to keep the original table index
    # this prevents mismatches when adding weights to the tables in BatchMPC.to_tudat()
    table = (
        pd.merge(
            how="left",
            left=table.reset_index(),
            right=observatories_table,
            left_on="observatory",
            right_on="Code",
        )
        .drop("Code", axis=1)
        .set_index("index")
    )

    # create column that modifies the JD with an approximate timezone.
    # the time JD.50 will then be the approximate midnight at that timezone.
    # if we then take the int, we can group by this number to get all observations that night.
    table = table.assign(epochJD_tz_int=lambda x: np.floor(x.epochJD + x.jd_tz))
    # table = table.assign(epochJD_tz_int2=lambda x: np.round(x.epochJD + x.jd_tz, 2))

    # Below are the weights applied as described per table in:
    # https://www.sciencedirect.com/science/article/pii/S0019103517301987

    # NON-CCD
    # ###################
    # TABLE 5: Non-CCD residuals
    # Conditions for Table 5
    # 1890-1-1 and 1950-1-1
    pre_1890 = table.epochJD <= 2411368.0
    between_1890_1950 = (table.epochJD > 2411368.0) & (table.epochJD <= 2433282.0)
    after_1950 = table.epochJD > 2433282.0

    photographic = table.note2.isin([np.nan, "P", "A", "N", "Z"])
    occultations = table.note2 == "E"
    hipparcos = table.note2 == "H"
    transit_circle = table.note2 == "T"
    encoder = table.note2 == "e"
    micrometer = table.note2 == "M"
    satellite = table.note2.isin(["S", "s"])
    multinormal_place = table.note2 == "n"

    # Apply Table 5
    table.inv_w = table.inv_w.mask((photographic & pre_1890), 10.0)
    table.inv_w = table.inv_w.mask((photographic & between_1890_1950), 5.0)
    table.inv_w = table.inv_w.mask((photographic & after_1950), 2.5)

    table.inv_w = table.inv_w.mask((occultations), 0.2)
    table.inv_w = table.inv_w.mask((hipparcos), 0.2)
    table.inv_w = table.inv_w.mask((transit_circle), 0.5)
    table.inv_w = table.inv_w.mask((encoder), 0.75)
    table.inv_w = table.inv_w.mask((micrometer), 2.0)
    table.inv_w = table.inv_w.mask((satellite), 1.5)
    table.inv_w = table.inv_w.mask((multinormal_place), 1.0)
    # ###################

    # CCD
    # ###################
    # TABLE 3: astroid observers
    # Table 3 conditions:
    # NOTE for now CCD will receive the same base weighting as CCD.
    # CMOS sensors for astronomy is a new development.
    ccd = table.note2.isin(["C", "c", "D"])
    cmos = table.note2.isin(["B"])
    ccd = ccd | cmos

    tab3_no_catalog = table.catalog.isin(["unknown", np.nan])

    # Apply Table 3:
    table.inv_w = table.inv_w.mask((ccd & ~tab3_no_catalog), 1.0)
    table.inv_w = table.inv_w.mask((ccd & tab3_no_catalog), 1.5)

    table.inv_w = table.inv_w.mask((table.observatory == "704"), 1.0)
    table.inv_w = table.inv_w.mask((table.observatory == "G96"), 0.5)
    table.inv_w = table.inv_w.mask((table.observatory == "F51"), 0.2)
    table.inv_w = table.inv_w.mask((table.observatory == "G45"), 0.6)
    table.inv_w = table.inv_w.mask((table.observatory == "699"), 0.8)
    table.inv_w = table.inv_w.mask((table.observatory == "D29"), 0.75)

    table.inv_w = table.inv_w.mask((table.observatory == "C51"), 1.0)
    table.inv_w = table.inv_w.mask((table.observatory == "E12"), 0.75)
    table.inv_w = table.inv_w.mask((table.observatory == "608"), 0.6)
    table.inv_w = table.inv_w.mask((table.observatory == "J75"), 1.0)
    # ###################

    # ###################
    # TABLE 2: epoch dependant residuals
    # Table 2 conditions:
    tab2_703_epoch = table.epochJD < 2456658.0  # 2014-1-1
    tab2_691_epoch = table.epochJD < 2452640.0  # 2003-1-1
    tab2_644_epoch = table.epochJD < 2452883.0  # 2003-9-1

    # Apply Table 2:
    table.inv_w = table.inv_w.mask((tab2_703_epoch & (table.observatory == "703")), 1.0)
    table.inv_w = table.inv_w.mask(
        (~tab2_703_epoch & (table.observatory == "703")), 0.8
    )
    table.inv_w = table.inv_w.mask((tab2_691_epoch & (table.observatory == "691")), 0.6)
    table.inv_w = table.inv_w.mask(
        (~tab2_691_epoch & (table.observatory == "691")), 0.5
    )
    table.inv_w = table.inv_w.mask((tab2_644_epoch & (table.observatory == "644")), 0.6)
    table.inv_w = table.inv_w.mask(
        (~tab2_644_epoch & (table.observatory == "644")), 0.4
    )
    # ###################

    # ###################
    # TABLE 4: Catalog dependant + NEO observers
    # Table 4 conditions:
    # NOTE these are the LCO observatories described in the paper
    LCO_original = [
        "K92",
        "K93",
        "Q63",
        "Q64",
        "V37",
        "W84",
        "W85",
        "W86",
        "W87",
        "K91",
        "E10",
        "F65",
    ]
    # NOTE these are additional LCO observatories, mostly online after publication
    # Aqawan and Clamshell (0.4m) style observatories have not been included due to their previous omssion in the paper.
    # Those thus assume the default inv_weight of 1.0
    # Since remaining 1m telescopes are duplicates of existing observatories, they are included.
    # Dates retrieved from: https://sbnmpc.astro.umd.edu/mpecwatch/obs.html
    # See also: https://lco.global/observatory/sites/mpccodes/
    LCO_new = [
        # "L09",  # online 2018 Aqawan
        # "Q58",  # online 2017 Clamshell
        # "Q59",  # online 2017 Clamshell
        # "T03",  # online 2017 Clamshell
        # "T04",  # online 2016 Clamshell
        # "V38",  # online 2018 Aqawan
        "V39",  # online 2019 1m
        # "W79",  # online 2018 Aqawan
        # "W89",  # online 2015 Aqawan
        # "Z17",  # online 2017 Aqawan
        # "Z21",  # online 2015 Aqawan
        "Z24",  # online 2021 1m
        "Z31",  # online 2021 1m
    ]

    LCO_obs = LCO_original + LCO_new
    tab4_LCO_observatories = table.observatory.isin(LCO_obs)
    # NOTE see catalog codes info: https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
    tab4_Catalog_UCAC4 = table.catalog == "q"
    tab4_Catalog_PPMXL = table.catalog == "t"
    tab4_Catalog_GAIA = table.catalog.isin(["U", "V", "W", "X", "3", "6"])
    tab4_Catalog_USNOB12 = table.catalog.isin(["o", "s"])

    tab4_G83_UCAC4_PPMXL = (
        (table.observatory == "G83") & tab4_Catalog_UCAC4 & tab4_Catalog_PPMXL
    )
    tab4_G83_GAIA = (table.observatory == "G83") & tab4_Catalog_GAIA

    tab4_Y28_GAIA_PPMXL = (
        (table.observatory == "Y28") & tab4_Catalog_PPMXL & tab4_Catalog_GAIA
    )
    tab4_568_USNOB = (table.observatory == "568") & tab4_Catalog_USNOB12
    tab4_568_GAIA = (table.observatory == "568") & tab4_Catalog_GAIA
    tab4_568_PPMXL = (table.observatory == "568") & tab4_Catalog_PPMXL
    tab4_T09_T12_T14_GAIA = (table.observatory == "568") & tab4_Catalog_GAIA
    tab4_309_UCAC4_PPMXL = (
        (table.observatory == "568") & tab4_Catalog_UCAC4 & tab4_Catalog_PPMXL
    )
    tab4_309_GAIA = (table.observatory == "568") & tab4_Catalog_GAIA

    # Apply Table 4:
    table.inv_w = table.inv_w.mask(table.observatory == "645", 0.3)
    table.inv_w = table.inv_w.mask(table.observatory == "673", 0.3)
    table.inv_w = table.inv_w.mask(table.observatory == "689", 0.5)
    table.inv_w = table.inv_w.mask(table.observatory == "950", 0.5)
    table.inv_w = table.inv_w.mask(table.observatory == "H01", 0.3)
    table.inv_w = table.inv_w.mask(table.observatory == "J04", 0.4)
    table.inv_w = table.inv_w.mask(tab4_G83_UCAC4_PPMXL, 0.3)
    table.inv_w = table.inv_w.mask(tab4_G83_GAIA, 0.2)
    table.inv_w = table.inv_w.mask(tab4_LCO_observatories, 0.4)
    table.inv_w = table.inv_w.mask(table.observatory == "W84", 0.5)  # LCO includes W84?

    table.inv_w = table.inv_w.mask(tab4_Y28_GAIA_PPMXL, 0.3)
    table.inv_w = table.inv_w.mask(tab4_568_USNOB, 0.5)
    table.inv_w = table.inv_w.mask(tab4_568_GAIA, 0.1)
    table.inv_w = table.inv_w.mask(tab4_568_PPMXL, 0.2)
    table.inv_w = table.inv_w.mask(tab4_T09_T12_T14_GAIA, 0.1)
    table.inv_w = table.inv_w.mask(tab4_309_UCAC4_PPMXL, 0.3)
    table.inv_w = table.inv_w.mask(tab4_309_GAIA, 0.2)
    # ###################

    # Transform residual into weight:
    table = table.assign(
        weight_pre=lambda x: 1
        / np.square(Quantity(x.inv_w, unit=u.arcsec).to(u.rad).value)
    )

    # Reduce weight if there are more than 4 observations that night:
    # NOTE For satellites, this is done per julian day.
    # calculate number of observations on one night at one observatory for one target:
    table = table.assign(
        observations_on_epoch=lambda x: x.groupby(
            ["epochJD_tz_int", "observatory", "number"]
        ).epochUTC.transform("count")
    )
    # divide weight by sqrt(N/4) if N > 4 else 1
    table = table.assign(
        mult_obs_deweight=lambda x: np.maximum(
            np.sqrt(x.observations_on_epoch / 4), 1.0
        )
    )

    table = table.assign(weight=lambda x: x.weight_pre / x.mult_obs_deweight)

    if return_full_table:
        return table
    else:
        return table.weight.to_numpy()


class BatchMPC:
    """This class provides an interface between observations
    in the Minor Planet Centre database and Tudat.

    Notes
    ----------
    Currently, transformations between reference frames are not implemented.
    As such observations are only returned in the J2000 Frame.

    Examples
    ----------
    Basic sequence of usage:

    Initialise and retrieve data:

    >>> MPCcodes = [1, 4] # Ceres and Vesta
    >>> batch = BatchMPC()
    >>> batch.get_observations(MPCcodes)
    >>> batch.filter(epoch_start=datetime.datetime(2016, 1, 1))

    Transform to Tudat format:
    >>> ...
    >>> bodies = environment_setup.create_system_of_bodies(body_settings)
    >>> observation_collection = batch.to_tudat(bodies=bodies, included_satellites=None)

    """

    def __init__(self) -> None:
        """Create an empty MPC batch."""
        self._table: pd.DataFrame = pd.DataFrame()
        self._observatories: List[str] = []
        self._space_telescopes: List[str] = []
        self._bands: List[str] = []
        self._MPC_codes: List[str] = []
        self._size: int = 0

        self._observatory_info: Union[pd.DataFrame, None] = None
        self._MPC_space_telescopes: List[str] = []

        self._get_station_info()

        self._epoch_start: float = 0.0
        self._epoch_end: float = 0.0

        self._bodies_created = {}

        self._EFCC18_applied = False
        self._custom_weights_set = False

        # for manual additions of table (from_pandas, from_astropy)
        self._req_cols = ["number", "epoch", "RA", "DEC", "band", "observatory"]

    def __copy__(self):
        new = BatchMPC()

        new._table = copy.deepcopy(self._table)
        new._refresh_metadata()

        return new

    def copy(self) -> "BatchMPC":
        """Create a copy of the batch, equivalent to copy.copy(batchMPC())

        Returns
        -------
        BatchMPC
            Copy of batch.
        """
        return copy.copy(self)

    # getters to make everything read-only
    @property
    def table(self) -> pd.DataFrame:
        """Pandas dataframe with observation data"""
        return self._table

    @property
    def observatories(self) -> List[str]:
        """List of observatories in batch"""
        return self._observatories

    @property
    def space_telescopes(self) -> List[str]:
        """List of satellite_observatories in batch"""
        return self._space_telescopes

    @property
    def MPC_objects(self) -> List[str]:
        """List of MPC objects"""
        return self._MPC_codes

    @property
    def size(self) -> int:
        """Number of observations in batch"""
        return self._size

    @property
    def bands(self) -> List[str]:
        """List of bands in batch"""
        return self._bands

    @property
    def epoch_start(self) -> float:
        """Epoch of oldest observation in batch in seconds since J2000 TDB"""
        return self._epoch_start

    @property
    def epoch_end(self) -> float:
        """Epoch of latest observation in batch in seconds since J2000 TDB"""
        return self._epoch_end

    @property
    def bodies_created(self) -> dict:
        """Dictionary with the bodies created by to_tudat and details."""
        return self._bodies_created

    def __len__(self):
        return self._size

    def __add__(self, other):
        temp = BatchMPC()

        temp._table = (
            pd.concat([self._table, other._table])
            .sort_values("epoch")
            .drop_duplicates()
        )

        temp._refresh_metadata()

        return temp

    # helper functions
    def _refresh_metadata(self) -> None:
        """Internal. Update batch metadata"""
        self._table.drop_duplicates()

        self._observatories = list(self._table.observatory.unique())
        self._space_telescopes = [
            x for x in self._observatories if x in self._MPC_space_telescopes
        ]
        if "bands" in self._table.columns:
            self._bands = list(self._table.band.unique())
        self._MPC_codes = list(self._table.number.unique())
        self._size = len(self._table)

        self._epoch_start = self._table.epochJ2000secondsTDB.min()
        self._epoch_end = self._table.epochJ2000secondsTDB.max()

    def _get_station_info(self) -> None:
        """Internal. Retrieve data on MPC listed observatories"""
        try:
            temp = MPC.get_observatory_codes().to_pandas()  # type: ignore
            # This query checks if Longitude is Nan: non-terretrial telescopes
            sats = list(temp.query("Longitude != Longitude").Code.values)
            self._observatory_info = temp
            self._MPC_space_telescopes = sats
        except Exception as e:
            print("An error occured while retrieving observatory data")
            print(e)

    def _add_observatory_positions(
        self, bodies: environment.SystemOfBodies, earth_name
    ) -> None:
        """Internal. Add observatory cartesian postions to station data"""
        temp = self._observatory_info

        if temp is None:
            # This error will probably never occur
            txt = """Observatory positions can not be assigned.
                  This is likely due to a failure in retrieving the observatories."""
            raise ValueError(txt)

        r_earth = bodies.get(earth_name).shape_model.average_radius

        # Add geocentric cartesian positions
        temp = (
            temp.assign(X=lambda x: x.cos * r_earth * np.cos(np.radians(x.Longitude)))
            .assign(Y=lambda x: x.cos * r_earth * np.sin(np.radians(x.Longitude)))
            .assign(Z=lambda x: x.sin * r_earth)
        )
        self._observatory_info = temp

    def _apply_EFCC18(
        self,
        bias_file=None,
        Nside=None,
        catalog_flags=None,
    ):
        """
        Internal, applies star catalog biases based on
        'Star catalog position and proper motion corrections in asteroid astrometry II: The Gaia era' by Eggl et al. (2018)
        """
        if self._EFCC18_applied:
            pass
        else:
            RA_corr, DEC_corr = get_biases_EFCC18(
                RA=self.table.RA.to_numpy(),
                DEC=self.table.DEC.to_numpy(),
                epochJ2000secondsTDB=self.table.epochJ2000secondsTDB.to_numpy(),
                catalog=self.table.catalog.to_numpy(),
                bias_file=bias_file,
                Nside=Nside,
                catalog_flags=catalog_flags,
            )

            self._table = self._table.assign(RA_EFCC18=lambda x: x.RA - RA_corr)
            self._table = self._table.assign(DEC_EFCC18=lambda x: x.DEC - DEC_corr)
            self._table = self._table.assign(corr_RA_EFCC18=lambda x: RA_corr)
            self._table = self._table.assign(corr_DEC_EFCC18=lambda x: DEC_corr)

            self._EFCC18_applied = True

    # methods for data retrievels
    def get_observations(
        self, MPCcodes: List[Union[str, int]], drop_misc_observations: bool = True
    ) -> None:
        """Retrieve all observations for a set of MPC listed objeccts.
        This method uses astroquery to retrieve the observations from the MPC.
        An internet connection is required, observations are cached for faster subsequent retrieval.
        Removes duplicate and irrelevant observation data by default (see `drop_misc_observations`).

        Parameters
        ----------
        MPCcodes : List[int]
            List of integer MPC object codes for minor planets or and comets.
        drop_misc_observations : List[int]
            Drops observations made by method: radar and offset (natural satellites).
            Drops observations made by roaming observers.
            Drops duplicate listings to denote first observation.
        """

        if not isinstance(MPCcodes, list):
            raise ValueError("MPCcodes parameter must be list of integers/strings")
        for code in MPCcodes:
            if not (isinstance(code, int) or isinstance(code, str)):
                raise ValueError(
                    "All codes in the MPCcodes parameter must be integers or string"
                )

            try:
                obs = MPC.get_observations(code).to_pandas()  # type: ignore

                # convert JD to J2000 and UTC, convert deg to rad
                obs = (
                    obs.assign(
                        epochJ2000secondsTDB=lambda x: (
                            Time(x.epoch, format="jd", scale="utc").tdb.value
                            - 2451545.0
                        )
                        * 86400
                    )
                    .assign(RA=lambda x: np.radians(x.RA))
                    .assign(DEC=lambda x: np.radians(x.DEC))
                    .assign(epochUTC=lambda x: Time(x.epoch, format="jd").to_datetime())
                )

                # drop miscellenous observations:
                if drop_misc_observations:
                    first_discoveries = ["x", "X"]
                    roaming = ["V", "v", "W", "w"]
                    radar = ["R", "r", "Q", "q"]
                    offset = ["O"]
                    observation_types_to_drop = (
                        first_discoveries + roaming + radar + offset
                    )
                    obs = obs.query("note2 != @observation_types_to_drop")

                # convert object mpc code to string
                if "comettype" in obs.columns:
                    # for the case where we have a comet
                    obs.loc[:, "number"] = obs.desig
                else:
                    obs.loc[:, "number"] = obs.number.astype(str)

                self._table = pd.concat([self._table, obs])

            except Exception as e:
                print(f"An error occured while retrieving observations of MPC: {code}")
                print(e)

        self._refresh_metadata()

    def _add_table(self, table: pd.DataFrame, in_degrees: bool = True):
        """Internal. Formats a manually entered table of observations, use from_astropy or from_pandas."""
        obs = table
        obs = obs.assign(
            epochJ2000secondsTDB=lambda x: (
                Time(x.epoch, format="jd", scale="utc").tdb.value - 2451545.0
            )
            * 86400
        ).assign(epochUTC=lambda x: Time(x.epoch, format="jd").to_datetime())
        if in_degrees:
            obs = obs.assign(RA=lambda x: np.radians(x.RA)).assign(
                DEC=lambda x: np.radians(x.DEC)
            )

        # convert object mpc code to string
        obs["number"] = obs.number.astype(str)
        self._table = pd.concat([self._table, obs])
        self._refresh_metadata()

    def from_astropy(
        self, table: astropy.table.QTable, in_degrees: bool = True, frame: str = "J2000"
    ):
        """Manually input an astropy table with observations into the batch. 
        Usefull when manual filtering from astroquery is required
        Table must contain the following columns:
        number - MPC code of the minor planet\\
        epoch - in julian days\\
        RA - right ascension in either degrees or radians (degrees is default)\\
        DEC - declination in either degrees or radians (degrees is default)\\
        band - band of the observation (currently unused)\\
        observatory - MPC code of the observatory

        Parameters
        ----------
        table : astropy.table.QTable
            Astropy table with the observations
        in_degrees : bool, optional
            if True RA and DEC are assumed in degrees, else radians, by default True
        frame : str, optional
            Reference frame. Please note that only J2000 is currently supported
            , by default "J2000"
        """
        if not (
            isinstance(table, astropy.table.QTable)
            or isinstance(table, astropy.table.Table)
        ):
            raise ValueError(
                "Table must be of type astropy.table.QTable or astropy.table.Table"
            )
        if frame != "J2000":
            txt = "Only observations in J2000 are supported currently"
            raise NotImplementedError(txt)

        # check if all mandatory names are present
        if not set(self._req_cols).issubset(set(table.colnames)):
            txt = f"Table must include a set of mandatory columns: {self._req_cols}"
            raise ValueError(txt)

        self._add_table(table=table.to_pandas(), in_degrees=in_degrees)

    def from_pandas(
        self, table: pd.DataFrame, in_degrees: bool = True, frame: str = "J2000"
    ):
        """Manually input an pandas dataframe with observations into the batch. 
        Usefull when manual filtering from astroquery is required
        Table must contain the following columns:
        number - MPC code of the minor planet\\
        epoch - in julian days\\
        RA - right ascension in either degrees or radians (degrees is default)\\
        DEC - declination in either degrees or radians (degrees is default)\\
        band - band of the observation (currently unused)\\
        observatory - MPC code of the observatory

        Parameters
        ----------
        table : astropy.table.QTable
            Astropy table with the observations
        in_degrees : bool, optional
            if True RA and DEC are assumed in degrees, else radians, by default True
        frame : str, optional
            Reference frame. Please not that only J2000 is currently supported
            , by default "J2000"
        """
        if not isinstance(table, pd.DataFrame):
            raise ValueError("Table must be of type pandas.DataFrame")
        if frame != "J2000":
            txt = "Only observations in J2000 are supported currently"
            raise NotImplementedError(txt)

        # check if all mandatory names are present
        if not set(self._req_cols).issubset(set(table.columns)):
            txt = f"Table must include a set of mandatory columns: {self._req_cols}"
            raise ValueError(txt)

        self._add_table(table=table, in_degrees=in_degrees)

    def set_weights(
        self,
        weights: Union[list, np.ndarray, pd.Series],
    ):
        """Manually set weights per observation. Weights are passed to
        observation collection when `.to_tudat()` is called. Set the
        `apply_weights_VFCC17` parameter in `.to_tudat()` to `False` to avoid
        overwriting. The order of the weights should match the order found in
        the `.table` parameter.

        Parameters
        ----------
        weights : Union[list, np.ndarray, pd.Series]
            Iterable with weights per observation.

        Raises
        ------
        ValueError
            If the size of weights does not match the number of observations in
            the batch table.
        """
        if len(weights) != len(self.table):
            raise ValueError(
                "Weights vector must have same length as number of observations (BatchMPC.table)."
            )
        if isinstance(weights, pd.Series):
            # this is to avoid errors with indexing in pandas.
            weights = weights.to_numpy()

        self._table.loc[:, "weight"] = weights
        self._custom_weights_set = True

    def filter(
        self,
        bands: Union[List[str], None] = None,
        catalogs: Union[List[str], None] = None,
        observation_types: Union[List[str], None] = None,
        observatories: Union[List[str], None] = None,
        observatories_exclude: Union[List[str], None] = None,
        epoch_start: Union[float, datetime.datetime, None] = None,
        epoch_end: Union[float, datetime.datetime, None] = None,
        in_place: bool = True,
    ) -> Union[None, "BatchMPC"]:
        """Filter out observations from the batch.

        Parameters
        ----------
        bands : Union[List[str], None], optional
            List of observation bands to keep in the batch, by default None
        catalogs : Union[List[str], None], optional
            List of star catalogs to keep in the batch. See 
            https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
            for the exact encodings used by MPC, by default None
        observation_types : Union[List[str], None], optional
            List of observation types to keep in the batch, e.g. CCD,
            Photographic etc. See the Note2 section of the MPC format description:
            https://www.minorplanetcenter.net/iau/info/OpticalObs.html
            for the exact encoding, by default None
        observatories : Union[List[str], None], optional
            List of observatories to keep in the batch, by default None
        observatories_exclude : Union[List[str], None], optional
            List of observatories to remove from the batch, by default None
        epoch_start : Union[float, datetime.datetime, None], optional
            Start date to include observations from, can be in python datetime in utc\
                 or the more conventional tudat seconds since j2000 in TDB if float,\
                     by default None
        epoch_end : Union[float, datetime.datetime, None], optional
            Final date to include observations from, can be in python datetime in utc\
                 or the more conventional tudat seconds since j2000 in TDB if float,\
                     by default None
        in_place : bool, optional
            If true, modify the current batch object.\
                  If false returns a new object that is\
                      filtered, currect batch remains unmodified, by default True

        Raises
        ------
        ValueError
            Is raised if bands, observatories, or observatories_exclude are not list or 
            None.
        ValueError
            Is raised if both observations_exclude and observatories are not None.

        Returns
        -------
        None or BatchMPC
            Returns a new instance of BatchMPC that is filtered.
        """

        # basic user input handling
        assert not isinstance(
            observatories, int
        ), "stations parameter must be of type 'str' or 'List[str]'"

        if not (isinstance(observatories, list) or (observatories is None)):
            raise ValueError("observatories parameter must be list of strings or None")

        if not (isinstance(catalogs, list) or (catalogs is None)):
            raise ValueError("catalogs parameter must be list of strings or None")

        if not (isinstance(observation_types, list) or (observation_types is None)):
            raise ValueError(
                "observation_types parameter must be list of strings or None"
            )

        if not (
            isinstance(observatories_exclude, list) or (observatories_exclude is None)
        ):
            raise ValueError(
                "observatories_exclude parameter must be list of strings or None"
            )

        if not (isinstance(bands, list) or (bands is None)):
            raise ValueError("bands parameter must be list of strings or None")

        if (observatories is not None) and (observatories_exclude is not None):
            txt = "Include or exclude observatories, not both at the same time."
            raise ValueError(txt)

        if in_place:
            if bands is not None:
                self._table = self._table.query("band == @bands")
            if catalogs is not None:
                self._table = self._table.query("catalog == @catalogs")
            if observation_types is not None:
                self._table = self._table.query("note2 == @observation_types")
            if observatories is not None:
                self._table = self._table.query("observatory == @observatories")
            if observatories_exclude is not None:
                self._table = self._table.query("observatory != @observatories_exclude")
            if epoch_start is not None:
                if isinstance(epoch_start, float) or isinstance(epoch_start, int):
                    self._table = self._table.query(
                        "epochJ2000secondsTDB >= @epoch_start"
                    )
                elif isinstance(epoch_start, datetime.datetime):
                    self._table = self._table.query("epochUTC >= @epoch_start")
            if epoch_end is not None:
                if isinstance(epoch_end, float) or isinstance(epoch_end, int):
                    self._table = self._table.query(
                        "epochJ2000secondsTDB <= @epoch_end"
                    )
                elif isinstance(epoch_end, datetime.datetime):
                    self._table = self._table.query("epochUTC <= @epoch_end")

            self._refresh_metadata()
            return None
        else:
            new = self.copy()
            new.filter(
                bands=bands,
                catalogs=catalogs,
                observation_types=observation_types,
                observatories=observatories,
                observatories_exclude=observatories_exclude,
                epoch_start=epoch_start,
                epoch_end=epoch_end,
                in_place=True,
            )
            return new

    def to_tudat(
        self,
        bodies: environment.SystemOfBodies,
        included_satellites: Union[Dict[str, str], None],
        station_body: str = "Earth",
        add_sbdb_gravity_model: bool = False,
        apply_weights_VFCC17: bool = True,
        apply_star_catalog_debias: bool = True,
        debias_kwargs: dict = dict(),
    ) -> estimation.ObservationCollection:
        """Converts the observations in the batch into a Tudat compatible format and
          sets up the relevant Tudat infrastructure to support estimation.
        This method does the following:\\
            1. (By Default) Applies star catalog debiasing and estimation
            weight calculation.\\
            2. Creates an empty body for each minor planet with their MPC code as a 
            name.\\
            3. Adds this body to the system of bodies inputted to the method.\\
            4. Retrieves the global position of the terrestrial observatories in 
            the batch and adds these stations to the Tudat environment.\\
            5. Creates link definitions between each unique terrestrial 
            observatory/ minor planet combination in the batch.\\
            6. (Optionally) creates a link definition between each 
            space telescope / minor planet combination in the batch. 
            This requires an addional input.\\
            7. Creates a `SingleObservationSet` object for each unique link that 
            includes all observations for that link.\\
            8. (By Default) Add the relevant weights to the `SingleObservationSet`
            per observation.\\
            8. Returns the observations


        Parameters
        ----------
        bodies : environment.SystemOfBodies
            SystemOfBodies containing at least the earth to allow for the placement of 
            terrestrial telescopes
        included_satellites : Union[Dict[str, str], None], optional
            A dictionary that links the name of a space telescope used by the user with 
            the observatory code in MPC. Used when utilising observations from space 
            telescopes like TESS and WISE. The keys should be the MPC observatory 
            codes. The values should be the bodies' in the user's code. The relevant 
            observatory code can be retrieved using 
            the .observatories_table() method, by default None
        station_body : str, optional
            Body to attach ground stations to. Does not need to be changed unless the 
            `Earth` body has been renamed, by default "Earth"
        station_body : bool, optional
            Adds a central_sbdb gravity model to the object, generated using JPL's small body database. 
            This option is only available for a limited number of bodies and raises an error if unavailable. 
            See tudatpy.numerical_simulation.environment_setup.gravity_field.central_sbdb for more info.
            Enabled if True, by default False
        apply_weights_VFCC17 : bool, optional
            Applies the weighting scheme as described in: "Statistical analysis of astrometric errors
            for the most productive asteroid surveys" by Veres et al. (2017). Overwrites custom weights
            set through the `.set_weights()` method if True, by default True
        apply_star_catalog_debias : bool, optional
            Applies star catalog debiasing as described in: "Star catalog position and proper motion corrections 
            in asteroid astrometry II: The Gaia era" by Eggl et al. (2018), by default True
        debias_kwargs : dict, optional
            Additional options when applying star catalog debiasing. A different debias file 
            can be set here. Options are set as kwargs using a dictionary, see data._biases.get_biases_EFCC18() 
            for more info, by default dict()

        Returns
        -------
        estimation.ObservationCollection
            ObservationCollection with the observations found in the batch

        Examples
        ----------
        Using Space Telescopes

        Create dictionary to link name in tudat with mpc code in batch:

        >> sats_dict = {"C57":"TESS"}
        >> obs = batch1.to_tudat(bodies, included_satellites=sats_dict)


        """
        # apply EFCC18 star catalog corrections:
        if apply_star_catalog_debias:
            self._apply_EFCC18(**debias_kwargs)
            RA_col = "RA_EFCC18"
            DEC_col = "DEC_EFCC18"
        else:
            RA_col = "RA"
            DEC_col = "DEC"

        # Calculate observation weights and update table:
        if apply_weights_VFCC17 and not self._custom_weights_set:
            temp_table:pd.DataFrame = get_weights_VFCC17(  # type:ignore
                mpc_table=self.table,
                return_full_table=True,
            )
            self._table.loc[:, "weight"] = temp_table.loc[:, "weight"]
        else:
            temp_table = self.table.copy()

        # start user input validation
        # Ensure that Earth is in the SystemOfBodies object
        try:
            bodies.get(station_body)
        except Exception as e:
            print(
                f"Body {station_body} is not in bodies, if you have renamed Earth, "
                + "set the station_body paramater to the new name."
            )
            raise e

        # Add positions of the observatories
        self._add_observatory_positions(bodies, station_body)

        # get satellites to include and exclude
        if included_satellites is not None:
            sat_obs_codes_included = list(included_satellites.keys())
            # this appears unused but is used in a pandas query:
            sat_obs_codes_excluded = list(
                set(self._space_telescopes) - set(sat_obs_codes_included)
            )

            # Ensure that the satellite is in the SystemOfBodies object
            for sat in list(included_satellites.values()):
                try:
                    bodies.get(sat)
                except Exception as e:
                    print(f"Body {sat} is not in bodies")
                    raise e
        else:
            sat_obs_codes_included = []
            # this appears unused but is used in a pandas query:
            sat_obs_codes_excluded = self._space_telescopes

        # end user input validation

        # get relevant stations positions
        tempStations = self._observatory_info.query("Code == @self.observatories").loc[
            :, ["Code", "X", "Y", "Z"]
        ]

        # add station positions to the observations
        observations_table = pd.merge(
            left=temp_table,  # type: ignore
            right=tempStations,
            left_on="observatory",
            right_on="Code",
            how="left",
        )
        # remove the observations from unused satellite observatories
        observations_table = observations_table.query("Code != @sat_obs_codes_excluded")

        # add asteroid bodies to SystemOfBodies object
        for body in self._MPC_codes:
            if not bodies.does_body_exist(body):
                bodies.create_empty_body(str(body))
                # save the name of the body added
                self._bodies_created[str(body)] = "empty body"

                if add_sbdb_gravity_model:
                    add_gravity_field_model(bodies, str(body), central_sbdb(str(body)))
                    self._bodies_created[str(body)] = "empty body + sbdb gravity"

        # add ground stations to the earth body
        for idx in range(len(tempStations)):
            station_name = tempStations.iloc[idx].Code

            # skip if it is a satellite observatory
            if station_name in self._space_telescopes:
                continue

            ground_station_settings = environment_setup.ground_station.basic_station(
                station_name=station_name,
                station_nominal_position=[
                    tempStations.iloc[idx].X,
                    tempStations.iloc[idx].Y,
                    tempStations.iloc[idx].Z,
                ],
            )

            # Check if station already exists
            if station_name not in bodies.get(station_body).ground_station_list:
                # Add the ground station to the environment
                environment_setup.add_ground_station(
                    bodies.get_body(station_body), ground_station_settings
                )

        # get unique combinations of mpc bodies and observatories
        unique_link_combos = (
            observations_table.loc[:, ["number", "observatory"]].drop_duplicates()
        ).values

        observation_set_list = []
        for combo in unique_link_combos:
            MPC_number = combo[0]
            station_name = combo[1]

            # CREATE LINKS
            link_ends = dict()

            # observed body link
            link_ends[observation.transmitter] = observation.body_origin_link_end_id(
                MPC_number
            )

            if station_name in sat_obs_codes_included:
                # link for a satellite
                sat_name = included_satellites[station_name]
                link_ends[observation.receiver] = observation.body_origin_link_end_id(
                    sat_name
                )
                link_definition = observation.link_definition(link_ends)
            else:
                # link for a ground station
                link_ends[observation.receiver] = (
                    observation.body_reference_point_link_end_id(
                        station_body, station_name
                    )
                )
                link_definition = observation.link_definition(link_ends)

            # get observations, angles and times for this specific link
            observations_for_this_link = observations_table.query(
                "number == @MPC_number"
            ).query("observatory == @station_name")

            observation_angles = observations_for_this_link.loc[
                :, [RA_col, DEC_col]
            ].to_numpy()

            observation_times = observations_for_this_link.loc[
                :, ["epochJ2000secondsTDB"]
            ].to_numpy()[:, 0]

            # create a set of obs for this link
            observation_set = estimation.single_observation_set(
                observation.angular_position_type,
                link_definition,
                observation_angles,
                observation_times,
                observation.receiver,
            )

            # apply weights if apply_weights is True or set_weights() has been used.
            if apply_weights_VFCC17 or self._custom_weights_set:
                observation_weights = observations_for_this_link.loc[
                    :, ["weight"]
                ].to_numpy()[:, 0]
                # this is to make sure the order is RA1, DEC1, 2, 2, 3, 3 etc.
                observation_weights = np.ravel(
                    [observation_weights, observation_weights], "F"
                )

                observation_set.weights_vector = observation_weights

            observation_set_list.append(observation_set)

        observation_collection = estimation.ObservationCollection(observation_set_list)
        return observation_collection

    def plot_observations_temporal(
        self,
        objects: Union[List[str], None] = None,
        figsize: Tuple[float] = (9.0, 6.0),
    ):
        """Generates a matplotlib figure with the declination and right ascension
        over time.


        Parameters
        ----------
        objects : Union[List[str], None], optional
            List of specific MPC objects in batch to plot, None to plot all
            , by default None
        projection : str, optional
            projection of the figure options are: 'aitoff', 'hammer',
            'lambert' and 'mollweide', by default "aitoff"
        figsize : Tuple[float], optional
            size of the matplotlib figure, by default (15, 7)

        Returns
        -------
        Matplotlib figure
            Matplotlib figure with the observations
        """
        fig, (axRA, axDEC) = plt.subplots(2, 1, figsize=figsize, sharex=True)

        if objects is None:
            objs = self.MPC_objects
        else:
            objs = [str(o) for o in objects]

        for obj in objs:
            tab = self._table.query("number == @obj")

            axRA.scatter(
                tab.epochUTC,
                np.degrees(tab.RA),
                label="MPC: " + obj,
                marker=".",
                # linestyle=None,
            )
            axDEC.scatter(
                tab.epochUTC,
                np.degrees(tab.DEC),
                label="MPC: " + obj,
                marker=".",
                # linestyle=None,
            )

        axRA.set_ylabel(r"Right Ascension $[\deg]$")
        axDEC.set_ylabel(r"Declination $[\deg]$")
        axDEC.set_xlabel(r"Time [year]")

        axRA.grid()
        axDEC.grid()
        axRA.set_title("Right Ascension")
        axDEC.set_title("Declination")
        axRA.legend(bbox_to_anchor=(1.01, 1), loc="upper left")

        buffer = 10
        axRA.set_ylim(0 - buffer, 360 + buffer)
        axDEC.set_ylim(-90 - buffer, 90 + buffer)

        fig.set_layout_engine("tight")

        return fig

    def plot_observations_sky(
        self,
        objects: Union[List[str], None] = None,
        projection: Union[str, None] = None,
        figsize: Tuple[float] = (14.0, 7.0),
    ):
        """Generates a matplotlib figure with the observations'
        right ascension and declination over time.

        Parameters
        ----------
        objects : Union[List[str], None], optional
            List of specific MPC objects in batch to plot, None to plot all
            , by default None
        projection : str, optional
            projection of the figure options are: 'aitoff', 'hammer',
            'lambert' and 'mollweide', by default "aitoff"
        figsize : Tuple[float], optional
            size of the matplotlib figure, by default (15, 7)

        Returns
        -------
        Matplotlib figure
            Matplotlib figure with the observations
        """
        fig, ax = plt.subplots(
            1, 1, subplot_kw={"projection": projection}, figsize=figsize
        )

        if objects is None:
            objs = self.MPC_objects
        else:
            objs = [str(o) for o in objects]

        markers = ["o", "+", "^"]

        vmin = mdates.date2num(self._table.query("number == @objs").epochUTC.min())
        vmax = mdates.date2num(self._table.query("number == @objs").epochUTC.max())

        for idx, obj in enumerate(objs):
            tab = self._table.query("number == @obj")

            a = plt.scatter(
                (tab.RA) - np.pi if projection is not None else np.degrees(tab.RA),
                (tab.DEC) if projection is not None else np.degrees(tab.DEC),
                s=30,
                marker=markers[int(idx % len(markers))],
                cmap=cm.plasma,
                label="MPC: " + obj,
                c=mdates.date2num(tab.epochUTC),
                vmin=vmin,
                vmax=vmax,
            )

        ax.legend(markerscale=2)
        ax.set_xlabel(r"Right Ascension $[\deg]$")
        ax.set_ylabel(r"Declination $[\deg]$")

        yticks = [
            f"{x}°"
            for x in (np.degrees(np.array(ax.get_yticks().tolist()))).astype(int)
        ]
        xticks = [
            f"{x}°"
            for x in (np.degrees(np.array(ax.get_xticks().tolist())) + 180).astype(int)
        ]
        if projection in ["aitoff", "hammer", "mollweide"]:
            ax.set_yticklabels(yticks)
            ax.set_xticklabels(xticks)
        elif projection == "lambert":
            ax.set_xticklabels(xticks)
        else:
            pass

        ticks = np.linspace(vmin, vmax, 7)
        labels = [mdates.num2date(t).strftime("%d %b %Y") for t in ticks]
        cbar = plt.colorbar(mappable=a, ax=ax, label="Time", ticks=ticks)

        cbar.ax.set_yticklabels(labels)

        startUTC = self._table.query("number == @objs").epochUTC.min()
        endUTC = self._table.query("number == @objs").epochUTC.max()
        ax.grid()
        fig.suptitle(f"{self.size} observations between {startUTC} and {endUTC}")

        fig.set_layout_engine("tight")

        return fig

    def summary(self):
        """Produce a quick summary of the content of the batch"""
        print()
        print("   Batch Summary:")
        print(f"1. Batch includes {len(self._MPC_codes)} minor planets:")
        print("  ", self.MPC_objects)
        satObs = len(self._table.query("observatory == @self.space_telescopes"))
        print(
            f"2. Batch includes {self.size} observations, including {satObs} "
            + "observations from space telescopes"
        )
        print(
            f"3. The observations range from {self._table.epochUTC.min()} "
            + f"to {self._table.epochUTC.max()}"
        )
        print(f"   In seconds TDB since J2000: {self.epoch_start} to {self.epoch_end}")
        print(
            f"   In Julian Days: "
            + f"{self._table.epoch.min()}"
            + f" to {self._table.epoch.max()}"
        )
        print(
            f"4. The batch contains observations from {len(self.observatories)} "
            + f"observatories, including {len(self.space_telescopes)} space telescopes"
        )
        print()

    def observatories_table(
        self,
        only_in_batch: bool = True,
        only_space_telescopes: bool = False,
        exclude_space_telescopes: bool = False,
        include_positions: bool = False,
    ) -> pd.DataFrame:
        """Returns a pandas DataFrame with information about all MPC observatories,
        Carthesian positions are only available after running the `to_tudat()` method.

        Parameters
        ----------
        only_in_batch : bool, optional
            Filter out observatories that are not in the batch, by default True
        only_space_telescopes : bool, optional
            Filter out all observatories except space telescopes, by default False
        only_space_telescopes : bool, optional
            Filter out all space telescopes, by default False
        include_positions : bool, optional
            Include cartesian positions of the terrestrial telescopes only available
            after running to_tudat(), by default False

        Returns
        -------
        pd.DataFrame
            Dataframe with information about the observatories.
        """
        temp = self._observatory_info
        temp2 = self._table

        count_observations = (
            temp2.groupby("observatory")
            .count()
            .rename(columns={"number": "count"})
            .reset_index(drop=False)
            .loc[:, ["observatory", "count"]]
        )

        temp = pd.merge(
            left=temp,
            right=count_observations,
            left_on="Code",
            right_on="observatory",
            how="left",
        )
        if only_in_batch:
            temp = temp.query("Code == @self.observatories")
        if only_space_telescopes:
            temp = temp.query("Code == @self._MPC_space_telescopes")
        if exclude_space_telescopes:
            temp = temp.query("Code != @self._MPC_space_telescopes")
        if not include_positions:
            temp = temp.loc[:, ["Code", "Name", "count"]]
        return temp
