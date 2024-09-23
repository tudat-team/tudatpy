import os
import re

import numpy as np
import pandas as pd

from typing import Union, Tuple, List

from astropy_healpix import HEALPix
from astropy.units import Quantity
import astropy.units as u
from astropy.time import Time

from astroquery.mpc import MPC

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
