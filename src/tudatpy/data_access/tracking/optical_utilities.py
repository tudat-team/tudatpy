import pandas as pd
import numpy as np
from astropy_healpix import HEALPix
from astropy.units import Quantity
import astropy.units as u
import os
import re
from tudatpy.astro import time_representation

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
    Nside: int | None = None,
    catalog_flags: list = DEFAULT_CATALOG_FLAGS,
) -> tuple[pd.DataFrame, int]:
    """Loads a healpix star catalog debias file and processes it into a dataframe. Automatically retrieves NSIDE parameter.

    Parameters
    ----------
    filepath : str
        Filepath of debias file.
    Nside : int | None, optional
        NSIDE value, to be left None in most cases as this is retrieved automatically by the function, by default None
    catalog_flags : list | None, optional
        list of catalog flags, should be left default in most cases, by default None

    Returns
    -------
    tuple[pd.DataFrame, int]
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
    bias_dataframe = bias_dataframe.stack(level=0, future_stack=True)

    return bias_dataframe, Nside


def get_biases_EFCC18(
    mpc_table: pd.DataFrame,
    bias_file: str | None = BIAS_LOWRES_FILE,
    Nside: int | None = None,
    catalog_flags: list[str] = DEFAULT_CATALOG_FLAGS,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate and return star catalog bias values as described in:
    "Star catalog position and proper motion corrections in asteroid astrometry II: The Gaia era" by Eggl et al. (2018).
    Uses the regular bias set by default. A high res version of the bias map can be retrieved from the paper.
    This can then be selected using the bias_file paramater.

    Parameters
    ----------
    mpc_table : pd.DataFrame
        Table retrieved by calling the mpc.BatchMPC.table property. Must contain
        'RA', 'DEC', 'epoch_seconds_UTC' and 'catalog' columns.
    bias_file : str | None, optional
        Optional bias file location, used to load in alternative debias coefficients. By default coefficients are retrieved from Tudat resources, by default None
    Nside : int | None, optional
        Optional Nside value, should be left None in most cases, by default None
    catalog_flags : list[str] | None, optional
        List of catalog values to use, should be left None in most cases, by default None

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Right Ascencion Corrections, Declination corrections
    """

    if bias_file is None:
        bias_file = BIAS_LOWRES_FILE

    RA = mpc_table["RA"].to_numpy()
    DEC = mpc_table["DEC"].to_numpy()
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


def create_augmented_optical_table(*args, **kwargs):
    from tudatpy.data_access.tracking.mpc import (
        create_augmented_optical_table as _create_augmented_optical_table,
    )

    return _create_augmented_optical_table(*args, **kwargs)


def optical_table_to_tracking_data(*args, **kwargs):
    from tudatpy.data_access.tracking.mpc import (
        optical_table_to_tracking_data as _optical_table_to_tracking_data,
    )

    return _optical_table_to_tracking_data(*args, **kwargs)


def read_astropy_optical_data(*args, **kwargs):
    from tudatpy.data_access.tracking.mpc import (
        read_astropy_optical_data as _read_astropy_optical_data,
    )

    return _read_astropy_optical_data(*args, **kwargs)


def read_pandas_optical_data(*args, **kwargs):
    from tudatpy.data_access.tracking.mpc import (
        read_pandas_optical_data as _read_pandas_optical_data,
    )

    return _read_pandas_optical_data(*args, **kwargs)


__all__ = [
    "BIAS_LOWRES_FILE",
    "DEFAULT_CATALOG_FLAGS",
    "create_augmented_optical_table",
    "get_biases_EFCC18",
    "load_bias_file",
    "optical_table_to_tracking_data",
    "read_astropy_optical_data",
    "read_pandas_optical_data",
]
