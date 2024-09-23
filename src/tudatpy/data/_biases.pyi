import typing
from astropy.time.core import Time
from astropy import units as u
from astropy.units.quantity import Quantity
from astropy_healpix.high_level import HEALPix
import astroquery.mpc.core
import numpy as np
import os as os
import pandas as pd
import re as re
__all__ = ['BIAS_LOWRES_FILE', 'DEFAULT_CATALOG_FLAGS', 'HEALPix', 'MPC', 'Quantity', 'Time', 'get_biases_EFCC18', 'get_weights_VFCC17', 'load_bias_file', 'np', 'os', 'pd', 're', 'u']

def get_biases_EFCC18(RA: typing.Union[float, numpy.ndarray, list], DEC: typing.Union[float, numpy.ndarray, list], epochJ2000secondsTDB: typing.Union[float, numpy.ndarray, list], catalog: typing.Union[str, numpy.ndarray, list], bias_file: typing.Optional[str]='/Users/alfonso/.tudat/resource/star_catalog_biases/debias_2018/bias.dat', Nside: typing.Optional[int]=None, catalog_flags: typing.List[str]=['a', 'b', 'c', 'd', 'e', 'g', 'i', 'j', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 't', 'u', 'v', 'w', 'L', 'N', 'Q', 'R', 'S', 'U', 'Y']) -> typing.Tuple[numpy.ndarray, numpy.ndarray]:
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

def get_weights_VFCC17(MPC_codes: typing.Union[pandas.core.series.Series, list, numpy.ndarray, NoneType]=None, epochUTC: typing.Union[pandas.core.series.Series, list, numpy.ndarray, NoneType]=None, observation_type: typing.Union[pandas.core.series.Series, list, numpy.ndarray, NoneType]=None, observatory: typing.Union[pandas.core.series.Series, list, numpy.ndarray, NoneType]=None, star_catalog: typing.Union[pandas.core.series.Series, list, numpy.ndarray, NoneType]=None, mpc_table: typing.Optional[pandas.core.frame.DataFrame]=None, return_full_table=False) -> typing.Union[numpy.ndarray, pandas.core.frame.DataFrame]:
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

def load_bias_file(filepath: str, Nside: typing.Optional[int]=None, catalog_flags: list=['a', 'b', 'c', 'd', 'e', 'g', 'i', 'j', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 't', 'u', 'v', 'w', 'L', 'N', 'Q', 'R', 'S', 'U', 'Y']) -> typing.Tuple[pandas.core.frame.DataFrame, int]:
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
BIAS_LOWRES_FILE: str = '/Users/alfonso/.tudat/resource/star_catalog_biases/debias_2018/bias.dat'
DEFAULT_CATALOG_FLAGS: list = ['a', 'b', 'c', 'd', 'e', 'g', 'i', 'j', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 't', 'u', 'v', 'w', 'L', 'N', 'Q', 'R', 'S', 'U', 'Y']
MPC: astroquery.mpc.core.MPCClass