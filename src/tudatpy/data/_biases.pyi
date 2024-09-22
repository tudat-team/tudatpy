import typing
from astropy.time.core import Time
from astropy import units as u
from astropy.units.quantity import Quantity
from astropy_healpix.high_level import HEALPix
import numpy as np
import os as os
import pandas as pd
import re as re
__all__ = ['BIAS_LOWRES_FILE', 'DEFAULT_CATALOG_FLAGS', 'HEALPix', 'Quantity', 'Time', 'get_biases_EFCC18', 'load_bias_file', 'np', 'os', 'pd', 're', 'u']

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