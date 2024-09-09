import typing
from itertools import combinations as comb
import numpy
import numpy as np
__all__ = ['anova_analysis', 'comb', 'get_orthogonal_array', 'get_yates_array', 'np']

def anova_analysis(objective_values, factorial_design_array: numpy.ndarray, no_of_factors: int, no_of_levels: int, level_of_interactions=2) -> tuple:
    """Function that performs an Analysis of Variance (ANOVA) for no_of_levels levels and no_of_factors factors.
	
	Parameters
	----------
	objective_values : list
		List of objective function values following from the factorial design array.
	factorial_design_array : np.array
		An array compatible with factorial design, generally a Yates array.
	no_of_factors : int
		Integer representing the number of design variables.
	no_of_levels : int
		Integer representing the number of levels that will be explored within each bound.
	level_of_interactions : int
		Integer either equal to 2 or 3, depending on what interactions you want to
		include
	
	Returns
	-------
	percentage_contribution_i : np.ndarray
		Array of percentage contributions of individual effects for each factor
	percentage_contribution_ii : np.ndarray
		Array of percentage contributions of linear effects for each factor
	percentage_contribution_iii : np.ndarray
		Array of percentage contributions of quadratic effects for each factor
	percentage_contribution_error : np.ndarray
		Float for percentage contribution of error
	"""

def get_orthogonal_array(no_of_factors: int, no_of_levels: int) -> numpy.ndarray:
    """Create orthogonal arrays from Latin Square in 4 successive steps:
	
	0) Take the column from the smaller array to create 2 new
	columns and 2x new rows,
	1) block 1 (1/2 rows): take old values 2x for new columns,
	2) block 2 (1/2 rows): take old values, use Latin-Square for new
	columns,
	3) column 1: divide experiments into groups of 1,2.
	
	Parameters
	----------
	no_of_factors : int
		Integer representing the number of design variables.
	no_of_levels : int
		Integer representing the number of levels that will be explored within each bound.
	
	Returns
	-------
	Lx : np.ndarray
		Orthogonal array containing the FFD experiments.
	"""

def get_yates_array(no_of_factors: int, no_of_levels: int) -> numpy.ndarray:
    """Function that creates a Yates array according to Yates algorithm
	
	Parameters
	----------
	no_of_factors : int
		Integer representing the number of design variables.
	no_of_levels : int
		Integer representing the number of levels that will be explored within each bound.
	
	Returns
	-------
	yates_array : np.ndarray
		Orthogonal array containing the FFD experiments.
	"""