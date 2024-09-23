import typing
'\n\nCopyright (c) 2010-2023, Delft University of Technology\nAll rigths reserved\n\nThis file is part of the Tudat. Redistribution and use in source and\nbinary forms, with or without modification, are permitted exclusively\nunder the terms of the Modified BSD license. You should have received\na copy of the license with this file. If not, please or visit:\nhttp://tudat.tudelft.nl/LICENSE.\n'
import matplotlib as matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FuncFormatter
import numpy
import numpy as np
from scipy import ndimage
from tudatpy.astro import time_conversion
from tudatpy import constants
__all__ = ['AutoMinorLocator', 'FuncFormatter', 'LinearSegmentedColormap', 'constants', 'matplotlib', 'ndimage', 'np', 'plot_porkchop', 'plot_porkchop_of_single_field', 'plt', 'time_conversion']

def plot_porkchop(departure_body: str, target_body: str, departure_epochs: numpy.ndarray, arrival_epochs: numpy.ndarray, delta_v: numpy.ndarray, C3: bool=False, total: bool=False, threshold: float=10, upscale: bool=False, number_of_levels: int=10, percent_margin: float=5, figsize: tuple[int, int]=(8, 8), show: bool=True, save: bool=False, filename: str='porkchop.png') -> None:
    """ΔV/C3 porkchop mission design plot.
	
	Parameters
	----------
	departure_body: str
		The name of the body from which the transfer is to be computed
	target_body: str
		The name of the body to which the transfer is to be computed
	departure_epochs: np.ndarray
		Discretized departure time window
	arrival_epochs: np.ndarray
		Discretized arrival time window
	delta_v: np.ndarray
		Array containing the ΔV of all coordinates of the grid of departure/arrival epochs
	C3: bool = False
		Whether to plot C3 (specific energy) instead of ΔV
	total: bool = False
		Whether to plot departure and arrival ΔV/C3, or only the total ΔV/C3. This option is only respected if the ΔV map obtained from
	threshold: float = 10
		Upper threshold beyond which ΔV/C3 is not plotted. This is useful to mask regions of the plot where the ΔV/C3 is too high to be of interest.
	upscale: bool = False
		Whether to use interpolation to increase the resolution of the plot. This is not always reliable, and the detail generated cannot be relied upon for analysis. Its only purpose is aesthetic improvement.
	number_of_levels: int = 10
		The number of levels in the ΔV/C3 contour plot
	percent_margin: float = 5
		Empty margin between the axes of the plot and the plotted data
	figsize: tuple[int, int] = (8, 8)
		Size of the figure
	show: bool bool = True
		Whether to show the plot
	save: bool = False
		Whether to save the plot
	filename: str = 'porkchop.png'
		The filename used for the saved plot
	"""

def plot_porkchop_of_single_field(departure_body: str, target_body: str, departure_epochs: numpy.ndarray, arrival_epochs: numpy.ndarray, delta_v: numpy.ndarray, C3: bool=False, threshold: float=10, upscale: bool=False, filled_contours: bool=True, plot_out_of_range_field: bool=True, plot_isochrones: bool=True, plot_global_optimum: bool=True, plot_minor_ticks: bool=True, number_of_levels: int=10, line_width: float=0.5, line_color: str='black', font_size: float=7, label: str=False, colorbar: bool=False, percent_margin: float=5, fig: matplotlib.figure.Figure=None, ax: matplotlib.axes._axes.Axes=None, figsize: tuple[int, int]=(8, 8), show: bool=True, save: bool=False, filename: str='porkchop.png') -> matplotlib.contour.QuadContourSet:
    """Create a ΔV/C3 porkchop mission design plot of single time window-ΔV field.
	
	Parameters
	----------
	departure_body: str
		The name of the body from which the transfer is to be computed
	target_body: str
		The name of the body to which the transfer is to be computed
	departure_epochs: np.ndarray
		Discretized departure time window
	arrival_epochs: np.ndarray
		Discretized arrival time window
	delta_v: np.ndarray
		Array containing the ΔV of all coordinates of the grid of departure/arrival epochs
	C3: bool = False
		Whether to plot C3 (specific energy) instead of ΔV
	threshold: float = 10
		Upper threshold beyond which ΔV/C3 is not plotted. This is useful to mask regions of the plot where the ΔV/C3 is too high to be of interest.
	upscale: bool = False
		Whether to use interpolation to increase the resolution of the plot. This is not always reliable, and the detail generated cannot be relied upon for analysis. Its only purpose is aesthetic improvement.
	filled_contours: bool = True
		Whether to plot filled contours or else just the contour lines
	plot_out_of_range_field: bool = Tru
		Whether to plot the out-of-range field (ΔV/C3 above the threshold) in a different color
	plot_isochrones: bool = True
		Whether to plot the isochrone lines (constant time of flight) on the plot
	plot_global_optimum: bool = True
		Whether to mark the global optimum with a cross
	plot_minor_ticks: bool = True
		Whether to show minor ticks on the axes
	number_of_levels: int = 10
		The number of levels in the ΔV/C3 contour plot
	line_width: float = 0.5
		Width of the contour plot lines
	line_color: str = 'black'
		Color of the contour plot lines
	font_size: float = 7
		Font size of the contour plot labels
	label: str = False
		Label used to identify the contour plot in legends
	colorbar: bool = False
		Whether to plot a colorbar
	percent_margin: float = 5
		Empty margin between the axes of the plot and the plotted data
	fig: matplotlib.figure.Figure = None
		Figure on which to plot
	ax:  matplotlib.axes._axes.Axes = None
		Axis on which to plot
	figsize: tuple[int, int] = (8, 8)
		Size of the figure
	show: bool bool = True
		Whether to show the plot
	save: bool = False
		Whether to save the plot
	filename: str = 'porkchop.png'
		The filename used for the saved plot
	
	Output
	------
	contour_lines: matplotlib.contour.QuadContourSet
		The contour plot Matplotlib object
	"""