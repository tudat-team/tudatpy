"""
This module contains functions and classes to read atmospheric S1/S2 tidal loading coefficients, either on a global
grid (for instance the Ray and Ponte, 2003, model used in the IERS Conventions 2010, Section 7.1.3) or for
individual stations.

The data is read without interpolation or conversion of its reference origin. It is used to create the atmospheric
loading deformation model with the functions
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_grid` and
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_stations`.
"""

from tudatpy.kernel.data_input.environment_data.atmospheric_loading import *
