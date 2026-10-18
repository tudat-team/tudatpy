"""
This module contains functions and classes to read the ocean pole tide loading grid of the IERS Conventions 2010
(Section 7.1.5).

The data is read without interpolation or conversion of its reference origin. It is used to create the corresponding
deformation model with the function
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.ocean_pole_tidal_from_data`. Ocean tidal loading
coefficients (BLQ files) are read by :func:`~tudatpy.dynamics.environment_setup.shape_deformation.ocean_tidal`.
"""

from tudatpy.kernel.data_input.environment_data.ocean_loading import *
