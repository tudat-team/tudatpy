"""
This module contains functions and classes to read precomputed station loading series (time series of station
displacements due to, for instance, non-tidal atmospheric, oceanic or hydrological loading), such as the EPHEDISP files
of the International Mass Loading Service.

The data is read without interpolation, station matching or conversion of its reference origin. It is used to create
the station loading deformation model with the function
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.station_loading`.
"""

from tudatpy.kernel.data_input.environment_data.station_loading import *
