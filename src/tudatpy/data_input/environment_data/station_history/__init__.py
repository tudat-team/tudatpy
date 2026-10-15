"""
This module contains functions and classes to read SINEX station files, such as the ITRF and IGS station coordinate
solutions and their extensions (discontinuities, post-seismic deformation, seasonal terms, and receiver, antenna and
eccentricity histories).

The records are read without creating a station model. They are used to create station motion models with the
functions in :ref:`ground_station` (for instance
:func:`~tudatpy.dynamics.environment_setup.ground_station.ign_station_history`).
For one position and velocity per DOMES monument, the simpler
:func:`~tudatpy.data_input.environment_data.ilrs.read_sinex_station_data` function can be used.
"""

from tudatpy.kernel.data_input.environment_data.station_history import *
