"""
This module contains functions and classes to read Earth orientation parameters (polar motion, UT1-UTC, length of day
and celestial pole offsets) from IERS C04 and ``finals2000A`` files.

The data is read without interpolation, and can be supplied to the settings of the high-accuracy Earth rotation model
(see :func:`~tudatpy.dynamics.environment_setup.rotation_model.gcrs_to_itrs`), which interpolates it.
"""

from tudatpy.kernel.data_input.environment_data.earth_orientation import *
