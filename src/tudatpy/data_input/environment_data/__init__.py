"""Environment-data readers and query interfaces.

This module groups helpers for external environment data sources, including
SPICE kernels, Horizons, SBDB, Space-Track, DIScos, ILRS/SINEX data, space
weather, comet-coma data, and mission-specific environment products.
"""

from tudatpy.kernel.data_input.environment_data import *

from . import coma, discos, horizons, ilrs, missions, sbdb, space_weather, spacetrack, spice
