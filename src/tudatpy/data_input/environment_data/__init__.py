"""Environment-data readers and query interfaces.

This module groups helpers for external environment data sources, including
SPICE kernels, Horizons, SBDB, Space-Track, DIScos, ILRS/SINEX data, space
weather, Earth orientation records, comet-coma data, and mission-specific environment products.
"""

from tudatpy.kernel.data_input.environment_data import *

from . import (
    atmospheric_loading,
    coma,
    discos,
    earth_orientation,
    horizons,
    ilrs,
    missions,
    ocean_loading,
    sbdb,
    space_weather,
    spacetrack,
    spice,
    terrestrial_frames,
)
