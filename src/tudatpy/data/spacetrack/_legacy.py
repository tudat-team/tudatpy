"""Environment conversions retained for the deprecated Space-Track API."""

from tudatpy.data_input.environment_data.spacetrack import OMMUtils as _OMMUtils
from tudatpy.dynamics.environment_setup.ephemeris.spacetrack_wrapper import (
    tle_to_tle,
    tle_to_tle_ephemeris,
)


class OMMUtils(_OMMUtils):
    tle_to_Tle_object = staticmethod(tle_to_tle)
    tle_to_TleEphemeris_object = staticmethod(tle_to_tle_ephemeris)
