from tudatpy.dynamics import environment
from tudatpy.kernel.dynamics.environment_setup.ephemeris import Tle


def tle_to_tle_ephemeris(tle_line_1: str, tle_line_2: str):
    """Convert a TLE line pair into a Tudat TleEphemeris object."""
    return environment.TleEphemeris("Earth", "J2000", Tle(tle_line_1, tle_line_2), False)


def tle_to_tle(tle_line_1: str, tle_line_2: str):
    """Convert a TLE line pair into a Tudat Tle object."""
    return Tle(tle_line_1, tle_line_2)
