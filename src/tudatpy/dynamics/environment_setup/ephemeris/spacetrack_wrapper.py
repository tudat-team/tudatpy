from tudatpy.dynamics import environment


def tle_to_tle_ephemeris(tle_line_1: str, tle_line_2: str):
    """Convert a TLE line pair into a Tudat TleEphemeris object."""
    return environment.TleEphemeris(
        "Earth", "J2000", environment.Tle(tle_line_1, tle_line_2), False
    )


def tle_to_tle(tle_line_1: str, tle_line_2: str):
    """Convert a TLE line pair into a Tudat Tle object."""
    return environment.Tle(tle_line_1, tle_line_2)
