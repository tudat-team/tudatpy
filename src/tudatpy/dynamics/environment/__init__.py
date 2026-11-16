from tudatpy.kernel.dynamics.environment import *


def __getattr__(name):

    if name == "Tle":
        from tudatpy._deprecation import deprecation_warning
        from tudatpy.dynamics.environment_setup.ephemeris import Tle

        deprecation_warning(
            "tudatpy.dynamics.environment.Tle", "tudatpy.dynamics.environment_setup.ephemeris.Tle"
        )
        return Tle
    raise AttributeError(f"module {__name__} has no attribute {name}")
