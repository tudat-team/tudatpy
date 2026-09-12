from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "SpaceTrackQuery": "tudatpy.data_input.environment_data.spacetrack.SpaceTrackQuery",
    "OMMUtils": "tudatpy.data_input.environment_data.spacetrack.OMMUtils",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
