from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "HorizonsBatch": "tudatpy.data_input.environment_data.horizons.HorizonsBatch",
    "HorizonsQuery": "tudatpy.data_input.environment_data.horizons.HorizonsQuery",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
