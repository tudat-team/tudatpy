from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "RadioBase": "tudatpy.data_input.tracking_data.tnf._converters.RadioBase",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
