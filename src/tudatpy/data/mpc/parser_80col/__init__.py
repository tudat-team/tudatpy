from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "parse_80cols_file": "tudatpy.data_input.tracking_data.obs_80_cols.parse_80cols_file",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
