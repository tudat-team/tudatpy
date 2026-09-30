from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "unpack_permanent_minor_planet": "tudatpy.data_input.tracking_data.obs_80_cols.unpackers.unpack_permanent_minor_planet",
    "unpack_provisional_minor_planet": "tudatpy.data_input.tracking_data.obs_80_cols.unpackers.unpack_provisional_minor_planet",
    "unpack_provisional_comet_or_satellite": "tudatpy.data_input.tracking_data.obs_80_cols.unpackers.unpack_provisional_comet_or_satellite",
    "unpack_permanent_natural_satellite": "tudatpy.data_input.tracking_data.obs_80_cols.unpackers.unpack_permanent_natural_satellite",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
