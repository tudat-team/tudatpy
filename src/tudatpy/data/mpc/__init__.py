from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "BatchMPC": "tudatpy.data_input.tracking_data.mpc.BatchMPC",
    "load_bias_file": "tudatpy.data_input.tracking_data.optical_utilities.load_bias_file",
    "get_biases_EFCC18": "tudatpy.data_input.tracking_data.optical_utilities.get_biases_EFCC18",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
