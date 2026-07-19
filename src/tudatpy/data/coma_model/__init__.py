from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "ComaPolyDataset": "tudatpy.data_input.environment_data.coma.ComaPolyDataset",
    "ComaStokesDataset": "tudatpy.data_input.environment_data.coma.ComaStokesDataset",
    "ComaWindDatasetCollection": "tudatpy.data_input.environment_data.coma.ComaWindDatasetCollection",
    "ComaModelFileProcessor": "tudatpy.data_input.environment_data.coma.ComaModelFileProcessor",
    "ComaWindModelFileProcessor": "tudatpy.data_input.environment_data.coma.ComaWindModelFileProcessor",
    "coma_model_file_processor": "tudatpy.data_input.environment_data.coma.coma_model_file_processor",
    "coma_wind_file_processor": "tudatpy.data_input.environment_data.coma.coma_wind_file_processor",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
