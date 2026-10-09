from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "get_first_failure_reason": "tudatpy.data_input.tracking_data.obs_80_cols.parsers.get_first_failure_reason",
    "parse_80cols_data": "tudatpy.data_input.tracking_data.obs_80_cols.parsers.parse_80cols_data",
    "parse_80cols_file": "tudatpy.data_input.tracking_data.obs_80_cols.parsers.parse_80cols_file",
    "identify_object": "tudatpy.data_input.tracking_data.obs_80_cols.parsers.identify_object",
    "enrich_observations": "tudatpy.data_input.tracking_data.obs_80_cols.parsers.enrich_observations",
    "parse_packed_permanent_designation": "tudatpy.data_input.tracking_data.obs_80_cols.parsers.parse_packed_permanent_designation",
    "parse_80cols_identification_fields": "tudatpy.data_input.tracking_data.obs_80_cols.parsers.parse_80cols_identification_fields",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
