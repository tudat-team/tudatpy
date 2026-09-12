import importlib
import warnings

_MIGRATION_TARGETS = {
    "tudatpy.data.mpc._legacy.BatchMPC": "tudatpy.data_input.tracking_data.mpc.BatchMPC",
    "tudatpy.data.mpc._legacy.get_biases_EFCC18": "tudatpy.data_input.tracking_data.optical_utilities.get_biases_EFCC18",
    "tudatpy.data.processTrk234._legacy.Trk234Processor": "tudatpy.data_input.tracking_data.tnf.read_tnf_data",
    "tudatpy.data.spacetrack._legacy.OMMUtils": "tudatpy.data_input.environment_data.spacetrack.OMMUtils",
    "tudatpy.dynamics.environment_setup.ephemeris.horizons_wrapper.HorizonsQuery": "tudatpy.data_input.environment_data.horizons.HorizonsQuery",
    "tudatpy.dynamics.environment_setup.ephemeris.horizons_wrapper.HorizonsBatch": "tudatpy.data_input.environment_data.horizons.HorizonsBatch",
}


for _converter in ("RadioBase", "DerivedDopplerConverter", "DerivedSraRangeConverter"):
    _MIGRATION_TARGETS[f"tudatpy.data.processTrk234._legacy_converters.{_converter}"] = (
        "tudatpy.data_input.tracking_data.tnf.read_tnf_data"
    )


def deprecated_getattr(module_name, aliases, name):
    if name not in aliases:
        raise AttributeError(f"module {module_name!r} has no attribute {name!r}")

    target_module_name, target_name = aliases[name].rsplit(".", 1)
    migration_target = _MIGRATION_TARGETS.get(aliases[name], aliases[name])
    warnings.warn(
        f"{module_name}.{name} is deprecated. Use {migration_target} instead.",
        DeprecationWarning,
        stacklevel=3,
    )
    return getattr(importlib.import_module(target_module_name), target_name)


def warn_custom_deprecation(module_name, name, message):
    warnings.warn(
        f"{module_name}.{name} is deprecated. {message}",
        DeprecationWarning,
        stacklevel=4,
    )


def deprecated_dir(module_globals, aliases):
    return sorted(set(module_globals) | set(aliases))
