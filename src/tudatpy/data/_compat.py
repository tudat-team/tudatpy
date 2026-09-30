import importlib
import sys
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


def _warn_deprecated(message):
    # Import machinery and nested legacy adapters must not hide the user's
    # call site. Keep this policy inside the removable compatibility package.
    frame = sys._getframe()
    while frame is not None:
        module = frame.f_globals.get("__name__", "")
        if not (
            module == "tudatpy.data"
            or module.startswith("tudatpy.data.")
            or module in ("_frozen_importlib", "_frozen_importlib_external")
            or module.startswith("importlib.")
        ):
            break
        frame = frame.f_back
    try:
        if frame is None:
            warnings.warn(message, DeprecationWarning, stacklevel=2)
        else:
            # Explicit location avoids double-counting import frames that
            # warnings.warn already skips on some Python versions. Reuse the
            # caller's registry to retain normal filtering and deduplication.
            warnings.warn_explicit(
                message,
                DeprecationWarning,
                frame.f_code.co_filename,
                frame.f_lineno,
                module=frame.f_globals.get("__name__", "<string>"),
                registry=frame.f_globals.setdefault("__warningregistry__", {}),
                module_globals=frame.f_globals,
            )
    finally:
        del frame


def deprecated_getattr(module_name, aliases, name):
    if name not in aliases:
        raise AttributeError(f"module {module_name!r} has no attribute {name!r}")

    target_module_name, target_name = aliases[name].rsplit(".", 1)
    migration_target = _MIGRATION_TARGETS.get(aliases[name], aliases[name])
    _warn_deprecated(
        f"{module_name}.{name} is deprecated. Use {migration_target} instead.",
    )
    return getattr(importlib.import_module(target_module_name), target_name)


def warn_custom_deprecation(module_name, name, message):
    _warn_deprecated(f"{module_name}.{name} is deprecated. {message}")


def deprecated_dir(module_globals, aliases):
    return sorted(set(module_globals) | set(aliases))
