import importlib
import warnings


def deprecated_getattr(module_name, aliases, name):
    if name not in aliases:
        raise AttributeError(f"module {module_name!r} has no attribute {name!r}")

    target_module_name, target_name = aliases[name].rsplit(".", 1)
    warnings.warn(
        f"{module_name}.{name} is deprecated. Use {target_module_name}.{target_name} instead.",
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
