"""Deprecated compatibility wrapper for observation I/O and simulation helpers."""

from importlib import import_module
import warnings

_TARGET_MODULE_NAME = "tudatpy.estimation.observations"
_target_module = import_module(_TARGET_MODULE_NAME)

__all__ = [name for name in dir(_target_module) if not name.startswith("_")]


def __getattr__(name):
    if not hasattr(_target_module, name):
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    warnings.warn(
        f"{__name__}.{name} is deprecated. Use {_TARGET_MODULE_NAME}.{name} instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return getattr(_target_module, name)


def __dir__():
    return sorted(set(globals()) | set(__all__))
