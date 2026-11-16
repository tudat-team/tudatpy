from importlib import import_module
import warnings

from tudatpy.kernel.estimation.observations_setup.ancillary_settings import *


def __getattr__(name):
    # Moved to tudatpy.dynamics.environment, kept here as deprecated aliases
    if name not in ("FrequencyBands", "dsn_default_turnaround_ratios", "cassini_turnaround_ratios"):
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    warnings.warn(
        f"{__name__}.{name} is deprecated. Use tudatpy.dynamics.environment.{name} instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return getattr(import_module("tudatpy.dynamics.environment"), name)


def __dir__():
    return sorted(
        set(globals())
        | {"FrequencyBands", "dsn_default_turnaround_ratios", "cassini_turnaround_ratios"}
    )
