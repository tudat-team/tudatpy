from tudatpy.kernel.estimation.observations_setup import *

# Select the Python package, which provides deprecated aliases on top of the kernel submodule: importing it binds it
# as attribute of this module, replacing the kernel submodule bound by the star import above
import tudatpy.estimation.observations_setup.ancillary_settings


def __getattr__(name):
    # Removable compatibility bridge for the former kernel submodule. Load it
    # lazily to avoid import cycles; its legacy function access emits the warning.
    if name == "observations_wrapper":
        from importlib import import_module

        return import_module(f"{__name__}.observations_wrapper")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | {"observations_wrapper"})
