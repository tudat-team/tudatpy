import importlib as _importlib

from tudatpy.kernel.estimation import *

# The wildcard import above exposes the kernel ``observations`` submodule as an
# attribute. Replace it with TudatPy's Python package so its helper functions are
# available through ``from tudatpy.estimation import observations``.
observations = _importlib.import_module("tudatpy.estimation.observations")

del _importlib
