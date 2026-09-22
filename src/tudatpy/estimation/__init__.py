from tudatpy.kernel.estimation import *

# Select the Python package so fresh parent imports retain its legacy wrapper.
from importlib import import_module as _import_module

observations_setup = _import_module("tudatpy.estimation.observations_setup")
del _import_module
