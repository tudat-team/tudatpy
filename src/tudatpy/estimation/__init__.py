from importlib import import_module

from tudatpy.kernel.estimation import *

observations = import_module(__name__ + ".observations")
