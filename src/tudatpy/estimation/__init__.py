from importlib import import_module

from tudatpy.kernel.estimation import *

# Select the Python packages so fresh parent imports retain their compatibility
# bridges instead of exposing the raw kernel modules.
observations = import_module(__name__ + ".observations")
observations_setup = import_module(__name__ + ".observations_setup")
