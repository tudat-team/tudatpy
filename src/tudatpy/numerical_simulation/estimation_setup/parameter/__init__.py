import warnings
warnings.warn(
    "tudatpy.numerical_simulation.estimation_setup.parameter is deprecated. Use tudatpy.dynamics.parameters_setup instead.",
    stacklevel=1
)

from tudatpy.kernel.dynamics.parameters_setup import *
