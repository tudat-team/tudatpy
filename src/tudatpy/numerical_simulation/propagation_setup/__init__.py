import warnings
warnings.warn(
    "tudatpy.numerical_simulation.propagation_setup is deprecated. Use tudatpy.dynamics.propagation_setup instead.",
    FutureWarning,
    stacklevel=1
)
from tudatpy.dynamics.propagation_setup import *