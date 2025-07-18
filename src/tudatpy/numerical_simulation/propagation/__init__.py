import warnings
warnings.warn(
    "tudatpy.numerical_simulation.propagation is deprecated. Use tudatpy.dynamics.propagation instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.dynamics.propagation import *