import warnings
warnings.warn(
    "tudatpy.numerical_simulation.environment is deprecated. Use tudatpy.dynamics.environment instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.dynamics.environment import *