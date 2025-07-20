import warnings
warnings.warn(
    "tudatpy.numerical_simulation.environment_setup is deprecated. Use tudatpy.dynamics.environment_setup instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.dynamics.environment_setup import *
