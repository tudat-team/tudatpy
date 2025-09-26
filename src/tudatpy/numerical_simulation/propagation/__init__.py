import warnings
warnings.warn(
    "tudatpy.numerical_simulation.propagation is deprecated as of v1.0 (see https://docs.tudat.space/en/latest/user-guide/project-updates/migration-guide.html). Use tudatpy.dynamics.propagation instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.dynamics.propagation import *