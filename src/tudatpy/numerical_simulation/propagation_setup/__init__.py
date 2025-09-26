import warnings
warnings.warn(
    "tudatpy.numerical_simulation.propagation_setup is deprecated as of v1.0 (see https://docs.tudat.space/en/latest/user-guide/project-updates/migration-guide.html). Use tudatpy.dynamics.propagation_setup instead.",
    FutureWarning,
    stacklevel=1
)
from tudatpy.dynamics.propagation_setup import *