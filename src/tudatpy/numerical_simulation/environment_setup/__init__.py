import warnings
warnings.warn(
    "tudatpy.numerical_simulation.environment_setup is deprecated as of v1.0 (see https://docs.tudat.space/en/latest/user-guide/project-updates/migration-guide.html). Use tudatpy.dynamics.environment_setup instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.dynamics.environment_setup import *
