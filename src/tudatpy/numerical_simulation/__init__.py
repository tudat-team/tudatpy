import warnings
warnings.warn(
    "tudatpy.numerical_simulation is deprecated as of v1.0 (see https://docs.tudat.space/en/latest/user-guide/project-updates/migration-guide.html). Use tudatpy.dynamics and/or tudatpy.estimation.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.dynamics import *
from tudatpy.kernel.estimation import *

from tudatpy.kernel.dynamics.simulator import * 
from tudatpy.kernel.estimation.estimation_analysis import * 
