import warnings
warnings.warn(
    "tudatpy.numerical_simulation.estimation is deprecated as of v1.0 (see https://docs.tudat.space/en/latest/user-guide/project-updates/migration-guide.html). Use tudatpy.estimation.observations, tudatpy.estimation.observations_setup and/or tudatpy.estimation.estimation_analysis instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.estimation.observations import *
from tudatpy.kernel.estimation.observations_setup import *
from tudatpy.kernel.estimation.estimation_analysis import *

from tudatpy.kernel.estimation.observations_setup.ancillary_settings import *
from tudatpy.kernel.estimation.observations_setup.observations_dependent_variables import *
from tudatpy.kernel.estimation.observations_setup.observations_simulation_settings import *
from tudatpy.kernel.estimation.observations_setup.observations_wrapper import *
from tudatpy.kernel.estimation.observations_setup.random_noise import *
from tudatpy.kernel.estimation.observations_setup.viability import *

from tudatpy.kernel.estimation.observations.observations_geometry import * 
from tudatpy.kernel.estimation.observations.observations_processing import * 



