import warnings
warnings.warn(
    "tudatpy.numerical_simulation.estimation_setup.observation is deprecated. Use tudatpy.estimation.observable_models_setup and/or tudatpy.estimation.observations_setup instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.estimation.observable_models_setup import *
from tudatpy.kernel.estimation.observable_models_setup.biases import *
from tudatpy.kernel.estimation.observable_models_setup.light_time_corrections import *
from tudatpy.kernel.estimation.observable_models_setup.links import *
from tudatpy.kernel.estimation.observable_models_setup.model_settings import *

from tudatpy.kernel.estimation.observations_setup import * 
from tudatpy.kernel.estimation.observations_setup.ancillary_settings import *
from tudatpy.kernel.estimation.observations_setup.observations_dependent_variables import *
from tudatpy.kernel.estimation.observations_setup.observations_simulation_settings import *
from tudatpy.kernel.estimation.observations_setup.observations_wrapper import *
from tudatpy.kernel.estimation.observations_setup.random_noise import *
from tudatpy.kernel.estimation.observations_setup.viability import *

