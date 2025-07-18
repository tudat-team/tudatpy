import warnings
warnings.warn(
    "tudatpy.numerical_simulation.estimation_setup is deprecated. " \
    "Features got distributed over tudatpy.estimation.observable_models_setup, " \
    " tudatpy.estimation.observable_models, "
    " tudatpy.estimation.observations_setup, " \
    " tudatpy.estimation.observations, " \
    " tudatpy.dynamics.parameters_setup, " \
    " tudatpy.dynamics.parameters, " \
    " and tudatpy.estimation.estimation_analysis instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.estimation.observable_models_setup import *
from tudatpy.kernel.estimation.observable_models import *
from tudatpy.kernel.estimation.observations_setup import *
from tudatpy.kernel.estimation.observations import *
from tudatpy.kernel.dynamics.parameters_setup import *
from tudatpy.kernel.dynamics.parameters import *
from tudatpy.kernel.estimation.estimation_analysis import *

from tudatpy.kernel.estimation.observable_models_setup.biases import *
from tudatpy.kernel.estimation.observable_models_setup.light_time_corrections import *
from tudatpy.kernel.estimation.observable_models_setup.links import *
from tudatpy.kernel.estimation.observable_models_setup.model_settings import *

from tudatpy.kernel.estimation.observations_setup.ancillary_settings import *
from tudatpy.kernel.estimation.observations_setup.observations_dependent_variables import *
from tudatpy.kernel.estimation.observations_setup.observations_simulation_settings import *
from tudatpy.kernel.estimation.observations_setup.observations_wrapper import *
from tudatpy.kernel.estimation.observations_setup.random_noise import *
from tudatpy.kernel.estimation.observations_setup.viability import *

from tudatpy.kernel.estimation.observable_models.observables_simulation import * 

from tudatpy.kernel.estimation.observations.observations_geometry import * 
from tudatpy.kernel.estimation.observations.observations_processing import * 

from tudatpy.kernel.estimation.estimation_analysis import * 

from tudatpy.kernel.astro.time_representation import * 

import sys
import importlib
def __getattr__(name):

    if name == "parameter":
        mod = importlib.import_module("tudatpy.numerical_simulation.estimation_setup.parameter")
        sys.modules["tudatpy.numerical_simulation.estimation_setup.parameter"] = mod
        return mod
    
    if name == "observation":
        mod = importlib.import_module("tudatpy.numerical_simulation.estimation_setup.observation")
        sys.modules["tudatpy.numerical_simulation.estimation_setup.observation"] = mod
        return mod
    
    raise AttributeError(f"module {__name__} has no attribute {name}")
