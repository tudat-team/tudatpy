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
