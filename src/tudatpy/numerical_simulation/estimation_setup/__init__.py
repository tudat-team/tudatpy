from tudatpy.kernel.numerical_simulation.estimation_setup import *
# this is kept for now because of a few remaining expose functions

import warnings
warnings.warn(
    "tudatpy.numerical_simulation.estimation_setup is deprecated. " \
    "Use tudatpy.estimation_refactoring.observable_models_setup, tudatpy.estimation_refactoring.observable_models, "
    " tudatpy.estimation_refactoring.observations_setup, tudatpy.estimation_refactoring.observations, " \
    "and/or tudatpy.estimation_refactoring.estimation_analysis instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.estimation_refactoring.observable_models_setup import *
from tudatpy.kernel.estimation_refactoring.observable_models import *
from tudatpy.kernel.estimation_refactoring.observations_setup import *
from tudatpy.kernel.estimation_refactoring.observations import *
from tudatpy.kernel.estimation_refactoring.estimation_analysis import *
