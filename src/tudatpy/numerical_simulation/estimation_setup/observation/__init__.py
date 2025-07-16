import warnings
warnings.warn(
    "tudatpy.numerical_simulation.estimation_setup.observation is deprecated. " \
    "Use tudatpy.estimation_refactoring.observable_models_setup and/or tudatpy.estimation_refactoring.observations_setup instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.estimation_refactoring.observable_models_setup import *
from tudatpy.kernel.estimation_refactoring.observations_setup import *
