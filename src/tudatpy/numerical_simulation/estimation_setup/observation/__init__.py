import warnings
warnings.warn(
    "tudatpy.numerical_simulation.estimation_setup.observation is deprecated. " \
    "Use tudatpy.estimation.observable_models_setup and/or tudatpy.estimation.observations_setup instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.estimation.observable_models_setup import *
from tudatpy.kernel.estimation.observations_setup import *
