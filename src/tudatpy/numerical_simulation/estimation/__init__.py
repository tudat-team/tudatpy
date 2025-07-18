from tudatpy.kernel.numerical_simulation.estimation import *
# the above is kept for now because of what remains in expose_estimation_propagated_covariance

import warnings
warnings.warn(
    "tudatpy.numerical_simulation.estimation is deprecated. " \
    "Use tudatpy.estimation_refactoring.observations, tudatpy.estimation_refactoring.observations_setup and/or " \
    "tudatpy.estimation_refactoring.estimation_analysis instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.estimation_refactoring.observations import *
from tudatpy.kernel.estimation_refactoring.observations_setup import *
from tudatpy.kernel.estimation_refactoring.estimation_analysis import *

