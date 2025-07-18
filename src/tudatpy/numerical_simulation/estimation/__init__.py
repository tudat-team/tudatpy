import warnings
warnings.warn(
    "tudatpy.numerical_simulation.estimation is deprecated. " \
    "Use tudatpy.estimation.observations, tudatpy.estimation.observations_setup and/or " \
    "tudatpy.estimation.estimation_analysis instead.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.estimation.observations import *
from tudatpy.kernel.estimation.observations_setup import *
from tudatpy.kernel.estimation.estimation_analysis import *

