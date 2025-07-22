import warnings
warnings.warn(
    "tudatpy.numerical_simulation is deprecated. Use tudatpy.dynamics and/or tudatpy.estimation.",
    FutureWarning,
    stacklevel=1
)

from tudatpy.kernel.dynamics import *
from tudatpy.kernel.estimation import *

from tudatpy.kernel.dynamics.simulator import * 
from tudatpy.kernel.estimation.estimation_analysis import * 
