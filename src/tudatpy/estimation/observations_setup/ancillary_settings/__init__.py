from tudatpy.kernel.estimation.observations_setup.ancillary_settings import *
from tudatpy._deprecation import keep_typo_backwards_compatible

import sys
keep_typo_backwards_compatible(sys.modules[__name__], "ancilliary", "ancillary")