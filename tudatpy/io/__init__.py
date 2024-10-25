from tudatpy.data._support import *

# Deprecation warning
from warnings import warn
warn(
    "Importing from the `tudatpy.io` module is not support since tudatpy 0.8 " 
    "and now raises an error two minor versions hence. To eliminate this error "
    "import from `tudatpy.data` instead.",
DeprecationWarning, 2)
