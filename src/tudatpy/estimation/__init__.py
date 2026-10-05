from tudatpy.kernel.estimation import *

# The wildcard import above binds the kernel ``observations`` submodule as an
# attribute of this package, which shadows TudatPy's Python subpackage of the
# same name. Dropping that attribute first is required: ``from . import X`` binds
# by attribute lookup on the parent and skips the submodule import when the
# attribute already exists, so without the ``del`` the kernel module would win
# and the helper functions would not be reachable through
# ``from tudatpy.estimation import observations``.
del observations
from . import observations
