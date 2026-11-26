# Copyright 2025 - Tudat Development Team
# This file is part of tudatpy.

__version__ = "1.1.0.dev5"


import os as _os

from ._openmp import is_openmp_available as is_openmp_available

try:
    from ._openmp_build import ENABLED as _openmp_enabled
except ModuleNotFoundError:
    # The unconfigured source tree has no build metadata.
    _openmp_enabled = False

if _openmp_enabled and "OMP_WAIT_POLICY" not in _os.environ and "GOMP_SPINCOUNT" not in _os.environ:
    _os.environ["GOMP_SPINCOUNT"] = "30000"
