# Copyright 2025 - Tudat Development Team
# This file is part of tudatpy.


def is_openmp_available() -> bool:
    """Return whether the installed TudatPy kernel was built with OpenMP support."""
    from . import kernel

    return kernel.is_openmp_available()
