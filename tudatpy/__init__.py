#from ._layer_propagation_setup import modify_propagation_setup
#from .kernel.simulation import propagation_setup
from ._version import *

#modify_propagation_setup(propagation_setup)


__all__ = [
    '__version__',
    'apps',
    'bodies',
    'io',
    'elements',
    'prototype',
    'kernel',
    'kernel.io',
    'kernel.constants',
    'kernel.astro',
    'kernel.interface',
    'kernel.math',
    'kernel.simulation',
    'kernel.aerodynamics',
    'kernel.unit_tests',
]

# Clean up namespace
# del modify_simulation_setup
