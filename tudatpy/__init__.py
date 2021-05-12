#from ._layer_propagation_setup import modify_propagation_setup
#from .kernel.simulation import propagation_setup
from ._version import *
from .kernel import constants
from .kernel import astro
from .kernel import interface
from .kernel import math
from .kernel import simulation

#modify_propagation_setup(propagation_setup)


__all__ = [
    '__version__',
    'apps',
    'bodies',
    'io',
    'prototype',
    'constants',
    'astro',
    'interface',
    'kernel',
    'math',
    'simulation'
]

# Clean up namespace
# del modify_simulation_setup
