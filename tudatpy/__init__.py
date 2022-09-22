from ._version import *
from .kernel import constants
from .kernel import astro
from .kernel import interface
from .kernel import math
from .kernel import utils
from .kernel import numerical_simulation
from .kernel import trajectory_design

numerical_simulation.propagation_setup.integrator = math.integrator_setup

#from .kernel import io

__all__ = [
    '__version__',
    'apps',
    'bodies',
    'utils',
    #'io',
    'prototype',
    'constants',
    'astro',
    'interface',
    'kernel',
    'math',
    'numerical_simulation'
]

# Clean up namespace
# del modify_simulation_setup
