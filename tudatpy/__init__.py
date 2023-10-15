from ._version import *
from tudatpy.kernel import constants
from tudatpy.kernel import astro
from tudatpy.kernel import interface
from tudatpy.kernel import math
from tudatpy.kernel import utils
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel import trajectory_design
#from tudatpy.kernel import io

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
