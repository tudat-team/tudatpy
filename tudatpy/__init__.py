from .core import _constants
from .core import _interpolators
from .core import _spice_interface
from .core import _ephemerides
from .core import _reference_frames
from .core import _aerodynamics
from .core import _basic_astrodynamics
from .core import _gravitation
from .core import _numerical_integrators
from .core import _propagators
from .core import _orbital_element_conversions
from .core import _simulation_setup
from .core import _io
from .core import _unit_tests
from ._layer_simulation_setup import modify_simulation_setup
import sys

sys.modules["tudatpy.constants"] = _constants
sys.modules["constants"] = _constants

sys.modules["tudatpy.interpolators"] = _interpolators
sys.modules["interpolators"] = _interpolators

sys.modules["tudatpy.spice_interface"] = _spice_interface
sys.modules["spice_interface"] = _spice_interface

sys.modules["tudatpy.ephemerides"] = _ephemerides
sys.modules["ephemerides"] = _ephemerides

sys.modules["tudatpy.reference_frames"] = _reference_frames
sys.modules["reference_frames"] = _reference_frames

sys.modules["tudatpy.aerodynamics"] = _aerodynamics
sys.modules["aerodynamics"] = _aerodynamics

sys.modules["tudatpy.basic_astrodynamics"] = _basic_astrodynamics
sys.modules["basic_astrodynamics"] = _basic_astrodynamics

sys.modules["tudatpy.gravitation"] = _gravitation
sys.modules["gravitation"] = _gravitation

sys.modules["tudatpy.numerical_integrators"] = _numerical_integrators
sys.modules["numerical_integrators"] = _numerical_integrators

sys.modules["tudatpy.propagators"] = _propagators
sys.modules["propagators"] = _propagators

sys.modules["tudatpy.orbital_element_conversions"] = _orbital_element_conversions
sys.modules["orbital_element_conversions"] = _orbital_element_conversions

modify_simulation_setup(_simulation_setup)
sys.modules["tudatpy.simulation_setup"] = _simulation_setup
sys.modules["simulation_setup"] = _simulation_setup

sys.modules["tudatpy.io"] = _io
sys.modules["io"] = _io

sys.modules["tudatpy.unit_tests"] = _unit_tests
sys.modules["unit_tests"] = _unit_tests

__all__ = [
    'elements',
    'io',
    'constants',
    'interpolators',
    'spice_interface',
    'ephemerides',
    'reference_frames',
    'aerodynamics',
    'basic_astrodynamics',
    'gravitation',
    'numerical_integrators',
    'propagators',
    'orbital_element_conversions'
    'simulation_setup',
    'unit_tests',
    'prototype'
]

# Clean up namespace
del modify_simulation_setup
# del _constants
# del _interpolators
# del _spice_interface
# del _basic_astrodynamics
# del _gravitation
# del _numerical_integrators
# del _propagators
# del _orbital_element_conversions
# del _simulation_setup
