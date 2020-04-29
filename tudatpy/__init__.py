from .core import _constants
from .core import _interpolators
from .core import _spice_interface
from .core import _basic_astrodynamics
from .core import _gravitation
from .core import _numerical_integrators
from .core import _propagators
from .core import _orbital_element_conversions
from .core import _simulation_setup
import sys

sys.modules["tudatpy.constants"] = _constants
sys.modules["tudatpy.interpolators"] = _interpolators
sys.modules["tudatpy.spice_interface"] = _spice_interface
sys.modules["tudatpy.basic_astrodynamics"] = _basic_astrodynamics
sys.modules["tudatpy.gravitation"] = _gravitation
sys.modules["tudatpy.numerical_integrators"] = _numerical_integrators
sys.modules["tudatpy.propagators"] = _propagators
sys.modules["tudatpy.orbital_element_conversions"] = _orbital_element_conversions
sys.modules["tudatpy.simulation_setup"] = _simulation_setup

__all__ = [
    'constants',
    'interpolators',
    'spice_interface',
    'basic_astrodynamics',
    'gravitation',
    'numerical_integrators',
    'propagators',
    'orbital_element_conversions'
    'simulation_setup'
]
