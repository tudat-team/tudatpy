from ._layer_simulation_setup import modify_simulation_setup
from .kernel import simulation_setup

modify_simulation_setup(simulation_setup)

__all__ = [
    'apps',
    'bodies',
    'io',
    'elements',
    'prototype',
    'kernel',
    'kernel.io',
    'kernel.constants',
    'kernel.interpolators',
    'kernel.spice_interface',
    'kernel.ephemerides',
    'kernel.reference_frames',
    'kernel.aerodynamics',
    'kernel.basic_astrodynamics',
    'kernel.gravitation',
    'kernel.numerical_integrators',
    'kernel.propagators',
    'kernel.orbital_element_conversions'
    'kernel.simulation_setup',
    'kernel.unit_tests',
]

# Clean up namespace
del modify_simulation_setup
