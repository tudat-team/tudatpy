from ._layer_simulation_setup import modify_simulation_setup
from .kernel import simulation_setup
from .kernel import get_spice_kernel_path
from .kernel import spice_interface

modify_simulation_setup(simulation_setup)


# Identical function from spice interface, added for static library fix.
def load_standard_spice_kernels(alternative_kernels: list = None) -> None:
    spice_interface.load_kernel(get_spice_kernel_path() + 'pck00010.tpc')
    spice_interface.load_kernel(get_spice_kernel_path() + 'gm_de431.tpc')
    if alternative_kernels:
        for alternative_kernel in alternative_kernels:
            spice_interface.load_kernel(alternative_kernel)
    else:
        spice_interface.load_kernel(get_spice_kernel_path() + 'tudat_merged_spk_kernel.bsp')
    spice_interface.load_kernel(get_spice_kernel_path() + 'naif0012.tls')


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
