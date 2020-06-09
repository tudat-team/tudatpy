import tudatpy.spice_interface as spice_interface
import unittest


class TestSpiceInterface(unittest.TestCase):
    """
    Testing of interfaced Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h
    """

    def test_spice_interface_load_spice_kernels(self):
        spice_interface.load_standard_spice_kernels()
