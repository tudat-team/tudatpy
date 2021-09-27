import tudatpy.kernel.interface.spice as spice_interface
# from tudatpy import load_standard_spice_kernels
# import pytest
# from pytest_bdd import scenario, given, when, then


# EXAMPLES: https://www.jetbrains.com/help/pycharm/pytest.html

def test_load_default_kernels():
    # load_standard_spice_kernels()
    spice.load_standard_kernels()


def test_clear_spice_kernels():
    spice.clear_kernels()

