import tudatpy.interface.spice as spice_interface

# from tudatpy import load_standard_spice_kernels
# import pytest
# from pytest_bdd import scenario, given, when, then


# EXAMPLES: https://www.jetbrains.com/help/pycharm/pytest.html


def test_load_default_kernels():
    # load_standard_spice_kernels()
    spice_interface.load_standard_kernels()


def test_clear_spice_kernels():
    spice_interface.clear_kernels()
