import tudatpy.kernel.spice_interface as spice_interface
# import pytest
# from pytest_bdd import scenario, given, when, then


# EXAMPLES: https://www.jetbrains.com/help/pycharm/pytest.html

def test_load_default_kernels():
    spice_interface.load_standard_spice_kernels()


def test_clear_spice_kernels():
    spice_interface.clear_spice_kernels()

