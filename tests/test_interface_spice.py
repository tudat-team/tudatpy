from tudatpy.interface import spice

# from tudatpy import load_standard_spice_kernels
# import pytest
# from pytest_bdd import scenario, given, when, then


# EXAMPLES: https://www.jetbrains.com/help/pycharm/pytest.html


def test_load_default_kernels() -> None:

    spice.clear_kernels()

    try:
        spice.load_standard_kernels()
        assert spice.get_total_count_of_kernels_loaded() == 13
    finally:
        spice.clear_kernels()
