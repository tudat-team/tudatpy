# IO IS DEPRECATED

# from tudatpy import io
# import pytest
# import os
# import sys


# def test_no_null_bytes():
#     """For osx
#     - https://github.com/tudat-team/tudatpy-feedstock/issues/1
#     - https://gitter.im/conda-forge/conda-forge.github.io/archives/2019/01/22
#     """
#     assert "\x00" not in kernel.io.get_resource_path()
#     assert "\x00" not in kernel.io.get_ephemeris_path()
#     assert "\x00" not in kernel.io.get_earth_orientation_path()
#     assert "\x00" not in kernel.io.get_quadrature_path()
#     assert "\x00" not in kernel.io.get_spice_kernel_path()
#     assert "\x00" not in kernel.io.get_atmosphere_tables_path()
#     assert "\x00" not in kernel.io.get_gravity_models_path()
#     assert "\x00" not in kernel.io.get_space_weather_path()


# def test_paths_exist():
#     """For linux which broke from osx's null byte fix."""
#     assert os.path.exists(kernel.io.get_resource_path())
#     assert os.path.exists(kernel.io.get_ephemeris_path())
#     assert os.path.exists(kernel.io.get_earth_orientation_path())
#     assert os.path.exists(kernel.io.get_quadrature_path())
#     assert os.path.exists(kernel.io.get_spice_kernel_path())
#     assert os.path.exists(kernel.io.get_atmosphere_tables_path())
#     assert os.path.exists(kernel.io.get_gravity_models_path())
#     assert os.path.exists(kernel.io.get_space_weather_path())


# def test_resource_paths_exist():
#     """For testing tudat resources exist in hidden directory."""
#     try:
#         if sys.platform == "win32" or sys.platform == "win64":
#             home_path = os.path.join(os.environ["HOMEDRIVE"], os.environ["HOMEPATH"])
#             print(f"Checking {home_path} for .tudat/resource existence.")
#         else:
#             home_path = os.environ["HOME"]
#             print(f"Checking {home_path} for .tudat/resource existence.")
#         assert os.path.exists(os.path.join(home_path, ".tudat/resource")) == True
#     except KeyError:
#         pytest.skip("Reason: CONDA_BUILD not found in env.")
