from tudatpy import kernel
import pytest
import os


def test_no_null_bytes():
    """ For osx
    - https://github.com/tudat-team/tudatpy-feedstock/issues/1
    - https://gitter.im/conda-forge/conda-forge.github.io/archives/2019/01/22
    """
    assert '\x00' not in kernel.paths.get_resource_path()
    assert '\x00' not in kernel.paths.get_ephemeris_path()
    assert '\x00' not in kernel.paths.get_earth_orientation_path()
    assert '\x00' not in kernel.paths.get_quadrature_path()
    assert '\x00' not in kernel.paths.get_spice_kernel_path()
    assert '\x00' not in kernel.paths.get_atmosphere_tables_path()
    assert '\x00' not in kernel.paths.get_gravity_models_path()
    assert '\x00' not in kernel.paths.get_space_weather_path()


def test_paths_exist():
    """ For linux which broke from osx's null byte fix.
    """
    assert os.path.exists(kernel.paths.get_resource_path())
    assert os.path.exists(kernel.paths.get_ephemeris_path())
    assert os.path.exists(kernel.paths.get_earth_orientation_path())
    assert os.path.exists(kernel.paths.get_quadrature_path())
    assert os.path.exists(kernel.paths.get_spice_kernel_path())
    assert os.path.exists(kernel.paths.get_atmosphere_tables_path())
    assert os.path.exists(kernel.paths.get_gravity_models_path())
    assert os.path.exists(kernel.paths.get_space_weather_path())


def test_paths_for_conda_prefix():
    """ For testing conda-builds to ensure that prefix replacement occurs.
    """
    try:
        if os.environ["CONDA_PREFIX"]:
            assert os.path.join(os.environ["CONDA_PREFIX"], "resource") in kernel.paths.get_resource_path()
            assert os.path.join(os.environ["CONDA_PREFIX"], "resource") in kernel.paths.get_ephemeris_path()
            assert os.path.join(os.environ["CONDA_PREFIX"], "resource") in kernel.paths.get_earth_orientation_path()
            assert os.path.join(os.environ["CONDA_PREFIX"], "resource") in kernel.paths.get_quadrature_path()
            assert os.path.join(os.environ["CONDA_PREFIX"], "resource") in kernel.paths.get_spice_kernel_path()
            assert os.path.join(os.environ["CONDA_PREFIX"], "resource") in kernel.paths.get_atmosphere_tables_path()
            assert os.path.join(os.environ["CONDA_PREFIX"], "resource") in kernel.paths.get_gravity_models_path()
            assert os.path.join(os.environ["CONDA_PREFIX"], "resource") in kernel.paths.get_space_weather_path()
    except KeyError:
        pytest.skip("Reason: CONDA_PREFIX not found in env.")
