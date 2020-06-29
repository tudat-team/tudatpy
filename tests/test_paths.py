from tudatpy import kernel


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
