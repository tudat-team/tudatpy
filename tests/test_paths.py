from tudatpy import kernel


def test_no_null_bytes():
    """ For osx
    - https://github.com/tudat-team/tudatpy-feedstock/issues/1
    - https://gitter.im/conda-forge/conda-forge.github.io/archives/2019/01/22
    """
    assert '\x00' not in kernel.get_tudat_path()
    assert '\x00' not in kernel.get_tudat_data_path()
    assert '\x00' not in kernel.get_ephemeris_data_files_path()
    assert '\x00' not in kernel.get_earth_orientation_data_files_path()
    assert '\x00' not in kernel.get_spice_kernel_path()
    assert '\x00' not in kernel.get_atmosphere_tables_path()
    assert '\x00' not in kernel.get_gravity_models_path()
    assert '\x00' not in kernel.get_space_weather_data_path()
