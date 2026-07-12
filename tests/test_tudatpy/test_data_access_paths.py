from pathlib import Path

from tudatpy.data_access import paths

PATH_FUNCTIONS = (
    paths.get_atmosphere_tables_path,
    paths.get_default_output_path,
    paths.get_earth_deformation_path,
    paths.get_earth_orientation_path,
    paths.get_ephemeris_path,
    paths.get_gravity_models_path,
    paths.get_nequick2_path,
    paths.get_quadrature_path,
    paths.get_resource_path,
    paths.get_space_weather_path,
    paths.get_spice_kernel_path,
    paths.get_station_location_path,
    paths.get_test_data_path,
    paths.get_tudat_data_path,
    paths.get_tudat_path,
)


def test_data_access_paths_do_not_contain_null_characters():
    for path_function in PATH_FUNCTIONS:
        assert "\x00" not in path_function()


def test_core_data_access_paths_exist():
    for path_function in (
        paths.get_resource_path,
        paths.get_ephemeris_path,
        paths.get_earth_orientation_path,
        paths.get_quadrature_path,
        paths.get_spice_kernel_path,
        paths.get_atmosphere_tables_path,
        paths.get_gravity_models_path,
        paths.get_space_weather_path,
        paths.get_station_location_path,
        paths.get_test_data_path,
    ):
        assert Path(path_function()).exists()
