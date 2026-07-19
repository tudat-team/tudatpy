from pathlib import Path

from tudatpy.data_input import resource_paths

PATH_FUNCTIONS = (
    resource_paths.get_atmosphere_tables_path,
    resource_paths.get_default_output_path,
    resource_paths.get_earth_deformation_path,
    resource_paths.get_earth_orientation_path,
    resource_paths.get_ephemeris_path,
    resource_paths.get_gravity_models_path,
    resource_paths.get_nequick2_path,
    resource_paths.get_quadrature_path,
    resource_paths.get_resource_path,
    resource_paths.get_space_weather_path,
    resource_paths.get_spice_kernel_path,
    resource_paths.get_station_location_path,
    resource_paths.get_test_data_path,
    resource_paths.get_tudat_data_path,
    resource_paths.get_tudat_path,
)


def test_data_input_resource_paths_do_not_contain_null_characters():
    for path_function in PATH_FUNCTIONS:
        assert "\x00" not in path_function()


def test_core_data_input_resource_paths_exist():
    for path_function in (
        resource_paths.get_resource_path,
        resource_paths.get_ephemeris_path,
        resource_paths.get_earth_orientation_path,
        resource_paths.get_quadrature_path,
        resource_paths.get_spice_kernel_path,
        resource_paths.get_atmosphere_tables_path,
        resource_paths.get_gravity_models_path,
        resource_paths.get_space_weather_path,
        resource_paths.get_station_location_path,
        resource_paths.get_test_data_path,
    ):
        assert Path(path_function()).exists()
