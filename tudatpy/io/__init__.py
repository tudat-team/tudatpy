from tudatpy.data._support import *

# Deprecation warning
from warnings import warn
warn(
    "Importing from the `tudatpy.io` module is deprecated since tudatpy 0.7 " 
    "and will raise an error two minor versions hence. To eliminate this warning "
    "import from `tudatpy.data` instead.",
DeprecationWarning, 2)

__all__ = [
    "save2txt",
    "save_time_history_to_file",
    "get_resource_path",
    "get_ephemeris_path",
    "get_earth_orientation_path",
    "get_quadrature_path",
    "get_spice_kernel_path",
    "get_atmosphere_tables_path",
    "get_gravity_models_path",
    "get_space_weather_path",
    "read_vector_history_from_file",
    "read_matrix_history_from_file",
    "missile_DATCOM_data",
    "DynamicCoefficientNames",
    "StaticCoefficientNames"
]
