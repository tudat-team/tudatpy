from tudatpy._deprecation import deprecation_warning

deprecation_warning(
    old_name="tudatpy.interface.spice", new_name="tudatpy.data_input.environment_data.spice"
)

from tudatpy.data_input.environment_data.spice import *
