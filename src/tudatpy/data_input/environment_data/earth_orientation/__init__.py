"""Raw Earth orientation data, independent of environment or estimation models.

Records use UTC Modified Julian Dates, radians, radians per second, and seconds.
Missing optional values are ``None``. Reading preserves raw coverage and does
not construct interpolators, rotation models, or station deformation models.
"""

from tudatpy.kernel.data_input.environment_data.earth_orientation import *
