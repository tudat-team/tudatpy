import sys

print(
    "Warning: tudatpy.interface.spice is deprecated and will be removed in a future "
    "release. Use tudatpy.data_input.environment_data.spice instead.",
    file=sys.stderr,
)

from tudatpy.data_input.environment_data.spice import *
