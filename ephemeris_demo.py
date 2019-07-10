from tudatpy.simulation.environment import get_default_body_settings
from tudatpy.simulation.environment import create_bodies
from tudatpy.simulation.environment import BodySettings
from tudatpy.simulation.environment import EphemerisType
from tudatpy.simulation.environment import EphemerisSettings
from tudatpy.spice import load_standard_spice_kernels

test = get_default_body_settings(["Earth", "Venus"], 0, 100, 0.01)
body_map = create_bodies(test)
for i in range(100):
    print(body_map["Earth"].ephemeris.get_cartesian_state(i))
print(type(test))
