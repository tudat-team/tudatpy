from tudatpy.simulation.environment import get_default_body_settings
from tudatpy.simulation.environment import create_bodies
from tudatpy.simulation.environment import BodySettings
from tudatpy.simulation.environment import EphemerisType
from tudatpy.simulation.environment import EphemerisSettings
# from tudatpy.simulation.environment import load_standard_spice_kernels
from tudatpy.spice import load_standard_spice_kernels

import unittest
import math

class TestBody(unittest.TestCase):

    def test_body_settings_basic(self):
        test = BodySettings()
        assert(math.isnan(test.constant_mass) == True)

    def test_base_ephemeris_type(self):
        test = EphemerisType
        assert(str(test.approximate_planet_positions) == "approximate_planet_positions")
        assert(str(test.direct_spice_ephemeris) == "direct_spice_ephemeris")
        assert(str(test.interpolated_spice) == "interpolated_spice")
        assert(str(test.constant_ephemeris) == "constant_ephemeris")
        assert(str(test.kepler_ephemeris) == "kepler_ephemeris")
        assert(str(test.custom_ephemeris) == "custom_ephemeris")

    def test_ephemeris_settings(self):
        test = EphemerisSettings(
            EphemerisType.constant_ephemeris,
            "SSB",
            "ECLIPJ2000"
        )
        assert(test.ephemeris_type == EphemerisType.constant_ephemeris)
        assert(test.frame_origin == "SSB")
        assert(test.frame_orientation == "ECLIPJ2000")
        assert(test.multi_arc_ephemeris == False)
        test.reset_frame_origin("O")
        assert(test.frame_origin == "O")
        test.reset_make_multi_arc_ephemeris(True)
        assert(test.multi_arc_ephemeris == True)

    def test_default_body_settings(self):
        load_standard_spice_kernels()
        test = get_default_body_settings(["Earth", "Venus"])
        body_map = create_bodies(test)
        print(body_map)
        print(type(test))