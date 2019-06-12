from tudatpy import constants
import unittest


class TestPhysicalConstants(unittest.TestCase):
    """
    Testing of interfaced Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h
    """

    def test_sea_level_gravitational_constant(self):
        assert constants.SEA_LEVEL_GRAVITATIONAL_ACCELERATION == 9.80665

    def test_julian_day(self):
        assert constants.JULIAN_DAY == 86400.0

    def test_julian_day_long(self):
        assert constants.JULIAN_DAY_LONG == 86400.0

    def test_julian_year_in_days(self):
        assert constants.JULIAN_YEAR_IN_DAYS == 365.25

    def test_julian_year_in_days_long(self):
        assert constants.JULIAN_YEAR_IN_DAYS == 365.25

    def test_julian_year(self):
        assert constants.JULIAN_YEAR == 31557600.0

    def test_sidereal_day(self):
        assert constants.SIDEREAL_DAY == 86164.09054

    def test_sidereal_year_in_days(self):
        assert constants.SIDEREAL_YEAR_IN_DAYS == 365.25636

    def test_sidereal_year(self):
        assert constants.SIDEREAL_YEAR == 31558149.504

    def test_speed_of_light(self):
        assert constants.SPEED_OF_LIGHT == 299792458.0

    def test_speed_of_light_long(self):
        assert constants.SPEED_OF_LIGHT_LONG == 299792458.0

    def test_gravitational_constant(self):
        assert constants.GRAVITATIONAL_CONSTANT == 6.67259e-11

    def test_astronomical_unit(self):
        assert constants.ASTRONOMICAL_UNIT == 149597870691.0

    def test_molar_gas_constant(self):
        assert constants.MOLAR_GAS_CONSTANT == 8.3144598

    def test_planck_constant(self):
        assert constants.PLANCK_CONSTANT == 6.62606957e-34

    def test_boltzmann_constant(self):
        assert constants.BOLTZMANN_CONSTANT == 1.3806488e-23

    def test_stefan_boltzmann_constant(self):
        assert constants.STEFAN_BOLTZMANN_CONSTANT == 5.6703726225913323e-08

    def test_inverse_square_speed_of_light(self):
        assert constants.INVERSE_SQUARE_SPEED_OF_LIGHT == 1.1126500560536185e-17

    def test_inverse_quartic_speed_of_light(self):
        assert constants.INVERSE_QUARTIC_SPEED_OF_LIGHT == 1.2379901472361205e-34

    def test_inverse_quintic_speed_of_light(self):
        assert constants.INVERSE_QUINTIC_SPEED_OF_LIGHT == 4.129490633270435e-43

    def test_vacuum_permeability(self):
        assert constants.VACUUM_PERMEABILITY == 1.2566370614359173e-06

    def test_vacuum_permittivity(self):
        assert constants.VACUUM_PERMITTIVITY == 8.85418781762039e-12

    def test_lg_time_rate_term(self):
        assert constants.LG_TIME_RATE_TERM == 6.969290134e-10

    def test_lg_time_rate_term_long(self):
        assert constants.LG_TIME_RATE_TERM_LONG == 6.969290134e-10


class TestCelestialBodyConstants(unittest.TestCase):
    """
    Testing of interfaced Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h
    """

    def test_earth_equatorial_radius(self):
        assert constants.EARTH_EQUATORIAL_RADIUS == 6378136.6

    def test_earth_flattening_factor(self):
        assert constants.EARTH_FLATTENING_FACTOR == 298.25642

    def test_earth_geodesy_normalized_j2(self):
        assert constants.EARTH_GEODESY_NORMALIZED_J2 == -0.484165143790815E-03

    def test_sun_gravitational_parameter(self):
        assert constants.SUN_GRAVITATIONAL_PARAMETER == 1.32712440018e20

    def test_mercury_gravitational_parameter(self):
        assert constants.MERCURY_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 6023600.0

    def test_venus_gravitational_parameter(self):
        assert constants.VENUS_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 408523.71

    def test_earth_gravitational_parameter(self):
        assert constants.EARTH_GRAVITATIONAL_PARAMETER == 3.986004418E14

    def test_moon_gravitational_parameter(self):
        assert constants.MOON_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / (
                328900.56 * (1.0 + 81.30059))

    def test_mars_gravitational_parameter(self):
        assert constants.MARS_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 3098708.0

    def test_jupiter_gravitational_parameter(self):
        assert constants.JUPITER_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 1047.3486

    def test_saturn_gravitational_parameter(self):
        assert constants.SATURN_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 3497.898

    def test_uranus_gravitational_parameter(self):
        assert constants.URANUS_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 22902.98

    def test_neptune_gravitational_parameter(self):
        assert constants.NEPTUNE_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 19412.24

    def test_pluto_gravitational_parameter(self):
        assert constants.PLUTO_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 1.35e8

    def test_planet_names(self):
        assert constants.planet_names == {0: "Sun",
                                          1: "Mercury",
                                          2: "Venus",
                                          3: "Earth",
                                          4: "Mars",
                                          5: "Jupiter",
                                          6: "Saturn",
                                          7: "Uranus",
                                          8: "Neptune",
                                          9: "Pluto"}

    def test_planet_id_numbers(self):
        assert constants.planet_id_numbers == {"Sun": 0,
                                               "Mercury": 1,
                                               "Venus": 2,
                                               "Earth": 3,
                                               "Mars": 4,
                                               "Jupiter": 5,
                                               "Saturn": 6,
                                               "Uranus": 7,
                                               "Neptune": 8,
                                               "Pluto": 9}
