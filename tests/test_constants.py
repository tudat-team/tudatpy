import tudatpy.kernel.constants as constants


def test_constant_values():
    assert constants.SEA_LEVEL_GRAVITATIONAL_ACCELERATION == 9.80665

    assert constants.JULIAN_DAY == 86400.0

    assert constants.JULIAN_DAY_LONG == 86400.0

    assert constants.JULIAN_YEAR_IN_DAYS == 365.25

    assert constants.JULIAN_YEAR_IN_DAYS == 365.25

    assert constants.JULIAN_YEAR == 31557600.0

    assert constants.SIDEREAL_DAY == 86164.09054

    assert constants.SIDEREAL_YEAR_IN_DAYS == 365.25636

    assert constants.SIDEREAL_YEAR == 31558149.504

    assert constants.SPEED_OF_LIGHT == 299792458.0

    assert constants.SPEED_OF_LIGHT_LONG == 299792458.0

    assert constants.GRAVITATIONAL_CONSTANT == 6.67259e-11

    assert constants.ASTRONOMICAL_UNIT == 149597870691.0

    assert constants.MOLAR_GAS_CONSTANT == 8.3144598

    assert constants.PLANCK_CONSTANT == 6.62606957e-34

    assert constants.BOLTZMANN_CONSTANT == 1.3806488e-23

    assert constants.STEFAN_BOLTZMANN_CONSTANT == 5.6703726225913323e-08

    assert constants.INVERSE_SQUARE_SPEED_OF_LIGHT == 1.1126500560536185e-17

    assert constants.INVERSE_QUARTIC_SPEED_OF_LIGHT == 1.2379901472361205e-34

    assert constants.INVERSE_QUINTIC_SPEED_OF_LIGHT == 4.129490633270435e-43

    assert constants.VACUUM_PERMEABILITY == 1.2566370614359173e-06

    assert constants.VACUUM_PERMITTIVITY == 8.85418781762039e-12

    assert constants.LG_TIME_RATE_TERM == 6.969290134e-10

    assert constants.LG_TIME_RATE_TERM_LONG == 6.969290134e-10

    assert constants.EARTH_EQUATORIAL_RADIUS == 6378136.6

    assert constants.EARTH_FLATTENING_FACTOR == 298.25642

    assert constants.EARTH_GEODESY_NORMALIZED_J2 == -0.484165143790815E-03

    assert constants.SUN_GRAVITATIONAL_PARAMETER == 1.32712440018e20

    assert constants.MERCURY_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 6023600.0

    assert constants.VENUS_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 408523.71

    assert constants.EARTH_GRAVITATIONAL_PARAMETER == 3.986004418E14

    assert constants.MOON_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / (
            328900.56 * (1.0 + 81.30059))

    assert constants.MARS_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 3098708.0

    assert constants.JUPITER_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 1047.3486

    assert constants.SATURN_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 3497.898

    assert constants.URANUS_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 22902.98

    assert constants.NEPTUNE_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 19412.24

    assert constants.PLUTO_GRAVITATIONAL_PARAMETER == constants.SUN_GRAVITATIONAL_PARAMETER / 1.35e8
