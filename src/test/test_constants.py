import core
import pytest


def test_sea_level_gravitational_constant():
    assert core.SEA_LEVEL_GRAVITATIONAL_ACCELERATION == 9.80665


def test_julian_day():
    assert core.JULIAN_DAY == 86400.0


def test_julian_day_long():
    assert core.JULIAN_DAY_LONG == 86400.0


def test_julian_year_in_days():
    assert core.JULIAN_YEAR_IN_DAYS == 365.25


def test_julian_year_in_days_long():
    assert core.JULIAN_YEAR_IN_DAYS == 365.25


def test_julian_year():
    assert core.JULIAN_YEAR == 31557600.0


def test_sidereal_day():
    assert core.SIDEREAL_DAY == 86164.09054


def test_sidereal_year_in_days():
    assert core.SIDEREAL_YEAR_IN_DAYS == 365.25636


def test_sidereal_year():
    assert core.SIDEREAL_YEAR == 31558149.504


def test_speed_of_light():
    assert core.SPEED_OF_LIGHT == 299792458.0


def test_speed_of_light_long():
    assert core.SPEED_OF_LIGHT_LONG == 299792458.0


def test_gravitational_constant():
    assert core.GRAVITATIONAL_CONSTANT == 6.67259e-11


def test_astronomical_unit():
    assert core.ASTRONOMICAL_UNIT == 149597870691.0


def test_molar_gas_constant():
    assert core.MOLAR_GAS_CONSTANT == 8.3144598


def test_planck_constant():
    assert core.PLANCK_CONSTANT == 6.62606957e-34


def test_boltzmann_constant():
    assert core.BOLTZMANN_CONSTANT == 1.3806488e-23


def test_stefan_boltzmann_constant():
    assert core.STEFAN_BOLTZMANN_CONSTANT == 5.6703726225913323e-08


def test_inverse_square_speed_of_light():
    assert core.INVERSE_SQUARE_SPEED_OF_LIGHT == 1.1126500560536185e-17


def test_inverse_quartic_speed_of_light():
    assert core.INVERSE_QUARTIC_SPEED_OF_LIGHT == 1.2379901472361205e-34


def test_inverse_quintic_speed_of_light():
    assert core.INVERSE_QUINTIC_SPEED_OF_LIGHT == 4.129490633270435e-43


def test_vacuum_permeability():
    assert core.VACUUM_PERMEABILITY == 1.2566370614359173e-06


def test_vacuum_permittivity():
    assert core.VACUUM_PERMITTIVITY == 8.85418781762039e-12


def test_lg_time_rate_term():
    assert core.LG_TIME_RATE_TERM == 6.969290134e-10


def test_lg_time_rate_term_long():
    assert core.LG_TIME_RATE_TERM_LONG == 6.969290134e-10
