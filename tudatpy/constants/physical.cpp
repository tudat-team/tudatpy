// TODO: Disable warnings according to detected compiler. GCC/clang.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"

#include <boost/python.hpp>
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace bp = boost::python;
namespace tpc = tudat::physical_constants;

BOOST_PYTHON_MODULE (physical) {

        // Definition of physical constants.
        bp::scope().attr("SEA_LEVEL_GRAVITATIONAL_ACCELERATION") = tpc::SEA_LEVEL_GRAVITATIONAL_ACCELERATION;
        bp::scope().attr("JULIAN_DAY") = tpc::JULIAN_DAY;
        bp::scope().attr("JULIAN_DAY_LONG") = tpc::JULIAN_DAY_LONG;
        bp::scope().attr("JULIAN_YEAR_IN_DAYS") = tpc::JULIAN_YEAR_IN_DAYS;
        bp::scope().attr("JULIAN_YEAR_IN_DAYS_LONG") = tpc::JULIAN_YEAR_IN_DAYS_LONG;
        bp::scope().attr("JULIAN_YEAR") = tpc::JULIAN_YEAR;
        bp::scope().attr("SIDEREAL_DAY") = tpc::SIDEREAL_DAY;
        bp::scope().attr("SIDEREAL_YEAR_IN_DAYS") = tpc::SIDEREAL_YEAR_IN_DAYS;
        bp::scope().attr("SIDEREAL_YEAR") = tpc::SIDEREAL_YEAR;
        bp::scope().attr("SPEED_OF_LIGHT") = tpc::SPEED_OF_LIGHT;
        bp::scope().attr("SPEED_OF_LIGHT_LONG") = tpc::SPEED_OF_LIGHT_LONG;
        bp::scope().attr("GRAVITATIONAL_CONSTANT") = tpc::GRAVITATIONAL_CONSTANT;
        bp::scope().attr("ASTRONOMICAL_UNIT") = tpc::ASTRONOMICAL_UNIT;
        bp::scope().attr("SPECIFIC_GAS_CONSTANT_AIR") = tpc::SPECIFIC_GAS_CONSTANT_AIR;
        bp::scope().attr("MOLAR_GAS_CONSTANT") = tpc::MOLAR_GAS_CONSTANT;
        bp::scope().attr("PLANCK_CONSTANT") = tpc::PLANCK_CONSTANT;
        bp::scope().attr("BOLTZMANN_CONSTANT") = tpc::BOLTZMANN_CONSTANT;
        bp::scope().attr("STEFAN_BOLTZMANN_CONSTANT") = tpc::STEFAN_BOLTZMANN_CONSTANT;
        bp::scope().attr("INVERSE_SQUARE_SPEED_OF_LIGHT") = tpc::INVERSE_SQUARE_SPEED_OF_LIGHT;
        bp::scope().attr("INVERSE_CUBIC_SPEED_OF_LIGHT") = tpc::INVERSE_CUBIC_SPEED_OF_LIGHT;
        bp::scope().attr("INVERSE_QUARTIC_SPEED_OF_LIGHT") = tpc::INVERSE_QUARTIC_SPEED_OF_LIGHT;
        bp::scope().attr("INVERSE_QUINTIC_SPEED_OF_LIGHT") = tpc::INVERSE_QUINTIC_SPEED_OF_LIGHT;
        bp::scope().attr("VACUUM_PERMEABILITY") = tpc::VACUUM_PERMEABILITY;
        bp::scope().attr("VACUUM_PERMITTIVITY") = tpc::VACUUM_PERMITTIVITY;
        bp::scope().attr("LG_TIME_RATE_TERM") = tpc::LG_TIME_RATE_TERM;
        bp::scope().attr("LG_TIME_RATE_TERM_LONG") = tpc::LG_TIME_RATE_TERM_LONG;
}
#pragma GCC diagnostic pop