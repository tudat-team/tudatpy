// TODO: Disable warnings according to detected compiler. GCC/clang.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"

#include <boost/python.hpp>
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

// TODO: Figure out why the below imports result in a symbol error.
// Includes for BodySettings.
//#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"

//#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
//
//#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createBodyShapeModel.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createRadiationPressureInterface.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"
//#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/args.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/default_call_policies.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/import.hpp>
#include <boost/python/init.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost/python/object.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/self.hpp>
#include <boost/python/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include "boost/python/numpy.hpp"

namespace bn = boost::python::numpy;

using namespace boost::python;
using namespace tudat::physical_constants;
using namespace tudat::ephemerides;

typedef std::function<Eigen::Vector6d()> stateVectorFunction;

//boost::python::object get_state_vector_function_pywrapper(const bn::ndarray& state)
//{
//    auto func = get_string_function(name);
//    auto call_policies = boost::python::default_call_policies();
//    typedef Eigen::Vector6d func_sig;
//    return boost::python::make_function(func, call_policies, func_sig());
//}

//stateVectorFunction

BOOST_PYTHON_MODULE (constants) {
    // Definition of constants.
    boost::python::scope().attr("SEA_LEVEL_GRAVITATIONAL_ACCELERATION") = SEA_LEVEL_GRAVITATIONAL_ACCELERATION;
    boost::python::scope().attr("JULIAN_DAY") = JULIAN_DAY;
    boost::python::scope().attr("JULIAN_DAY_LONG") = JULIAN_DAY_LONG;
    boost::python::scope().attr("JULIAN_YEAR_IN_DAYS") = JULIAN_YEAR_IN_DAYS;
    boost::python::scope().attr("JULIAN_YEAR_IN_DAYS_LONG") = JULIAN_YEAR_IN_DAYS_LONG;
    boost::python::scope().attr("JULIAN_YEAR") = JULIAN_YEAR;
    boost::python::scope().attr("SIDEREAL_DAY") = SIDEREAL_DAY;
    boost::python::scope().attr("SIDEREAL_YEAR_IN_DAYS") = SIDEREAL_YEAR_IN_DAYS;
    boost::python::scope().attr("SIDEREAL_YEAR") = SIDEREAL_YEAR;
    boost::python::scope().attr("SPEED_OF_LIGHT") = SPEED_OF_LIGHT;
    boost::python::scope().attr("SPEED_OF_LIGHT_LONG") = SPEED_OF_LIGHT_LONG;
    boost::python::scope().attr("GRAVITATIONAL_CONSTANT") = GRAVITATIONAL_CONSTANT;
    boost::python::scope().attr("ASTRONOMICAL_UNIT") = ASTRONOMICAL_UNIT;
    boost::python::scope().attr("SPECIFIC_GAS_CONSTANT_AIR") = SPECIFIC_GAS_CONSTANT_AIR;
    boost::python::scope().attr("MOLAR_GAS_CONSTANT") = MOLAR_GAS_CONSTANT;
    boost::python::scope().attr("PLANCK_CONSTANT") = PLANCK_CONSTANT;
    boost::python::scope().attr("BOLTZMANN_CONSTANT") = BOLTZMANN_CONSTANT;
    boost::python::scope().attr("STEFAN_BOLTZMANN_CONSTANT") = STEFAN_BOLTZMANN_CONSTANT;
    boost::python::scope().attr("INVERSE_SQUARE_SPEED_OF_LIGHT") = INVERSE_SQUARE_SPEED_OF_LIGHT;
    boost::python::scope().attr("INVERSE_CUBIC_SPEED_OF_LIGHT") = INVERSE_CUBIC_SPEED_OF_LIGHT;
    boost::python::scope().attr("INVERSE_QUARTIC_SPEED_OF_LIGHT") = INVERSE_QUARTIC_SPEED_OF_LIGHT;
    boost::python::scope().attr("INVERSE_QUINTIC_SPEED_OF_LIGHT") = INVERSE_QUINTIC_SPEED_OF_LIGHT;
    boost::python::scope().attr("VACUUM_PERMEABILITY") = VACUUM_PERMEABILITY;
    boost::python::scope().attr("VACUUM_PERMITTIVITY") = VACUUM_PERMITTIVITY;
    boost::python::scope().attr("LG_TIME_RATE_TERM") = LG_TIME_RATE_TERM;
    boost::python::scope().attr("LG_TIME_RATE_TERM_LONG") = LG_TIME_RATE_TERM_LONG;

//    class_<Ephemeris>("Ephemeris", init<std::string, std::string>())
//            .def("reference_frame_origin", &Ephemeris::getReferenceFrameOrigin)
//            .def("reference_frame_orientation", &Ephemeris::getReferenceFrameOrientation);


//    boost::python::def("state_vector_function", get_state_vector_function_pywrapper);


    class_<ConstantEphemeris >("ConstantEphemeris",
                                                 init<
                                                         std::function<Eigen::Vector6d()>,
                                                         std::string,
                                                         std::string
                                                 >())


//    .add_property("get_cartesian_state", &ConstantEphemeris::getCartesianState)
//    .add_property("update_constant_state", &ConstantEphemeris::updateConstantState)
    ;

//    class_<BodySettings>("BodySettings")
//        .add_property("constant_mass", &BodySettings::constantMass)
//        .add_property("atmosphere_settings", &BodySettings::atmosphereSettings)
//        .add_property("ephemeris_settings", &BodySettings::ephemerisSettings)
//        .add_property("gravity_field_settings", &BodySettings::gravityFieldSettings)
//        .add_property("rotation_model_settings", &BodySettings::rotationModelSettings)
//        .add_property("shape_model_settings", &BodySettings::shapeModelSettings)
//        .add_property("radiation_pressure_settings", &BodySettings::radiationPressureSettings)
//        .add_property("aerodynamic_coefficient_settings", &BodySettings::aerodynamicCoefficientSettings)
//        .add_property("gravity_field_variation_settings", &BodySettings::gravityFieldVariationSettings)
//        .add_property("ground_station_settings", &BodySettings::groundStationSettings)
//        ;
}

#pragma GCC diagnostic pop