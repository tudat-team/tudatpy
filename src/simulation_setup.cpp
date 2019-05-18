//
// Created by ggarrett on 18-5-19.
//

#include <boost/python.hpp>
#include "../../tudat/Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"

using namespace boost::python;
using namespace tudat::simulation_setup;

BOOST_PYTHON_MODULE(simulation_setup)
        {

            class_<BodySettings>("BodySettings")
                    .add_property("constant_mass", &BodySettings::constantMass)
                    .add_property("atmosphere_settings", &BodySettings::atmosphereSettings)
                    .add_property("ephemeris_settings", &BodySettings::ephemerisSettings)
                    .add_property("gravity_field_settings", &BodySettings::gravityFieldSettings)
                    .add_property("rotation_model_settings", &BodySettings::rotationModelSettings)
                    .add_property("shape_model_settings", &BodySettings::shapeModelSettings)
                    .add_property("radiation_pressure_settings", &BodySettings::radiationPressureSettings)
                    .add_property("aerodynamic_coefficient_settings", &BodySettings::aerodynamicCoefficientSettings)
                    .add_property("gravity_field_variation_settings", &BodySettings::gravityFieldVariationSettings)
                    .add_property("ground_station_settings", &BodySettings::groundStationSettings)
                    ;

        }
