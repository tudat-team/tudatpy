//
// Created by ggarrett on 18-5-19.
//

#include <boost/python.hpp>
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodyShapeModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRadiationPressureInterface.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"

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
