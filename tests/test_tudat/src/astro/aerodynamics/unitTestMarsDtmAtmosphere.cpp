/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include "fstream"
#include "iostream"

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

using namespace tudat::aerodynamics;


namespace tudat
{
namespace unit_tests
{

    using namespace tudat;
using namespace aerodynamics;
using namespace simulation_setup;
using namespace numerical_integrators;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace propagators;
using namespace basic_mathematics;
using namespace basic_astrodynamics;

BOOST_AUTO_TEST_SUITE(test_mars_dtm_atmosphere)

BOOST_AUTO_TEST_CASE(testMarsDtmAtmosphere)
{
    // Define tolerance for equality
    double tolerance = 1.0E-15;
    std::shared_ptr<AtmosphereSettings> marsDtmAtmosphereSettings;
    marsDtmAtmosphereSettings = std::make_shared<MarsDtmAtmosphereSettings>();

    // Create a minimal SystemOfBodies for the atmosphere model
    SystemOfBodies bodies;
    bodies.createEmptyBody("Mars");

    std::shared_ptr<aerodynamics::AtmosphereModel> marsAtmosphereModel = createAtmosphereModel(
            marsDtmAtmosphereSettings, "Mars", bodies);
    std::shared_ptr<MarsDtmAtmosphereModel> atmosphereModel =
            std::dynamic_pointer_cast<MarsDtmAtmosphereModel>(
                    createAtmosphereModel(marsDtmAtmosphereSettings, "Mars", bodies));
    int alt_km = 400E3;
    int time = 86400;
    double latitude = unit_conversions::convertDegreesToRadians(15.0);
    double longitude = unit_conversions::convertDegreesToRadians(10.0);
    double rho = atmosphereModel->getDensity(alt_km, longitude, latitude, time);

    // reference density for this specific inputs
    double rho_ref = 2.57862799E-14;
    double rho_diff = std::abs(rho - rho_ref);
    // Check density
    BOOST_CHECK_SMALL(rho - rho_ref, tolerance);
}

// Test is Mars DTM model is being used properly (no crashes)
BOOST_AUTO_TEST_CASE(testMarsDtmAtmosphereInPropagation)
{
    using namespace tudat;
    using namespace aerodynamics;
    using namespace simulation_setup;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create Earth object
    BodyListSettings defaultBodySettings =
        getDefaultBodySettings( { "Mars" } );
    defaultBodySettings.at( "Mars" )->atmosphereSettings = marsDtmAtmosphereSettings( );
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Create vehicle object.
    double vehicleMass = 5.0E3;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

    // Set aerodynamic coefficients.
    std::shared_ptr<AerodynamicCoefficientSettings> aerodynamicCoefficientSettings =
        std::make_shared<ConstantAerodynamicCoefficientSettings>(
            2.0, 4.0, Eigen::Vector3d::Zero( ), Eigen::Vector3d::UnitX( ), Eigen::Vector3d::Zero( ),
            negative_aerodynamic_frame_coefficients, negative_aerodynamic_frame_coefficients );
    bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ));

    // Define acceleration model settings.
    SelectedAccelerationMap accelerationMap;
    std::vector<std::string> bodiesToPropagate;
    std::vector<std::string> centralBodies;
    std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared<AccelerationSettings>( aerodynamic ));
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Mars" );

    // Set initial state
    Eigen::Vector6d systemInitialState;
    systemInitialState << 3378.0E3, 3000E3, 4200E3, 0.0, 0.0, 0.0;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Set variables to save
    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> > dependentVariables;
    dependentVariables.push_back(
        std::make_shared<SingleDependentVariableSaveSettings>(
            altitude_dependent_variable, "Vehicle", "Mars" ));
    dependentVariables.push_back(
        std::make_shared<BodyAerodynamicAngleVariableSaveSettings>(
            "Vehicle", reference_frames::longitude_angle ));
    dependentVariables.push_back(
        std::make_shared<BodyAerodynamicAngleVariableSaveSettings>(
            "Vehicle", reference_frames::latitude_angle ));
    dependentVariables.push_back(
        std::make_shared<SingleDependentVariableSaveSettings>(
            local_density_dependent_variable, "Vehicle", "Mars" ));
    dependentVariables.push_back(
        std::make_shared<SingleAccelerationDependentVariableSaveSettings>(
            aerodynamic, "Vehicle", "Mars", 0 ));

    // Set propagation/integration settings
    std::shared_ptr<PropagationTimeTerminationSettings> terminationSettings =
        std::make_shared<propagators::PropagationTimeTerminationSettings>( 1000.0 );
    std::shared_ptr<IntegratorSettings<> > integratorSettings =
        std::make_shared<IntegratorSettings<> >
            ( rungeKutta4, 0.0, 10.0 );
    std::shared_ptr<tudat::propagators::TranslationalStatePropagatorSettings<double> >
        translationalPropagatorSettings =
        std::make_shared<TranslationalStatePropagatorSettings < double> >
        ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, 0.0,
            integratorSettings, terminationSettings,
            cowell, dependentVariables );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator<> dynamicsSimulator(
        bodies, translationalPropagatorSettings );

    BOOST_CHECK_EQUAL( dynamicsSimulator.getSingleArcPropagationResults( )->getPropagationIsPerformed( ), true );
    BOOST_CHECK_EQUAL( dynamicsSimulator.getSingleArcPropagationResults( )->integrationCompletedSuccessfully( ), true );

}



BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests

}//namespace tudat
