/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#if TUDAT_BUILD_WITH_MCD

#include "tudat/astro/aerodynamics/mcdAtmosphereModel.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
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
using namespace basic_astrodynamics;
using namespace propagators;
using namespace basic_mathematics;

BOOST_AUTO_TEST_SUITE( test_mcd_atmosphere )

// Test basic MCD atmosphere model functionality
BOOST_AUTO_TEST_CASE( testMcdAtmosphere )
{
    // Define tolerance for equality
    double tolerance = 1.0E-6;

    // Create MCD atmosphere model with default settings
    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel = std::make_shared< aerodynamics::McdAtmosphereModel >( );

    // Test parameters
    double altitude = 100.0E3;  // 100 km
    double longitude = unit_conversions::convertDegreesToRadians( 10.0 );
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double time = 0.0;  // J2000 epoch

    // Get atmospheric properties
    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Check that values are reasonable (using placeholder values)
    BOOST_CHECK( density > 0.0 );
    BOOST_CHECK( density < 1.0 );  // Should be less than surface density
    BOOST_CHECK( pressure > 0.0 );
    BOOST_CHECK( pressure < 1000.0 );    // Should be less than typical surface pressure
    BOOST_CHECK( temperature > 100.0 );  // Above absolute zero
    BOOST_CHECK( temperature < 300.0 );  // Below typical max temperature

    std::cout << "MCD Atmosphere Test Results:" << std::endl;
    std::cout << "  Density: " << density << " kg/m^3" << std::endl;
    std::cout << "  Pressure: " << pressure << " Pa" << std::endl;
    std::cout << "  Temperature: " << temperature << " K" << std::endl;
}

// Test MCD atmosphere with different dust scenarios
BOOST_AUTO_TEST_CASE( testMcdAtmosphereDustScenarios )
{
    double altitude = 50.0E3;  // 50 km
    double longitude = 0.0;
    double latitude = 0.0;
    double time = 0.0;

    // Test different dust scenarios
    for( int dustScenario = 1; dustScenario <= 8; dustScenario++ )
    {
        std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
                std::make_shared< aerodynamics::McdAtmosphereModel >( "", dustScenario, 0, 0.0, 0.0, 0 );

        double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );

        BOOST_CHECK( density > 0.0 );

        std::cout << "Dust scenario " << dustScenario << " - Density: " << density << " kg/m^3" << std::endl;
    }
}
// Test MCD atmosphere in propagation
BOOST_AUTO_TEST_CASE( testMcdAtmosphereInPropagation )
{
    using namespace tudat;
    using namespace aerodynamics;
    using namespace simulation_setup;
    using namespace numerical_integrators;
    using namespace basic_astrodynamics;
    using namespace propagators;

    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Create Mars with MCD atmosphere
    BodyListSettings defaultBodySettings = getDefaultBodySettings( { "Mars" } );
    defaultBodySettings.at( "Mars" )->atmosphereSettings = mcdAtmosphereSettings( );
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Create vehicle
    double vehicleMass = 5.0E3;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

    // Set aerodynamic coefficients
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >( 2.0,
                                                                        4.0,
                                                                        Eigen::Vector3d::Zero( ),
                                                                        Eigen::Vector3d::UnitX( ),
                                                                        Eigen::Vector3d::Zero( ),
                                                                        negative_aerodynamic_frame_coefficients,
                                                                        negative_aerodynamic_frame_coefficients );
    bodies.at( "Vehicle" )
            ->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );

    // Define accelerations
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Mars" );

    // Set initial state
    Eigen::Vector6d systemInitialState;
    systemInitialState << 3500.0E3, 0.0, 0.0, 0.0, 3500.0, 0.0;

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Set dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
            std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Vehicle", "Mars" ) );
    dependentVariables.push_back(
            std::make_shared< SingleDependentVariableSaveSettings >( local_density_dependent_variable, "Vehicle", "Mars" ) );

    // Set propagation settings
    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );

    std::shared_ptr< IntegratorSettings<> > integratorSettings = std::make_shared< IntegratorSettings<> >( rungeKutta4, 0.0, 10.0 );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                accelerationModelMap,
                                                                                bodiesToPropagate,
                                                                                systemInitialState,
                                                                                0.0,
                                                                                integratorSettings,
                                                                                terminationSettings,
                                                                                cowell,
                                                                                dependentVariables );

    // Create simulation object and propagate
    SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, translationalPropagatorSettings );

    BOOST_CHECK_EQUAL( dynamicsSimulator.getSingleArcPropagationResults( )->getPropagationIsPerformed( ), true );
    BOOST_CHECK_EQUAL( dynamicsSimulator.getSingleArcPropagationResults( )->integrationCompletedSuccessfully( ), true );

    std::cout << "MCD atmosphere propagation test completed successfully." << std::endl;
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat

#endif  // TUDAT_BUILD_WITH_MCD