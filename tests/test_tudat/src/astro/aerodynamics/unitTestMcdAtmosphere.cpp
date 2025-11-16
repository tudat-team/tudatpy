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

// Helper function to convert date to seconds since J2000
double convertDateToJ2000( int day, int month, int year, int hour, int min, int sec )
{
    double julianDay = basic_astrodynamics::convertCalendarDateToJulianDay( year, month, day, hour, min, static_cast< double >( sec ) );

    // J2000 epoch is JD 2451545.0
    double daysSinceJ2000 = julianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000;
    return daysSinceJ2000 * physical_constants::JULIAN_DAY;
}

// Test Case 1: High-res mode, 20km altitude, perturbation=1
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase1 )
{
    // Input parameters from INPUT_K1.txt
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 20000.0;  // meters above areoid
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    // Create MCD atmosphere: dust scenario=1, perturbation=1, high-res=1
    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 1, 1, 5.0, 0.0, 1 );

    // Get atmospheric properties
    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );
    double zonalWind = atmosphereModel->getZonalWind( );
    double meridionalWind = atmosphereModel->getMeridionalWind( );

    // Expected values from REF_OUTPUT_K1
    double expectedPressure = 7.57E+01;         // Pa
    double expectedDensity = 2.23E-03;          // kg/m^3
    double expectedTemperature = 1.77E+02;      // K
    double expectedZonalWind = -3.42E+01;       // m/s
    double expectedMeridionalWind = -7.82E+00;  // m/s

    // Tolerance for comparison (adjust based on expected accuracy)
    double relativeTolerance = 0.05;  // 5% relative error

    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( zonalWind, expectedZonalWind, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( meridionalWind, expectedMeridionalWind, relativeTolerance * 100.0 );

    std::cout << "Test Case 1 Results:" << std::endl;
    std::cout << "  Pressure: " << pressure << " Pa (expected: " << expectedPressure << ")" << std::endl;
    std::cout << "  Density: " << density << " kg/m^3 (expected: " << expectedDensity << ")" << std::endl;
    std::cout << "  Temperature: " << temperature << " K (expected: " << expectedTemperature << ")" << std::endl;
}

// Test Case 2: High-res mode, 20km altitude, perturbation=2
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase2 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 20000.0;
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 1, 2, 5.0, 0.0, 1 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Expected from REF_OUTPUT_K2
    double expectedPressure = 7.78E+01;
    double expectedDensity = 2.29E-03;
    double expectedTemperature = 1.78E+02;

    double relativeTolerance = 0.05;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );

    std::cout << "Test Case 2 completed." << std::endl;
}

// Test Case 3: High-res mode, 50km altitude, perturbation=3 with gravity waves
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase3 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 50000.0;
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 1, 3, 5.0, 16000.0, 1 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Expected from REF_OUTPUT_K3
    double expectedPressure = 1.80E+00;
    double expectedDensity = 7.31E-05;
    double expectedTemperature = 1.29E+02;

    double relativeTolerance = 0.05;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );

    std::cout << "Test Case 3 completed." << std::endl;
}

// Test Case 4: High-res mode, 150km altitude, perturbation=1
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase4 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 150000.0;
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 1, 1, 5.0, 0.0, 1 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Expected from REF_OUTPUT_K4
    double expectedPressure = 4.79E-06;
    double expectedDensity = 1.21E-10;
    double expectedTemperature = 1.93E+02;

    double relativeTolerance = 0.05;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );

    std::cout << "Test Case 4 completed." << std::endl;
}

// Test Case 5: High-res mode, 150km altitude, perturbation=2
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase5 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 150000.0;
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 1, 2, 5.0, 0.0, 1 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Expected from REF_OUTPUT_K5
    double expectedPressure = 4.63E-06;
    double expectedDensity = 1.19E-10;
    double expectedTemperature = 1.90E+02;

    double relativeTolerance = 0.05;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );

    std::cout << "Test Case 5 completed." << std::endl;
}

// Test Case 6: Low-res mode, 150km altitude, perturbation=3 with gravity waves
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase6 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 150000.0;
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 1, 3, 5.0, 16000.0, 0 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Expected from REF_OUTPUT_K6
    double expectedPressure = 4.79E-06;
    double expectedDensity = 9.47E-11;
    double expectedTemperature = 2.26E+02;

    double relativeTolerance = 0.05;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );

    std::cout << "Test Case 6 completed." << std::endl;
}

// Test Case 7: High-res mode, 20km altitude, dust scenario=7 (warm), perturbation=1
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase7 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 20000.0;
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 7, 1, 5.0, 0.0, 1 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Expected from REF_OUTPUT_K7
    double expectedPressure = 7.83E+01;
    double expectedDensity = 2.25E-03;
    double expectedTemperature = 1.82E+02;

    double relativeTolerance = 0.05;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );

    std::cout << "Test Case 7 completed." << std::endl;
}

// Test Case 8: High-res mode, 20km altitude, dust scenario=7 (warm), perturbation=2
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase8 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 20000.0;
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 7, 2, 5.0, 0.0, 1 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Expected from REF_OUTPUT_K8
    double expectedPressure = 7.98E+01;
    double expectedDensity = 2.29E-03;
    double expectedTemperature = 1.82E+02;

    double relativeTolerance = 0.05;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );

    std::cout << "Test Case 8 completed." << std::endl;
}

// Test Case 9: Low-res mode, 150km altitude, dust scenario=1, perturbation=1
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase9 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 150000.0;
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 1, 1, 5.0, 0.0, 0 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Expected from REF_OUTPUT_K9
    double expectedPressure = 4.79E-06;
    double expectedDensity = 1.21E-10;
    double expectedTemperature = 1.93E+02;

    double relativeTolerance = 0.05;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );

    std::cout << "Test Case 9 completed." << std::endl;
}

// Test Case 10: Low-res mode, 150km altitude, dust scenario=7 (warm), perturbation=2
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase10 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 150000.0;
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel =
            std::make_shared< aerodynamics::McdAtmosphereModel >( "", 7, 2, 5.0, 0.0, 0 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );

    // Expected from REF_OUTPUT_K10
    double expectedPressure = 8.26E-06;
    double expectedDensity = 1.71E-10;
    double expectedTemperature = 2.38E+02;

    double relativeTolerance = 0.05;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( density, expectedDensity, relativeTolerance * 100.0 );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, relativeTolerance * 100.0 );

    std::cout << "Test Case 10 completed." << std::endl;
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