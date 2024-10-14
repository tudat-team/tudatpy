/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
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


#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/basics/testMacros.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"

namespace tudat
{
namespace unit_tests
{

using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace spice_interface;
using namespace input_output;
using namespace reference_frames;
using namespace propagators;

BOOST_AUTO_TEST_SUITE( test_environment_updater )

//! Test set up of point mass gravitational accelerations, both direct and third-body.
BOOST_AUTO_TEST_CASE( test_centralGravityEnvironmentUpdate )
{

    double initialTime = 86400.0;
    double finalTime = 2.0 * 86400.0;

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Get settings for celestial bodies
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 10.0 * 86400.0 ), "Earth" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Sun", 0.0,10.0 * 86400.0 ), "Sun" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Moon", 0.0,10.0 * 86400.0 ), "Moon" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Mars", 0.0,10.0 * 86400.0 ), "Mars" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Venus", 0.0,10.0 * 86400.0 ), "Venus" );

    // Create settings for sh gravity field.
    double gravitationalParameter = 398600.4418E9;
    Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );
    Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );
    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                gravitationalParameter, 6378.0E3, cosineCoefficients, sineCoefficients,
                "IAU_Earth" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    // Define variables used in tests.
    SelectedAccelerationMap accelerationSettingsMap;
    std::map< std::string, std::string > centralBodies;
    std::vector< std::string > propagatedBodyList;
    std::vector< std::string > centralBodyList;

    std::unordered_map< IntegratedStateType, Eigen::VectorXd > integratedStateToSet;
    double testTime = 0.0;

    {

        // Define (arbitrary) test state.
        Eigen::VectorXd testState = ( Eigen::VectorXd( 6 ) << 1.44E6, 2.234E8, -3343.246E7, 1.2E4, 1.344E3, -22.343E3 ).finished( );
        integratedStateToSet[ translational_state ] = testState;
        testTime = 2.0 * 86400.0;

        // Test if environment is updated for inly central gravity accelerations
        {
            // Define accelerations
            accelerationSettingsMap[ "Moon" ][ "Sun" ].push_back(
                        std::make_shared< AccelerationSettings >( point_mass_gravity ) );
            accelerationSettingsMap[ "Moon" ][ "Earth" ].push_back(
                        std::make_shared< AccelerationSettings >( point_mass_gravity ) );

            // Define origin of integration to be barycenter.
            centralBodies[ "Moon" ] = "SSB";
            propagatedBodyList.push_back( "Moon" );
            centralBodyList.push_back( centralBodies[ "Moon" ] );

            // Create accelerations
            AccelerationMap accelerationsMap = createAccelerationModelsMap(
                        bodies, accelerationSettingsMap, centralBodies );

            // Create environment update settings.
            std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                            "Moon", centralBodies[ "Moon" ], bodies, initialTime ), finalTime );
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                    createEnvironmentUpdaterSettings< double >( propagatorSettings, bodies );

            // Test update settings
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_translational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_translational_state_update ).size( ), 2 );

            // Create and call updater.
            std::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                    createEnvironmentUpdaterForDynamicalEquations< double, double >(
                        propagatorSettings, bodies );
            updater->updateEnvironment( testTime, integratedStateToSet );

            // Test if Earth, Sun and Moon are updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Earth" )->getState( ),
                        bodies.at( "Earth" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Sun" )->getState( ),
                        bodies.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Moon" )->getState( ), testState,
                        std::numeric_limits< double >::epsilon( ) );

            // Test if Mars is not updated
            bool exceptionIsCaught = false;
            try
            {
                bodies.at( "Mars" )->getState( );
            }
            catch( ... )
            {
                exceptionIsCaught = true;
            }
            BOOST_CHECK_EQUAL( exceptionIsCaught, true );

            // Update environment to new time, and state from environment.
            updater->updateEnvironment(
                        0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ), { translational_state } );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Earth" )->getState( ),
                        bodies.at( "Earth" )->getEphemeris( )->getCartesianState( 0.5 * testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Sun" )->getState( ),
                        bodies.at( "Sun" )->getEphemeris( )->getCartesianState( 0.5 * testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Moon" )->getState( ),
                        bodies.at( "Moon" )->getEphemeris( )->getCartesianState( 0.5 * testTime ),
                        std::numeric_limits< double >::epsilon( ) );
        }

        // Test third body acceleration updates.
        {
            accelerationSettingsMap.clear( );
            centralBodies.clear( );
            propagatedBodyList.clear( );
            centralBodyList.clear( );

            // Set acceleration models.
            accelerationSettingsMap[ "Moon" ][ "Sun" ].push_back(
                        std::make_shared< AccelerationSettings >( point_mass_gravity ) );
            accelerationSettingsMap[ "Moon" ][ "Mars" ].push_back(
                        std::make_shared< AccelerationSettings >( point_mass_gravity ) );

            // Define origin of integration
            centralBodies[ "Moon" ] = "Earth";
            propagatedBodyList.push_back( "Moon" );
            centralBodyList.push_back( centralBodies[ "Moon" ] );

            // Create accelerations
            AccelerationMap accelerationsMap = createAccelerationModelsMap(
                        bodies, accelerationSettingsMap, centralBodies );

            // Create environment update settings.
            std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                            "Moon", centralBodies[ "Moon" ], bodies, initialTime ), finalTime );
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                    createEnvironmentUpdaterSettings< double >( propagatorSettings, bodies );

            // Test update settings
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_translational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_translational_state_update ).size( ), 3 );

            // Create and call updater.
            std::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                    createEnvironmentUpdaterForDynamicalEquations< double, double >(
                        propagatorSettings, bodies );
            updater->updateEnvironment( testTime, integratedStateToSet );

            // Test if Earth, Sun, Mars and Moon are updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Earth" )->getState( ),
                        bodies.at( "Earth" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Sun" )->getState( ),
                        bodies.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Moon" )->getState( ), testState,
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Mars" )->getState( ),
                        bodies.at( "Mars" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );

            // Test if Mars is not updated
            bool exceptionIsCaught = false;
            try
            {
                bodies.at( "Venus" )->getState( );
            }
            catch( ... )
            {
                exceptionIsCaught = true;
            }
            BOOST_CHECK_EQUAL( exceptionIsCaught, true );


            // Update environment to new time, and state from environment.
            updater->updateEnvironment(
                        0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ), { translational_state } );

        }

        // Test spherical harmonic acceleration update
        {
            accelerationSettingsMap.clear( );
            centralBodies.clear( );
            propagatedBodyList.clear( );
            centralBodyList.clear( );

            // Set acceleration models.
            accelerationSettingsMap[ "Moon" ][ "Sun" ].push_back(
                        std::make_shared< AccelerationSettings >( point_mass_gravity ) );
            accelerationSettingsMap[ "Moon" ][ "Earth" ].push_back(
                        std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );
            accelerationSettingsMap[ "Moon" ][ "Mars" ].push_back(
                        std::make_shared< AccelerationSettings >( point_mass_gravity ) );

            // Define origin of integration
            centralBodies[ "Moon" ] = "Earth";
            propagatedBodyList.push_back( "Moon" );
            centralBodyList.push_back( centralBodies[ "Moon" ] );

            // Create accelerations
            AccelerationMap accelerationsMap = createAccelerationModelsMap(
                        bodies, accelerationSettingsMap, centralBodies );

            // Create environment update settings.
            std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                            "Moon", centralBodies[ "Moon" ], bodies, initialTime ), finalTime );
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                    createEnvironmentUpdaterSettings< double >( propagatorSettings, bodies );

            // Test update settings
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 3 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_translational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_translational_state_update ).size( ), 3 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_rotational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_rotational_state_update ).size( ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( spherical_harmonic_gravity_field_update ), 1 );

            // Create and call updater.
            std::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                    createEnvironmentUpdaterForDynamicalEquations< double, double >(
                        propagatorSettings, bodies );
            updater->updateEnvironment( testTime, integratedStateToSet );

            // Test if Earth, Sun, Mars and Moon are updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Earth" )->getState( ),
                        bodies.at( "Earth" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Sun" )->getState( ),
                        bodies.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Moon" )->getState( ), testState,
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Mars" )->getState( ),
                        bodies.at( "Mars" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );


            // Test if Mars is not updated
            bool exceptionIsCaught = false;
            try
            {
                bodies.at( "Venus" )->getState( );
            }
            catch( ... )
            {
                exceptionIsCaught = true;
            }
            BOOST_CHECK_EQUAL( exceptionIsCaught, true );

            // Test if Earth rotation is updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Earth" )->getCurrentRotationToGlobalFrame( ).toRotationMatrix( ),
                        bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ).toRotationMatrix( ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Earth" )->getCurrentRotationToLocalFrame( ).toRotationMatrix( ),
                        bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToTargetFrame( testTime ).toRotationMatrix( ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Earth" )->getCurrentRotationMatrixDerivativeToGlobalFrame( ),
                        bodies.at( "Earth" )->getRotationalEphemeris( )->getDerivativeOfRotationToBaseFrame( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodies.at( "Earth" )->getCurrentRotationMatrixDerivativeToLocalFrame( ),
                        bodies.at( "Earth" )->getRotationalEphemeris( )->getDerivativeOfRotationToTargetFrame( testTime ),
                        std::numeric_limits< double >::epsilon( ) );


            // Test if Mars is not updated
            exceptionIsCaught = false;
            try
            {
                bodies.at( "Mars" )->getCurrentRotationToLocalFrame( ).toRotationMatrix( );
            }
            catch( ... )
            {
                exceptionIsCaught = true;
            }
            BOOST_CHECK_EQUAL( exceptionIsCaught, true );

            exceptionIsCaught = false;
            try
            {
                bodies.at( "Mars" )->getCurrentRotationMatrixDerivativeToGlobalFrame( );
            }
            catch( ... )
            {
                exceptionIsCaught = true;
            }
            BOOST_CHECK_EQUAL( exceptionIsCaught, true );

            // Update environment to new time, and state from environment.
            updater->updateEnvironment(
                        0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ), { translational_state } );

        }
    }
}

//! Dummy function to calculate mass as a function of time.
double getBodyMass( const double time )
{
    return 1000.0 - 500.0 * ( time / ( 10.0 * 86400.0 ) );
}

//! Test radiation pressure acceleration
BOOST_AUTO_TEST_CASE( test_NonConservativeForceEnvironmentUpdate )
{
    double initialTime = 86400.0;
    double finalTime = 2.0 * 86400.0;

    using namespace tudat::simulation_setup;
    using namespace tudat;

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Get settings for celestial bodies
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 10.0 * 86400.0 ), "Earth" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Sun", 0.0,10.0 * 86400.0 ), "Sun" );

    // Get settings for vehicle
    double area = 2.34;
    double coefficient = 1.2;
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->radiationPressureSettings[ "Sun" ] =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >( "Sun", area, coefficient );
    bodySettings.at( "Vehicle" )->ephemerisSettings =
            std::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( ) << 7000.0E3, 0.05, 0.3, 0.0, 0.0, 0.0 ).finished( ),
                0.0, spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
                getApolloCoefficientInterface( ) );
    bodies.at( "Vehicle" )->setBodyMassFunction( &getBodyMass );
    


    // Define test time and state.
    double testTime = 2.0 * 86400.0;
    std::unordered_map< IntegratedStateType, Eigen::VectorXd > integratedStateToSet;
    Eigen::VectorXd testState = 1.1 * bodies.at( "Vehicle" )->getEphemeris( )->getCartesianState( testTime ) +
            bodies.at( "Earth" )->getEphemeris( )->getCartesianState( testTime );
    integratedStateToSet[ translational_state ] = testState;

    {
        // Define settings for accelerations
        SelectedAccelerationMap accelerationSettingsMap;
        accelerationSettingsMap[ "Vehicle" ][ "Sun" ].push_back(
                    std::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );

        // Define origin of integration
        std::map< std::string, std::string > centralBodies;
        centralBodies[ "Vehicle" ] = "Earth";
        std::vector< std::string > propagatedBodyList;
        propagatedBodyList.push_back( "Vehicle" );
        std::vector< std::string > centralBodyList;
        centralBodyList.push_back( centralBodies[ "Vehicle" ] );

        // Create accelerations
        AccelerationMap accelerationsMap = createAccelerationModelsMap(
                    bodies, accelerationSettingsMap, centralBodies );

        // Create environment update settings.
        std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                        "Vehicle", centralBodies[ "Vehicle" ], bodies, initialTime ), finalTime );
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                createEnvironmentUpdaterSettings< double >( propagatorSettings, bodies );

        BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 4 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_translational_state_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_translational_state_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( radiation_source_model_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( radiation_source_model_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( cannonball_radiation_pressure_target_model_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( cannonball_radiation_pressure_target_model_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_mass_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_mass_update ).size( ), 1 );


        std::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                createEnvironmentUpdaterForDynamicalEquations< double, double >(
                    propagatorSettings, bodies );

        updater->updateEnvironment( testTime, integratedStateToSet );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodies.at( "Sun" )->getState( ),
                    bodies.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodies.at( "Vehicle" )->getState( ), testState,
                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION(
                    bodies.at( "Vehicle" )->getBodyMass( ), getBodyMass( testTime ),
                    std::numeric_limits< double >::epsilon( ) );


        updater->updateEnvironment(
                    0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ), { translational_state } );
    }

    {
        // Define settings for accelerations
        SelectedAccelerationMap accelerationSettingsMap;
        accelerationSettingsMap[ "Vehicle" ][ "Sun" ].push_back(
                    std::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
        accelerationSettingsMap[ "Vehicle" ][ "Earth" ].push_back(
                    std::make_shared< AccelerationSettings >( aerodynamic ) );

        // Define origin of integration
        std::map< std::string, std::string > centralBodies;
        centralBodies[ "Vehicle" ] = "Earth";
        std::vector< std::string > propagatedBodyList;
        propagatedBodyList.push_back( "Vehicle" );
        std::vector< std::string > centralBodyList;
        centralBodyList.push_back( centralBodies[ "Vehicle" ] );

        // Create accelerations
        AccelerationMap accelerationsMap = createAccelerationModelsMap(
                    bodies, accelerationSettingsMap, centralBodies );


        // Define orientation angles.
        std::shared_ptr< aerodynamics::AtmosphericFlightConditions > vehicleFlightConditions =
                std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodies.at( "Vehicle" )->getFlightConditions( ) );
        double angleOfAttack = 35.0 * mathematical_constants::PI / 180.0;
        double angleOfSideslip = -0.00322;
        double bankAngle = 2.323432;

        std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > vehicleRotationModel =
                createAerodynamicAngleBasedRotationModel(
                                "Vehicle", "Earth", bodies,
                                "ECLIPJ2000", "VehicleFixed" );
        vehicleRotationModel->setAerodynamicAngleFunction(
                    [=](const double ){ return( Eigen::Vector3d( ) << angleOfAttack, angleOfSideslip, bankAngle ).finished( ); } );
        bodies.at( "Vehicle" )->setRotationalEphemeris( vehicleRotationModel );


        std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                        "Vehicle", centralBodies[ "Vehicle" ], bodies, initialTime ), finalTime );
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                createEnvironmentUpdaterSettings< double >( propagatorSettings, bodies );

        // Test update settings
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 6 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_translational_state_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_translational_state_update ).size( ), 2 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( radiation_source_model_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( radiation_source_model_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( cannonball_radiation_pressure_target_model_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( cannonball_radiation_pressure_target_model_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_mass_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_mass_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_rotational_state_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_rotational_state_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( vehicle_flight_conditions_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( vehicle_flight_conditions_update ).size( ), 1 );

        // Create and call updater.
        vehicleRotationModel->setIsBodyInPropagation( true );
        std::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                createEnvironmentUpdaterForDynamicalEquations< double, double >(
                    propagatorSettings, bodies );
        updater->updateEnvironment( testTime, integratedStateToSet );

        // Test if Earth, Sun and Vehicle states are updated.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodies.at( "Earth" )->getState( ),
                    bodies.at( "Earth" )->getEphemeris( )->getCartesianState( testTime ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodies.at( "Sun" )->getState( ),
                    bodies.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodies.at( "Vehicle" )->getState( ), testState,
                    std::numeric_limits< double >::epsilon( ) );

        // Test if Earth rotation is updated.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodies.at( "Earth" )->getCurrentRotationToGlobalFrame( ).toRotationMatrix( ),
                    bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ).toRotationMatrix( ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodies.at( "Earth" )->getCurrentRotationToLocalFrame( ).toRotationMatrix( ),
                    bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToTargetFrame( testTime ).toRotationMatrix( ),
                    std::numeric_limits< double >::epsilon( ) );

        // Test if body mass is updated
        BOOST_CHECK_CLOSE_FRACTION(
                    bodies.at( "Vehicle" )->getBodyMass( ), getBodyMass( testTime ),
                    std::numeric_limits< double >::epsilon( ) );

        // Check if flight conditions update has been called
        BOOST_CHECK_EQUAL(
                    ( vehicleFlightConditions->getCurrentAirspeed( ) == vehicleFlightConditions->getCurrentAirspeed( ) ), 1 );
        BOOST_CHECK_EQUAL(
                    ( vehicleFlightConditions->getCurrentAltitude( ) == vehicleFlightConditions->getCurrentAltitude( ) ), 1 );
        BOOST_CHECK_EQUAL(
                    ( vehicleFlightConditions->getCurrentDensity( ) == vehicleFlightConditions->getCurrentDensity( ) ), 1 );
        BOOST_CHECK_EQUAL( vehicleFlightConditions->getCurrentTime( ), testTime );

        // Check if radiation pressure update is updated.
        std::shared_ptr< electromagnetism::RadiationPressureTargetModel > radiationPressureInterface =
                bodies.at( "Vehicle" )->getRadiationPressureTargetModel( );
//        BOOST_CHECK_EQUAL(
//                    ( radiationPressureInterface->getCurrentTime( ) == radiationPressureInterface->getCurrentTime( ) ), 1 );
//        BOOST_CHECK_EQUAL(
//                    ( radiationPressureInterface->getRadiationPressure( ) ==
//                radiationPressureInterface->getRadiationPressure( ) ), 1 );
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                    radiationPressureInterface->getCurrentSolarVector( ),
//                    ( bodies.at( "Sun" )->getPosition( ) - bodies.at( "Vehicle" )->getPosition( ) ),
//                    std::numeric_limits< double >::epsilon( ) );

        updater->updateEnvironment(
                    0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ), { translational_state } );


    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat


