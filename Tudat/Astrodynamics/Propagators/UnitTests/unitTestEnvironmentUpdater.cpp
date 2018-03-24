/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/assign/list_of.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h"

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
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = getDefaultSingleBodySettings( "Earth", 0.0, 10.0 * 86400.0 );
    bodySettings[ "Sun" ] = getDefaultSingleBodySettings( "Sun", 0.0,10.0 * 86400.0 );
    bodySettings[ "Moon" ] = getDefaultSingleBodySettings( "Moon", 0.0,10.0 * 86400.0 );
    bodySettings[ "Mars" ] = getDefaultSingleBodySettings( "Mars", 0.0,10.0 * 86400.0 );
    bodySettings[ "Venus" ] = getDefaultSingleBodySettings( "Venus", 0.0,10.0 * 86400.0 );

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
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                gravitationalParameter, 6378.0E3, cosineCoefficients, sineCoefficients,
                "IAU_Earth" );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

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
        integratedStateToSet[ transational_state ] = testState;
        testTime = 2.0 * 86400.0;

        // Test if environment is updated for inly central gravity accelerations
        {
            // Define accelerations
            accelerationSettingsMap[ "Moon" ][ "Sun" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
            accelerationSettingsMap[ "Moon" ][ "Earth" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );

            // Define origin of integration to be barycenter.
            centralBodies[ "Moon" ] = "SSB";
            propagatedBodyList.push_back( "Moon" );
            centralBodyList.push_back( centralBodies[ "Moon" ] );

            // Create accelerations
            AccelerationMap accelerationsMap = createAccelerationModelsMap(
                        bodyMap, accelerationSettingsMap, centralBodies );

            // Create environment update settings.
            boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                            "Moon", centralBodies[ "Moon" ], bodyMap, initialTime ), finalTime );
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                    createEnvironmentUpdaterSettings< double >( propagatorSettings, bodyMap );

            // Test update settings
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_transational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_transational_state_update ).size( ), 2 );

            // Create and call updater.
            boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                    createEnvironmentUpdaterForDynamicalEquations< double, double >(
                        propagatorSettings, bodyMap );
            updater->updateEnvironment( testTime, integratedStateToSet );

            // Test if Earth, Sun and Moon are updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getState( ),
                        bodyMap.at( "Earth" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Sun" )->getState( ),
                        bodyMap.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Moon" )->getState( ), testState,
                        std::numeric_limits< double >::epsilon( ) );

            // Test if Mars is not updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getState( ), Eigen::Vector6d::Zero( ),
                        std::numeric_limits< double >::epsilon( ) );

            // Update environment to new time, and state from environment.
            updater->updateEnvironment(
                        0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ),
                        boost::assign::list_of( transational_state ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getState( ),
                        bodyMap.at( "Earth" )->getEphemeris( )->getCartesianState( 0.5 * testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Sun" )->getState( ),
                        bodyMap.at( "Sun" )->getEphemeris( )->getCartesianState( 0.5 * testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Moon" )->getState( ),
                        bodyMap.at( "Moon" )->getEphemeris( )->getCartesianState( 0.5 * testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getState( ), Eigen::Vector6d::Zero( ),
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
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
            accelerationSettingsMap[ "Moon" ][ "Mars" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );

            // Define origin of integration
            centralBodies[ "Moon" ] = "Earth";
            propagatedBodyList.push_back( "Moon" );
            centralBodyList.push_back( centralBodies[ "Moon" ] );

            // Create accelerations
            AccelerationMap accelerationsMap = createAccelerationModelsMap(
                        bodyMap, accelerationSettingsMap, centralBodies );

            // Create environment update settings.
            boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                            "Moon", centralBodies[ "Moon" ], bodyMap, initialTime ), finalTime );
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                    createEnvironmentUpdaterSettings< double >( propagatorSettings, bodyMap );

            // Test update settings
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_transational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_transational_state_update ).size( ), 3 );

            // Create and call updater.
            boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                    createEnvironmentUpdaterForDynamicalEquations< double, double >(
                        propagatorSettings, bodyMap );
            updater->updateEnvironment( testTime, integratedStateToSet );

            // Test if Earth, Sun, Mars and Moon are updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getState( ),
                        bodyMap.at( "Earth" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Sun" )->getState( ),
                        bodyMap.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Moon" )->getState( ), testState,
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getState( ),
                        bodyMap.at( "Mars" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );

            // Test if Venus is not updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Venus" )->getState( ), Eigen::Vector6d::Zero( ),
                        std::numeric_limits< double >::epsilon( ) );

            // Update environment to new time, and state from environment.
            updater->updateEnvironment(
                        0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ),
                        boost::assign::list_of( transational_state ) );

        }

        // Test spherical harmonic acceleration update
        {
            accelerationSettingsMap.clear( );
            centralBodies.clear( );
            propagatedBodyList.clear( );
            centralBodyList.clear( );

            // Set acceleration models.
            accelerationSettingsMap[ "Moon" ][ "Sun" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
            accelerationSettingsMap[ "Moon" ][ "Earth" ].push_back(
                        boost::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );
            accelerationSettingsMap[ "Moon" ][ "Mars" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );

            // Define origin of integration
            centralBodies[ "Moon" ] = "Earth";
            propagatedBodyList.push_back( "Moon" );
            centralBodyList.push_back( centralBodies[ "Moon" ] );

            // Create accelerations
            AccelerationMap accelerationsMap = createAccelerationModelsMap(
                        bodyMap, accelerationSettingsMap, centralBodies );

            // Create environment update settings.
            boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                            "Moon", centralBodies[ "Moon" ], bodyMap, initialTime ), finalTime );
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                    createEnvironmentUpdaterSettings< double >( propagatorSettings, bodyMap );

            // Test update settings
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 3 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_transational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_transational_state_update ).size( ), 3 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_rotational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_rotational_state_update ).size( ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( spherical_harmonic_gravity_field_update ), 1 );

            // Create and call updater.
            boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                    createEnvironmentUpdaterForDynamicalEquations< double, double >(
                        propagatorSettings, bodyMap );
            updater->updateEnvironment( testTime, integratedStateToSet );

            // Test if Earth, Sun, Mars and Moon are updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getState( ),
                        bodyMap.at( "Earth" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Sun" )->getState( ),
                        bodyMap.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Moon" )->getState( ), testState,
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getState( ),
                        bodyMap.at( "Mars" )->getEphemeris( )->getCartesianState( testTime ),
                        std::numeric_limits< double >::epsilon( ) );

            // Test if Venus is not updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Venus" )->getState( ), Eigen::Vector6d::Zero( ),
                        std::numeric_limits< double >::epsilon( ) );

            // Test if Earth rotation is updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getCurrentRotationToGlobalFrame( ).toRotationMatrix( ),
                        bodyMap.at( "Earth" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ).toRotationMatrix( ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getCurrentRotationToLocalFrame( ).toRotationMatrix( ),
                        bodyMap.at( "Earth" )->getRotationalEphemeris( )->getRotationToTargetFrame( testTime ).toRotationMatrix( ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getCurrentRotationMatrixDerivativeToGlobalFrame( ),
                        bodyMap.at( "Earth" )->getRotationalEphemeris( )->getDerivativeOfRotationToBaseFrame( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getCurrentRotationMatrixDerivativeToLocalFrame( ),
                        bodyMap.at( "Earth" )->getRotationalEphemeris( )->getDerivativeOfRotationToTargetFrame( testTime ),
                        std::numeric_limits< double >::epsilon( ) );

            // Test if Mars rotation is not updated
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getCurrentRotationToLocalFrame( ).toRotationMatrix( ),
                        Eigen::Matrix3d::Identity( ), std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getCurrentRotationToGlobalFrame( ).toRotationMatrix( ),
                        Eigen::Matrix3d::Identity( ), std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getCurrentRotationMatrixDerivativeToGlobalFrame( ),
                        Eigen::Matrix3d::Zero( ), std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getCurrentRotationMatrixDerivativeToLocalFrame( ),
                        Eigen::Matrix3d::Zero( ), std::numeric_limits< double >::epsilon( ) );

            // Update environment to new time, and state from environment.
            updater->updateEnvironment(
                        0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ),
                        boost::assign::list_of( transational_state ) );

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
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = getDefaultSingleBodySettings( "Earth", 0.0, 10.0 * 86400.0 );
    bodySettings[ "Sun" ] = getDefaultSingleBodySettings( "Sun", 0.0,10.0 * 86400.0 );

    // Get settings for vehicle
    double area = 2.34;
    double coefficient = 1.2;
    bodySettings[ "Vehicle" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Vehicle" ]->radiationPressureSettings[ "Sun" ] =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >( "Sun", area, coefficient );
    bodySettings[ "Vehicle" ]->ephemerisSettings =
            boost::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( ) << 7000.0E3, 0.05, 0.3, 0.0, 0.0, 0.0 ).finished( ),
                0.0, spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( bodySettings );
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                getApolloCoefficientInterface( ) );
    bodyMap[ "Vehicle" ]->setBodyMassFunction( &getBodyMass );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Define test time and state.
    double testTime = 2.0 * 86400.0;
    std::unordered_map< IntegratedStateType, Eigen::VectorXd > integratedStateToSet;
    Eigen::VectorXd testState = 1.1 * bodyMap[ "Vehicle" ]->getEphemeris( )->getCartesianState( testTime ) +
            bodyMap.at( "Earth" )->getEphemeris( )->getCartesianState( testTime );
    integratedStateToSet[ transational_state ] = testState;

    {
        // Define settings for accelerations
        SelectedAccelerationMap accelerationSettingsMap;
        accelerationSettingsMap[ "Vehicle" ][ "Sun" ].push_back(
                    boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );

        // Define origin of integration
        std::map< std::string, std::string > centralBodies;
        centralBodies[ "Vehicle" ] = "Earth";
        std::vector< std::string > propagatedBodyList;
        propagatedBodyList.push_back( "Vehicle" );
        std::vector< std::string > centralBodyList;
        centralBodyList.push_back( centralBodies[ "Vehicle" ] );

        // Create accelerations
        AccelerationMap accelerationsMap = createAccelerationModelsMap(
                    bodyMap, accelerationSettingsMap, centralBodies );

        // Create environment update settings.
        boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                        "Vehicle", centralBodies[ "Vehicle" ], bodyMap, initialTime ), finalTime );
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                createEnvironmentUpdaterSettings< double >( propagatorSettings, bodyMap );

        BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 3 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_transational_state_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_transational_state_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( radiation_pressure_interface_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( radiation_pressure_interface_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_mass_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_mass_update ).size( ), 1 );


        boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                createEnvironmentUpdaterForDynamicalEquations< double, double >(
                    propagatorSettings, bodyMap );

        updater->updateEnvironment( testTime, integratedStateToSet );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodyMap.at( "Sun" )->getState( ),
                    bodyMap.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodyMap.at( "Vehicle" )->getState( ), testState,
                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION(
                    bodyMap.at( "Vehicle" )->getBodyMass( ), getBodyMass( testTime ),
                    std::numeric_limits< double >::epsilon( ) );


        updater->updateEnvironment(
                    0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ),
                    boost::assign::list_of( transational_state ) );
    }

    {
        // Define settings for accelerations
        SelectedAccelerationMap accelerationSettingsMap;
        accelerationSettingsMap[ "Vehicle" ][ "Sun" ].push_back(
                    boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
        accelerationSettingsMap[ "Vehicle" ][ "Earth" ].push_back(
                    boost::make_shared< AccelerationSettings >( aerodynamic ) );

        // Define origin of integration
        std::map< std::string, std::string > centralBodies;
        centralBodies[ "Vehicle" ] = "Earth";
        std::vector< std::string > propagatedBodyList;
        propagatedBodyList.push_back( "Vehicle" );
        std::vector< std::string > centralBodyList;
        centralBodyList.push_back( centralBodies[ "Vehicle" ] );

        // Create accelerations
        AccelerationMap accelerationsMap = createAccelerationModelsMap(
                    bodyMap, accelerationSettingsMap, centralBodies );


        // Define orientation angles.
        boost::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions =
                bodyMap[ "Vehicle" ]->getFlightConditions( );
        double angleOfAttack = 35.0 * mathematical_constants::PI / 180.0;
        double angleOfSideslip = -0.00322;
        double bankAngle = 2.323432;
        vehicleFlightConditions->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                    boost::lambda::constant( angleOfAttack ),
                    boost::lambda::constant( angleOfSideslip ),
                    boost::lambda::constant( bankAngle ) );

        boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                        "Vehicle", centralBodies[ "Vehicle" ], bodyMap, initialTime ), finalTime );
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                createEnvironmentUpdaterSettings< double >( propagatorSettings, bodyMap );

        // Test update settings
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 5 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_transational_state_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_transational_state_update ).size( ), 2 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( radiation_pressure_interface_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( radiation_pressure_interface_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_mass_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_mass_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_rotational_state_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_rotational_state_update ).size( ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( vehicle_flight_conditions_update ), 1 );
        BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( vehicle_flight_conditions_update ).size( ), 1 );

        // Create and call updater.
        boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                createEnvironmentUpdaterForDynamicalEquations< double, double >(
                    propagatorSettings, bodyMap );
        updater->updateEnvironment( testTime, integratedStateToSet );

        // Test if Earth, Sun and Vehicle states are updated.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodyMap.at( "Earth" )->getState( ),
                    bodyMap.at( "Earth" )->getEphemeris( )->getCartesianState( testTime ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodyMap.at( "Sun" )->getState( ),
                    bodyMap.at( "Sun" )->getEphemeris( )->getCartesianState( testTime ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodyMap.at( "Vehicle" )->getState( ), testState,
                    std::numeric_limits< double >::epsilon( ) );

        // Test if Earth rotation is updated.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodyMap.at( "Earth" )->getCurrentRotationToGlobalFrame( ).toRotationMatrix( ),
                    bodyMap.at( "Earth" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ).toRotationMatrix( ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    bodyMap.at( "Earth" )->getCurrentRotationToLocalFrame( ).toRotationMatrix( ),
                    bodyMap.at( "Earth" )->getRotationalEphemeris( )->getRotationToTargetFrame( testTime ).toRotationMatrix( ),
                    std::numeric_limits< double >::epsilon( ) );

        // Test if body mass is updated
        BOOST_CHECK_CLOSE_FRACTION(
                    bodyMap.at( "Vehicle" )->getBodyMass( ), getBodyMass( testTime ),
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
        boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface =
                bodyMap.at( "Vehicle" )->getRadiationPressureInterfaces( ).at( "Sun" );
        BOOST_CHECK_EQUAL(
                    ( radiationPressureInterface->getCurrentTime( ) == radiationPressureInterface->getCurrentTime( ) ), 1 );
        BOOST_CHECK_EQUAL(
                    ( radiationPressureInterface->getCurrentRadiationPressure( ) ==
                radiationPressureInterface->getCurrentRadiationPressure( ) ), 1 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    radiationPressureInterface->getCurrentSolarVector( ),
                    ( bodyMap.at( "Sun" )->getPosition( ) - bodyMap.at( "Vehicle" )->getPosition( ) ),
                    std::numeric_limits< double >::epsilon( ) );

        updater->updateEnvironment(
                    0.5 * testTime, std::unordered_map< IntegratedStateType, Eigen::VectorXd >( ),
                    boost::assign::list_of( transational_state ) );


    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat


