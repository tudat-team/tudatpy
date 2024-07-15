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

#include <memory>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/simulation.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat;
using namespace tudat::propagators;
using namespace tudat::simulation_setup;
using namespace tudat::numerical_integrators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;

BOOST_AUTO_TEST_SUITE( test_propagation_radiation_pressure )

// Test custom state propagation, linearly decreasing with time
BOOST_AUTO_TEST_CASE( testMultiTypeRadiationPressurePropagation )
{
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::map< double, Eigen::VectorXd > > singleBodyNumericalResults;
    std::vector< std::map< double, Eigen::VectorXd > > singleBodyDependentVariableResults;

    for( unsigned int multiBodyTest = 0; multiBodyTest < 3; multiBodyTest++ )
    {
        for( unsigned int test = 0; test < 6; test++ )
        {
            std::vector< std::string > bodyNames = { "Earth", "Moon", "Sun" };

            // Specify initial time
            double initialEphemerisTime = 1.0E7;
            double finalEphemerisTime = initialEphemerisTime + 3.0 * 3600.0;

            // Get initial state vector as input to integration.
            // Set Keplerian elements for Capsule.
            Eigen::Vector6d initialStateInKeplerianElements1;
            initialStateInKeplerianElements1( semiMajorAxisIndex ) = spice_interface::getAverageRadius( "Earth" ) + 800.0E3;
            initialStateInKeplerianElements1( eccentricityIndex ) = 0.005;
            initialStateInKeplerianElements1( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
            initialStateInKeplerianElements1( argumentOfPeriapsisIndex )
                = unit_conversions::convertDegreesToRadians( 235.7 );
            initialStateInKeplerianElements1( longitudeOfAscendingNodeIndex )
                = unit_conversions::convertDegreesToRadians( 23.4 );
            initialStateInKeplerianElements1( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

            Eigen::Vector6d initialStateInKeplerianElements2;
            initialStateInKeplerianElements2( semiMajorAxisIndex ) = spice_interface::getAverageRadius( "Earth" ) + 1200.0E3;
            initialStateInKeplerianElements2( eccentricityIndex ) = 0.001;
            initialStateInKeplerianElements2( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 35.3 );
            initialStateInKeplerianElements2( argumentOfPeriapsisIndex )
                = unit_conversions::convertDegreesToRadians( 235.7 );
            initialStateInKeplerianElements2( longitudeOfAscendingNodeIndex )
                = unit_conversions::convertDegreesToRadians( 93.4 );
            initialStateInKeplerianElements2( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 239.87 );


            // Create bodies needed in simulation
            BodyListSettings bodySettings = getDefaultBodySettings(
                bodyNames, "Earth", "ECLIPJ2000" );
            std::vector< std::string > vehiclesToAdd = {"Vehicle"};
            if( multiBodyTest > 0 )
            {
                vehiclesToAdd.push_back( "Vehicle2" );
            }

            double bodyMass = 300.0;
            for( unsigned int i = 0; i < vehiclesToAdd.size( ); i++ )
            {
                bodySettings.addSettings( vehiclesToAdd.at( i ) );
                bodySettings.get( vehiclesToAdd.at( i ) )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Earth", "ECLIPJ2000", vehiclesToAdd.at( i ) + "_Fixed" );
                bodySettings.get( vehiclesToAdd.at( i ) )->bodyExteriorPanelSettings_ = bodyWingPanelledGeometry( 2.0, 1.0, 4.0, 20.0, 0.0, 0.0, 0.0, 0.0, false, false );
                if( vehiclesToAdd.at( i ) == "Vehicle" )
                {
                    bodySettings.get( vehiclesToAdd.at( i ))->ephemerisSettings = keplerEphemerisSettings(
                        initialStateInKeplerianElements1, initialEphemerisTime,
                        spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );
                }
                else
                {
                    bodySettings.get( vehiclesToAdd.at( i ))->ephemerisSettings = keplerEphemerisSettings(
                        initialStateInKeplerianElements2, initialEphemerisTime,
                        spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );
                }
                bodySettings.get( vehiclesToAdd.at( i ) )->rigidBodyPropertiesSettings = constantRigidBodyPropertiesSettings( bodyMass );
            }


            bool usePaneledEarth = true;
            if( test < 3 )
            {
                bodySettings.get("Earth")->radiationSourceModelSettings =
                    extendedRadiationSourceModelSettings( {
                        albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke, "Sun"),
                        delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke, "Sun" ) }, {6, 12}, {"Moon"} );
            }
            else
            {
                usePaneledEarth = false;
                bodySettings.get("Earth")->radiationSourceModelSettings =
                    isotropicPointRadiationSourceModelSettings( constantLuminosityModelSettings( 50000000000000000.) );

            }

            bodySettings.at("Moon")->radiationSourceModelSettings =
                extendedRadiationSourceModelSettings({
                    albedoPanelRadiosityModelSettings(SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1, "Sun"),
                    angleBasedThermalPanelRadiosityModelSettings(95, 385, 0.95, "Sun")
                }, {6, 12 }, {"Earth"});
            // Crate bodies
            std::map<std::string, std::vector<std::string>> sourceToTargetOccultingBodies;
            sourceToTargetOccultingBodies[ "Sun" ] = { "Earth" };
            sourceToTargetOccultingBodies[ "Moon" ] = { "Earth" };

            SystemOfBodies bodies = createSystemOfBodies( bodySettings );
            for( unsigned int i = 0; i < vehiclesToAdd.size( ); i++ )
            {
                addRadiationPressureTargetModel( bodies, vehiclesToAdd.at( i ), paneledRadiationPressureTargetModelSettingsWithOccultationMap( sourceToTargetOccultingBodies ) );
                addRadiationPressureTargetModel( bodies, vehiclesToAdd.at( i ), cannonballRadiationPressureTargetModelSettingsWithOccultationMap( 3.0, 1.5, sourceToTargetOccultingBodies ) );
            }

            // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
            SelectedAccelerationMap accelerationMap;
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
            bool usePaneledVehicleForEarth = true;
            bool usePaneledVehicleForSun = true;

            if( test % 3 == 0 )
            {
                accelerationsOfVehicle[ "Sun" ].push_back( radiationPressureAcceleration( paneled_target ));
                accelerationsOfVehicle[ "Earth" ].push_back( radiationPressureAcceleration( paneled_target ));
            }
            else if( test % 3 == 1 )
            {
                usePaneledVehicleForEarth = false;
                usePaneledVehicleForSun = false;
                accelerationsOfVehicle[ "Sun" ].push_back( radiationPressureAcceleration( cannonball_target ));
                accelerationsOfVehicle[ "Earth" ].push_back( radiationPressureAcceleration( cannonball_target ));
            }
            else if( test % 3 == 2 )
            {
                usePaneledVehicleForEarth = false;
                accelerationsOfVehicle[ "Sun" ].push_back( radiationPressureAcceleration( paneled_target ));
                accelerationsOfVehicle[ "Earth" ].push_back( radiationPressureAcceleration( cannonball_target ));
            }
            accelerationsOfVehicle[ "Earth" ].push_back( pointMassGravityAcceleration(  ) );
            std::vector< std::string > centralBodies;
            for( unsigned int i = 0; i < vehiclesToAdd.size( ); i++ )
            {
                accelerationMap[ vehiclesToAdd.at( i ) ] = accelerationsOfVehicle;
                centralBodies.push_back( "Earth" );
            }
            // Define list of bodies to propagate
            std::vector< std::string > bodiesToIntegrate = vehiclesToAdd;
            if( multiBodyTest == 2 )
            {
                std::reverse(bodiesToIntegrate.begin(),bodiesToIntegrate.end());
            }

            // Create acceleration models and propagation settings.
            AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );

            // Define numerical integrator settings.
            std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, initialEphemerisTime, 60.0 );


            // Convert apollo state from Keplerian elements to Cartesian elements.
            Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                bodiesToIntegrate, centralBodies, bodies, initialEphemerisTime );

            // Create dependent variables
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
            dependentVariablesList.push_back( singleAccelerationDependentVariable( radiation_pressure, "Vehicle", "Sun" ) );
            dependentVariablesList.push_back( singleAccelerationDependentVariable( radiation_pressure, "Vehicle", "Earth" ) );

            dependentVariablesList.push_back( vehiclePanelInertialSurfaceNormals( "Vehicle" ) );
            dependentVariablesList.push_back( vehiclePanelBodyFixedSurfaceNormals( "Vehicle" ) );
            if( usePaneledVehicleForSun )
            {
                dependentVariablesList.push_back( perVehiclePanelRadiationPressureForce( "Vehicle", "Sun" ));
            }
            dependentVariablesList.push_back( relativePositionDependentVariable( "Vehicle", "Sun" ) );
            dependentVariablesList.push_back( inertialToBodyFixedRotationMatrixVariable( "Vehicle" ) );
            if( usePaneledVehicleForEarth )
            {
                dependentVariablesList.push_back( perVehiclePanelRadiationPressureForce( "Vehicle", "Earth" ));
            }
            if( usePaneledEarth )
            {
                dependentVariablesList.push_back( paneledRadiationSourcePerPanelIrradiance( "Vehicle", "Earth" ));
                dependentVariablesList.push_back( paneledRadiationSourceGeometry( "Vehicle", "Earth" ));
            }
            dependentVariablesList.push_back( relativePositionDependentVariable( "Vehicle", "Earth" ) );


            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, initialEphemerisTime, integratorSettings,
                      std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ), cowell, dependentVariablesList );

            propagatorSettings->getPrintSettings( )->setPrintDependentVariableData( true );

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, propagatorSettings );

            std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );
            std::shared_ptr< SingleArcDependentVariablesInterface< double > >  dependentVariablesInterface =
                dynamicsSimulator.getSingleArcPropagationResults( )->getSingleArcDependentVariablesInterface( );
            std::pair< int, int > inertialPanelNormalIndices = dependentVariablesInterface->getSingleDependentVariableIndices( vehiclePanelInertialSurfaceNormals( "Vehicle" ) );
            std::pair< int, int > bodyFixedPanelNormalIndices = dependentVariablesInterface->getSingleDependentVariableIndices( vehiclePanelBodyFixedSurfaceNormals( "Vehicle" ) );
            std::pair< int, int > inertialToBodyFixedIndices = dependentVariablesInterface->getSingleDependentVariableIndices( inertialToBodyFixedRotationMatrixVariable( "Vehicle" ) );
            std::pair< int, int > vehicleSunRelativePositionIndices = dependentVariablesInterface->getSingleDependentVariableIndices( relativePositionDependentVariable( "Vehicle", "Sun" ) );

            for( auto it : dependentVariableHistory )
            {
                Eigen::VectorXd currentBodyFixedNormals = it.second.segment( bodyFixedPanelNormalIndices.first, bodyFixedPanelNormalIndices.second );
                Eigen::VectorXd currentInertialNormals = it.second.segment( inertialPanelNormalIndices.first, inertialPanelNormalIndices.second );

                BOOST_CHECK_EQUAL( currentBodyFixedNormals.rows( ), currentInertialNormals.rows( ) );
                BOOST_CHECK_EQUAL( currentBodyFixedNormals.rows( ), 8 * 3 );
                for( unsigned int i = 0; i < 6; i++ )
                {
                    Eigen::Vector3d currentNormal = currentBodyFixedNormals.segment( 3 * i, 3 );
                    Eigen::Vector3d expectedNormal;
                    switch( i )
                    {
                    case 0:
                        expectedNormal = Eigen::Vector3d::UnitX( );
                        break;
                    case 1:
                        expectedNormal = -Eigen::Vector3d::UnitX( );
                        break;
                    case 2:
                        expectedNormal = Eigen::Vector3d::UnitY( );
                        break;
                    case 3:
                        expectedNormal = -Eigen::Vector3d::UnitY( );
                        break;
                    case 4:
                        expectedNormal = Eigen::Vector3d::UnitZ( );
                        break;
                    case 5:
                        expectedNormal = -Eigen::Vector3d::UnitZ( );
                        break;
                    }
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( expectedNormal( j ) - currentNormal( j ) ), 1.0E-15 );
                    }
                }

                Eigen::Matrix3d currentInertialToBodyFixedMatrix  = getMatrixFromVectorRotationRepresentation(
                    it.second.segment( inertialToBodyFixedIndices.first, inertialToBodyFixedIndices.second ) );
                Eigen::VectorXd recomputedCurrentBodyFixedNormals = Eigen::VectorXd( currentInertialNormals.rows( ) );
                for( unsigned int i = 0; i < currentInertialNormals.rows( ) / 3; i++ )
                {
                    recomputedCurrentBodyFixedNormals.segment( 3 * i, 3 ) = currentInertialToBodyFixedMatrix * currentInertialNormals.segment( 3 * i, 3 );
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( recomputedCurrentBodyFixedNormals( 3 * i + j ) - currentBodyFixedNormals( 3 * i + j ) ), 1.0E-14 );
                    }
                }

                Eigen::Vector3d directionFromSun = it.second.segment( vehicleSunRelativePositionIndices.first, vehicleSunRelativePositionIndices.second ).normalized( );
                for( unsigned int i = 6; i < 8; i++ )
                {
                    Eigen::Vector3d inertialSolarPanelNormal = currentInertialNormals.segment( 3 * i, 3 );
                    if( i == 6 )
                    {
                        inertialSolarPanelNormal *= -1.0;
                    }
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( inertialSolarPanelNormal( j ) - directionFromSun( j ) ), 1.0E-14 );
                    }
                }

                if( usePaneledVehicleForSun )
                {
                    std::pair< int, int > perPanelSunForceIndices = dependentVariablesInterface->getSingleDependentVariableIndices( perVehiclePanelRadiationPressureForce( "Vehicle", "Sun" ) );
                    Eigen::VectorXd perPanelSunForce = it.second.segment( perPanelSunForceIndices.first, perPanelSunForceIndices.second );
                    BOOST_CHECK_EQUAL( perPanelSunForce.rows( ), 8 * 3 );

                    Eigen::Vector3d totalSunForce = Eigen::Vector3d::Zero( );
                    for( unsigned int i = 0; i < 8; i++ )
                    {
                        Eigen::Vector3d currentPanelForce = perPanelSunForce.segment( 3 * i, 3 );
                        Eigen::Vector3d currentPanelForceNormalized = perPanelSunForce.segment( 3 * i, 3 ).normalized( );
                        Eigen::Vector3d inertialPanelForce = currentInertialToBodyFixedMatrix.transpose( ) * currentPanelForce;
                        Eigen::Vector3d inertialPanelForceNormalized = currentInertialToBodyFixedMatrix.transpose( ) * currentPanelForceNormalized;

                        totalSunForce += inertialPanelForce;
                        if( currentPanelForce.norm( ) > 0.0 )
                        {
                            for ( unsigned int j = 0; j < 3; j++ )
                            {
                                BOOST_CHECK_SMALL( std::fabs( inertialPanelForceNormalized( j ) - directionFromSun( j )), 1.0E-14 );
                            }
                        }
                    }
                    Eigen::Vector3d totalSunAcceleration = it.second.segment( 0, 3 );
                    Eigen::Vector3d recomputedSunAcceleration = totalSunForce / bodyMass;
                    for ( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( totalSunAcceleration( j ) - recomputedSunAcceleration( j )), ( 1.0E-14 * totalSunAcceleration.norm( ) ) );
                    }
                }

                if( usePaneledVehicleForEarth )
                {
                    std::pair< int, int > perPanelEarthForceIndices = dependentVariablesInterface->getSingleDependentVariableIndices( perVehiclePanelRadiationPressureForce( "Vehicle", "Earth" ) );
                    Eigen::VectorXd perPanelEarthForce = it.second.segment( perPanelEarthForceIndices.first, perPanelEarthForceIndices.second );
                    BOOST_CHECK_EQUAL( perPanelEarthForce.rows( ), 8 * 3 );

                    Eigen::Vector3d totalEarthForce = Eigen::Vector3d::Zero( );
                    for( unsigned int i = 0; i < 8; i++ )
                    {
                        Eigen::Vector3d currentPanelForce = perPanelEarthForce.segment( 3 * i, 3 );
                        Eigen::Vector3d inertialPanelForce = currentInertialToBodyFixedMatrix.transpose( ) * currentPanelForce;
            //            Eigen::Vector3d inertialPanelForceNormalized = currentInertialToBodyFixedMatrix.transpose( ) * currentPanelForceNormalized;

                        totalEarthForce += inertialPanelForce;
                    }

                    Eigen::Vector3d totalEarthAcceleration = it.second.segment( 3, 3 );
                    Eigen::Vector3d recomputedEarthAcceleration = totalEarthForce / bodyMass;
                    for ( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( totalEarthAcceleration( j ) - recomputedEarthAcceleration( j )), ( 1.0E-14 * totalEarthAcceleration.norm( ) ) );
                    }
                }


            }
            if( multiBodyTest == 0 )
            {
                singleBodyNumericalResults.push_back( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ) );
                singleBodyDependentVariableResults.push_back( dynamicsSimulator.getDependentVariableHistory( ) );
            }
            else
            {
                auto dependentVariables = dynamicsSimulator.getDependentVariableHistory( );

                BOOST_CHECK_EQUAL( dependentVariables.size( ), singleBodyDependentVariableResults.at( test ).size( ) );
                for( auto it : singleBodyDependentVariableResults.at( test ) )
                {
                    BOOST_CHECK_EQUAL( dependentVariables.count( it.first ), 1 );
                    Eigen::VectorXd testVector1 = it.second;
                    Eigen::VectorXd testVector2 = dependentVariables.at( it.first );
                    BOOST_CHECK_EQUAL( testVector1.rows( ), testVector2.rows( ) );

//                    std::cout<<( testVector1 - testVector2 ).transpose( )<<std::endl;
                    for( unsigned int j = 0; j < testVector1.rows( ); j++ )
                    {
                        if( testVector1( j ) == 0.0 )
                        {
                            BOOST_CHECK_EQUAL( testVector2( j ), 0.0 );
                        }
                        else
                        {
                            BOOST_CHECK_CLOSE_FRACTION( testVector1( j ), testVector2( j ), std::numeric_limits< double >::epsilon( ) * 10.0 );
                        }
                    }
                }

                auto numericalResults = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

                BOOST_CHECK_EQUAL( numericalResults.size( ), singleBodyNumericalResults.at( test ).size( ) );
                for( auto it : singleBodyNumericalResults.at( test ) )
                {
                    Eigen::VectorXd testVector1 = it.second;
                    Eigen::VectorXd testVector2 = numericalResults.at( it.first );
//                    std::cout<<( testVector1  ).transpose( )<<std::endl;
//                    std::cout<<( testVector2  ).transpose( )<<std::endl<<std::endl;

                    if( multiBodyTest == 1 )
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testVector1, testVector2.segment( 0, 6 ), ( std::numeric_limits< double >::epsilon( ) * 10.0 ) );
                    }
                    else
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testVector1, testVector2.segment( 6, 6 ), ( std::numeric_limits< double >::epsilon( ) * 10.0 ) );
                    }
                }


            }
            //
        //    input_output::writeDataMapToTextFile(
        //        dynamicsSimulator.getDependentVariableHistory( ),
        //        "radiationPressureDependentVariables.dat",
        //        "/home/dominic/Downloads/",
        //        "",
        //        std::numeric_limits< double >::digits10,
        //        std::numeric_limits< double >::digits10,
        //        "," );

        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
