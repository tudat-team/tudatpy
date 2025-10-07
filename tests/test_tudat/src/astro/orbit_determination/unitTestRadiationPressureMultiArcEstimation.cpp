/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"

//namespace tudat
//{
//
//namespace unit_tests
//{

//BOOST_AUTO_TEST_SUITE( test_radiation_pressure_estimation )

// Using declarations.
using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::ground_stations;
using namespace tudat::observation_models;
using namespace tudat::ephemerides;
using namespace tudat::electromagnetism;
using namespace tudat;


int main( )
//BOOST_AUTO_TEST_CASE( test_PanelledRadiationPressureEstimation )
{
    spice_interface::loadStandardSpiceKernels( );

    for( int testCase = 0; testCase < 3; testCase++ )
    {
        // Define bodies in simulation
        std::vector< std::string > bodyNames;
        bodyNames.push_back( "Earth" );
        bodyNames.push_back( "Sun" );

        // Specify initial time
        double initialEphemerisTime = double( 1.0E7 );
        double finalEphemerisTime = initialEphemerisTime + 0.5 * 86400.0;

        // Create bodies needed in simulation
        BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - 3600.0, finalEphemerisTime + 3600.0 );
        setSimpleRotationSettingsFromSpice( bodySettings, "Earth", initialEphemerisTime );

        // Add Earth radiation pressure models
        bodySettings.at( "Earth" )->radiationSourceModelSettings =
                                        extendedRadiationSourceModelSettings({
                                            albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke, "Sun"),
                                            delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke, "Sun")
                                        }, {6, 12});


        bodySettings.addSettings( "Vehicle" );
        bodySettings.at( "Vehicle" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "Earth", "ECLIPJ2000" );
        bodySettings.at( "Vehicle" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Create spacecraft object.
        bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

        Eigen::Vector7d rotationalStateVehicle;
        rotationalStateVehicle.segment( 0, 4 ) =
                linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
        rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero( );
        bodies.at( "Vehicle" )
                ->setRotationalEphemeris(
                        std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000", "VehicleFixed" ) );
        {
            std::shared_ptr< FullPanelledBodySettings > panelSettings;
            if( testCase < 2 )
            {
                std::vector< std::shared_ptr< BodyPanelSettings > > panelSettingsList;

                double specular1 = 0.35, specular2 = 0.35, diffuse1 = 0.25, diffuse2 = 0.25;
                if( testCase == 1 )
                {
                    specular1 += 0.1;
                    specular2 -= 0.1;
                    diffuse1 += 0.05;
                    diffuse2 -= 0.05;
                }
                panelSettingsList.push_back( bodyPanelSettings( frameFixedPanelGeometry( Eigen::Vector3d::UnitX( ), 9.9 ),
                                                                specularDiffuseBodyPanelReflectionLawSettings( specular1, diffuse1, false ),
                                                                "Bus" ) );
                panelSettingsList.push_back( bodyPanelSettings( frameFixedPanelGeometry( Eigen::Vector3d::UnitY( ), 9.9 ),
                                                                specularDiffuseBodyPanelReflectionLawSettings( specular2, diffuse2, false ),
                                                                "Bus" ) );
                panelSettings = fullPanelledBodySettings( panelSettingsList );
            }
            else
            {
                panelSettings = bodyWingPanelledGeometry( 2., 3., 4., 0., 0.35, 0.25, 0.35, 0.25, false, false );
            }
            addBodyExteriorPanelledShape( panelSettings, "Vehicle", bodies );
            bodies.at( "Vehicle" )
                    ->setRadiationPressureTargetModels( { createRadiationPressureTargetModel(
                            std::make_shared< PaneledRadiationPressureTargetModelSettings >( ), "Vehicle", bodies ) } );
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set accelerations on Vehicle that are to be taken into account.
        SelectedAccelerationMap accelerationSettingsList;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

        accelerationsOfVehicle[ "Earth" ] = { sphericalHarmonicAcceleration( 2, 2 ), radiationPressureAcceleration( ) };
        accelerationsOfVehicle[ "Sun" ] = { pointMassGravityAcceleration( ), radiationPressureAcceleration( ) };

        accelerationSettingsList[ "Vehicle" ] = accelerationsOfVehicle;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToIntegrate = { "Vehicle" };
        std::vector< std::string > centralBodies = { "Earth" };

        // Create acceleration models
        AccelerationMap accelerationModelMap =
                createAccelerationModelsMap( bodies, accelerationSettingsList, bodiesToIntegrate, centralBodies );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set Keplerian elements for Vehicle.
        Eigen::Vector6d vehicleInitialStateInKeplerianElements;
        vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
        vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
        vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
        vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
        vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

        // Set initial state.
        double earthGravitationalParameter = getBodyGravitationalParameter( bodies, "Earth" );
        Eigen::Matrix< double, 6, 1 > systemInitialState =
                convertKeplerianToCartesianElements( vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

        vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7300.0E3;
        vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 95.3 );
        Eigen::Matrix< double, 6, 1 > systemInitialState2 =
                convertKeplerianToCartesianElements( vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

        // Create integrator settings
        std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKuttaFixedStepSettings( 40.0, rungeKuttaFehlberg78 );

        // Create propagator settings
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettingsArc1 =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies,
                        accelerationModelMap,
                        bodiesToIntegrate,
                        systemInitialState,
                        initialEphemerisTime,
                        integratorSettings,
                        propagationTimeTerminationSettings( finalEphemerisTime ) );

        // Create propagator settings
        double arcLength = 0.5 * 86400.0;
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettingsArc2 =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies,
                        accelerationModelMap,
                        bodiesToIntegrate,
                        systemInitialState2,
                        initialEphemerisTime + arcLength,
                        integratorSettings,
                        propagationTimeTerminationSettings( finalEphemerisTime + arcLength ) );

        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
        propagatorSettingsList.push_back( propagatorSettingsArc1 );
        propagatorSettingsList.push_back( propagatorSettingsArc2 );

        std::shared_ptr< MultiArcPropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList );


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////    DEFINE PARAMETERS THAT ARE TO BE ESTIMATED      ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define list of parameters to estimate.
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames = getInitialStateParameterSettings< double >( propagatorSettings, bodies );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", specular_reflectivity, "Bus" ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", diffuse_reflectivity, "Bus" ) );

        // Create parameters
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate( parameterNames, bodies );

        // Print identifiers and indices of parameters to terminal.
        printEstimatableParameterEntries( parametersToEstimate );

        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
        LinkEnds linkEnds;
        linkEnds[ observed_body ] = LinkEndId( "Vehicle", "" );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, linkEnds ) );

        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager =
                OrbitDeterminationManager< double, double >( bodies, parametersToEstimate, observationSettingsList, propagatorSettings );

        // Compute list of observation times.
        std::vector< double > baseTimeList;
        double observationTime = initialEphemerisTime + 1000.0;
        double observationInterval = 60.0;
        while( observationTime < finalEphemerisTime +  arcLength - 1000.0 )
        {
            baseTimeList.push_back( observationTime );
            observationTime += observationInterval;
        }

        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
        measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings<> >(
                position_observable, linkEnds, baseTimeList, observed_body ) );

        // Simulate observations.
        std::shared_ptr< ObservationCollection<> > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////    PERTURB PARAMETER VECTOR AND ESTIMATE PARAMETERS     ///////////        //
/////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Perturb parameter estimate.
        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
                parametersToEstimate->template getFullParameterValues< double >( );
        Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
        Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

        // Perturbe initial state estimate.
        parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
        parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );
        parameterPerturbation.segment( 6, 3 ) = Eigen::Vector3d::Constant( 1.0 );
        parameterPerturbation.segment( 9, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );
//        parameterPerturbation.segment( 12, 1 ) = Eigen::Vector1d::Constant( 0.1 );
//        parameterPerturbation.segment( 13, 1 ) = Eigen::Vector1d::Constant( 0.1 );

        initialParameterEstimate += parameterPerturbation;
        parametersToEstimate->resetParameterValues( initialParameterEstimate );
        Eigen::Matrix< double, Eigen::Dynamic, 1 > testValues = parametersToEstimate->template getFullParameterValues< double >( );
        std::cout << "TEST A:" << truthParameters.transpose( ) << std::endl;
        std::cout << "TEST B:" << testValues.transpose( ) << std::endl;

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( observationsAndTimes );

        estimationInput->defineEstimationSettings( true, true, true, true, true );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 4 ) );

        // Perform estimation
        std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

        double rmsResidual = linear_algebra::getVectorEntryRootMeanSquare( observationsAndTimes->getConcatenatedResiduals( ) );
        BOOST_CHECK_SMALL( rmsResidual, 1.0E-4 );
        Eigen::VectorXd parameterEstimate = estimationOutput->parameterEstimate_;
        std::cout << "parameter estimate: " << ( parameterEstimate ).transpose( ) << std::endl;

        Eigen::VectorXd estimationError = parameterEstimate - truthParameters;
//        for( int i = 0; i < 3; i++ )
//        {
//            BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 2.0E-4 );
//            BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 2.0E-7 );
//        }
//        BOOST_CHECK_SMALL( std::fabs( estimationError( 6 ) ), 2.0E-4 );
//        BOOST_CHECK_SMALL( std::fabs( estimationError( 6 ) ), 2.0E-4 );

        Eigen::Matrix< double, Eigen::Dynamic, 1 > newTestValues = parametersToEstimate->template getFullParameterValues< double >( );
        std::cout << "TEST C:" << newTestValues.transpose( ) << std::endl;

        std::cout << "estimation error: " << ( estimationError ).transpose( ) << std::endl;
        std::cout << "truth value: " << ( truthParameters ).transpose( ) << std::endl;
    }
}

//BOOST_AUTO_TEST_SUITE_END( )
//
//}  // namespace unit_tests
//
//}  // namespace tudat
