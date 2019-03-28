/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

//#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/simulateObservations.h"
#include "Tudat/SimulationSetup/EstimationSetup/orbitDeterminationManager.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EstimationSetup/podProcessing.h"


//namespace tudat
//{
//namespace unit_tests
//{
//BOOST_AUTO_TEST_SUITE( test_multi_arc_state_estimation )

//Using declarations.
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

template< typename ObservationScalarType = double , typename TimeType = double , typename StateScalarType  = double >
Eigen::VectorXd  executeMultiBodyMultiArcParameterEstimation(
        const int linkArcs )
{
    //Load spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Specify simulation time
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = initialEphemerisTime + 4.0 * 86400.0;

    //Define setting for celestial bodies
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - 86400.0, finalEphemerisTime + 86400.0 );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // CReate vehicles
    std::vector< std::string > vehicleNames = { "Borzi1" };//, "Borzi2" };
    int numberOfVehicles = vehicleNames.size( );
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        bodyMap[ vehicleNames.at( i ) ] = std::make_shared< Body >( );
        bodyMap[ vehicleNames.at( i ) ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                                           std::map< double, std::shared_ptr< Ephemeris > >( ), "Earth", "ECLIPJ2000" ) );
    }

    // Finalize environment creation
    setGlobalFrameBodyEphemerides( bodyMap, "Earth", "ECLIPJ2000" );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ vehicleNames.at( i ) ] = accelerationsOfVehicle;
    }

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate = vehicleNames;
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );
    std::vector< std::string > centralBodies;
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodies.push_back( "Earth" );
    }

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToIntegrate, centralBodies );


    // Define integration arc limits
    std::vector< double > integrationArcStartTimes;
    std::vector< double > integrationArcEndTimes;
    std::vector< double > integrationArcLimits;

    double integrationStartTime = initialEphemerisTime;
    double integrationEndTime = finalEphemerisTime;
    double arcDuration = 1.01 * 86400.0;
    double arcOverlap = 300.0;
    double currentStartTime = integrationStartTime;
    double currentEndTime = integrationStartTime + arcDuration;
    do
    {
        integrationArcLimits.push_back( currentStartTime );
        integrationArcEndTimes.push_back( currentEndTime );
        integrationArcStartTimes.push_back( currentStartTime );
        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < integrationEndTime );

    integrationArcLimits.push_back( currentStartTime + arcOverlap );


    // Create initial state container per arc
    std::vector< Eigen::VectorXd > allBodiesPerArcInitialStates;
    allBodiesPerArcInitialStates.resize( integrationArcStartTimes.size( ) );
    for( int j = 0; j < integrationArcStartTimes.size( ); j++ )
    {
        allBodiesPerArcInitialStates[ j ] = Eigen::VectorXd::Zero( 6 * numberOfVehicles );
    }

    // Create initial state container per body
    std::vector< Eigen::VectorXd > singleBodyForAllArcsInitialStates;
    singleBodyForAllArcsInitialStates.resize( numberOfVehicles );
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        singleBodyForAllArcsInitialStates[ i ] = Eigen::VectorXd::Zero( 6 * integrationArcStartTimes.size( ) );
    }

    // Define (semi-arbitrary) initial states for vehicles
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        Eigen::Vector6d nominalInitialStateInKeplerianElements;
        nominalInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
        nominalInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        nominalInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        nominalInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
        nominalInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
        nominalInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

        for( int j = 0; j < integrationArcStartTimes.size( ); j++ )
        {
            Eigen::Vector6d currentInitialStateInKeplerianElements = nominalInitialStateInKeplerianElements;
            currentInitialStateInKeplerianElements( inclinationIndex ) = 10.0 * static_cast< double >( i );
            currentInitialStateInKeplerianElements( trueAnomalyIndex ) = 10.0 * static_cast< double >( j );
            double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
            const Eigen::Vector6d currentInitialState = convertKeplerianToCartesianElements(
                        currentInitialStateInKeplerianElements, earthGravitationalParameter );
            allBodiesPerArcInitialStates[ j ].segment( i * 6, 6 ) = currentInitialState;
            singleBodyForAllArcsInitialStates[ i ].segment( j * 6, 6 ) = currentInitialState;
        }
    }


    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, TimeType( initialEphemerisTime  ), 30.0 );

    // Define propagator settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > propagatorSettingsList;
    for( unsigned int i = 0; i < integrationArcStartTimes.size( ); i++ )
    {
        propagatorSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                    ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                      allBodiesPerArcInitialStates.at( i ),
                      integrationArcEndTimes.at( i ), cowell, std::shared_ptr< DependentVariableSaveSettings >( ), 60.0 ) );
    }
    std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< StateScalarType > >( propagatorSettingsList );


    // Set parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                                      vehicleNames.at( i ), singleBodyForAllArcsInitialStates.at( i ), integrationArcStartTimes, "Earth" ) );
    }
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodyMap );


    // Define links and observations in simulation.
    std::vector< LinkEnds > linkEndsList;
    observation_models::ObservationSettingsMap observationSettingsMap;
    linkEndsList.resize( numberOfVehicles );
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        linkEndsList[ i ][ observed_body ] = std::make_pair( vehicleNames.at( i ), "" );
        observationSettingsMap.insert( std::make_pair( linkEndsList[ i ], std::make_shared< ObservationSettings >(
                                                           position_observable ) ) );
    }


    // Create orbit determination object.
    OrbitDeterminationManager< ObservationScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< ObservationScalarType, TimeType >(
                bodyMap, parametersToEstimate,
                observationSettingsMap, integratorSettings, propagatorSettings );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );


    // Define times at which observations are to be simulated
    TimeType observationTime;
    int numberOfObservationsPerArc = 100;
    double timeBuffer = 300.0;

    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservationsPerArc * integrationArcStartTimes.size( ) );
    for( unsigned int i = 0; i < integrationArcLimits.size( ) - 1; i++ )
    {
        double currentTimeStep = ( integrationArcLimits[ i + 1 ] - integrationArcLimits[ i ] - 2.0 * timeBuffer ) /
                static_cast< double >( numberOfObservationsPerArc - 1 );
        observationTime = integrationArcLimits[ i ] + timeBuffer;
        for( int j = 0; j < numberOfObservationsPerArc; j++ )
        {
            initialObservationTimes[ j + i * numberOfObservationsPerArc ] = observationTime;
            observationTime += currentTimeStep;
        }
    }


    // Define observation simulation settings
    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > measurementSimulationInput;
    std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > singleObservableSimulationInput;
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        singleObservableSimulationInput[ linkEndsList[ i ] ] = std::make_pair( initialObservationTimes, observed_body );
    }
    measurementSimulationInput[ position_observable ] = singleObservableSimulationInput;

    // Simulate observations
    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;
    PodInputDataType observationsAndTimes = simulateObservations< ObservationScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( )  );

    // Perturb initial states
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    for( unsigned int i = 0; i < numberOfNumericalBodies * integrationArcStartTimes.size( ); i++ )
    {
        initialParameterEstimate[ 0 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 1 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 2 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 3 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 4 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 5 + 6 * i ] += 1.0E-5;
    }
    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    // Define POD input
    std::shared_ptr< PodInput< ObservationScalarType, TimeType > > podInput =
            std::make_shared< PodInput< ObservationScalarType, TimeType > >(
                observationsAndTimes, ( initialParameterEstimate ).rows( ) );

    // Estimate parameters
    std::shared_ptr< PodOutput< StateScalarType, TimeType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( 3 ) );

//    std::cout<<"Initial states difference : "<<
//               ( propagatorSettings->getInitialStates( ) - originalInitiaStates ).transpose( )<<std::endl<<std::endl;
//    std::cout<<"Initial estimated states difference : "<<
//               ( parametersToEstimate->template getFullParameterValues< double >( ) - originalParameterEstimate ).transpose( )<<std::endl<<std::endl;

    std::string outputFolder = "/home/dominic/Software/tudatBundleTest/tudatBundle/tudatApplications/master_thesis/Output/";

    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "earthOrbitEstimationInformationMatrix.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
                                     "earthOrbitEstimationInformationMatrixNormalization.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( podOutput->weightsMatrixDiagonal_,
                                     "earthOrbitEstimationWeightsDiagonal.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( podOutput->residuals_,
                                     "earthOrbitEstimationResiduals.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                     "earthOrbitEstimationCorrelations.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ),
                                     "earthOrbitResidualHistory.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ),
                                     "earthOrbitParameterHistory.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "earthOrbitObservationMeasurements.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
                                     "earthOrbitObservationTimes.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
                                     "earthOrbitObservationLinkEnds.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedObservableTypes( podInput->getObservationsAndTimes( ) ) ),
                                     "earthOrbitObservationObservableTypes.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "earthOrbitObservationMeasurements.dat", 16,
                                     outputFolder );
    //    input_output::writeMatrixToFile( estimationError,
    //                                     "earthOrbitObservationTrueEstimationError.dat", 16,
    //                                     outputFolder );
    input_output::writeMatrixToFile( podOutput->getFormalErrorVector( ),
                                     "earthOrbitObservationFormalEstimationError.dat", 16,
                                     outputFolder );

    return ( podOutput->parameterEstimate_ - truthParameters ).template cast< double >( );
}

//BOOST_AUTO_TEST_CASE( test_MultiArcStateEstimation )
int main( )
{
    // Execute test for linked arcs and separate arcs.
    //    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
        Eigen::VectorXd parameterError = executeMultiBodyMultiArcParameterEstimation< double, double, double >(
                    false );
        int numberOfEstimatedArcs = ( parameterError.rows( ) - 3 ) / 6;

        std::cout <<"Estimation error: "<< parameterError.transpose( ) << std::endl;
        //        for( int i = 0; i < numberOfEstimatedArcs; i++ )
        //        {
        //            for( unsigned int j = 0; j < 3; j++ )
        //            {
        //                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j ) ), 1E-4 );
        //                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j + 3 ) ), 1.0E-10  );
        //            }
        //        }
    }

}

//BOOST_AUTO_TEST_SUITE_END( )

//}

//}
