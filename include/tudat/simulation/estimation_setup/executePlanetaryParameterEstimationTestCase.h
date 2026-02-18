/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef EXECUTEPLANETARYPARAMETERESTIMATIONTESTCASE_H
#define EXECUTEPLANETARYPARAMETERESTIMATIONTESTCASE_H

#include <memory>
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include <utility>

#include <Eigen/Core>

#include "tudat/basics/timeType.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/estimationInterfacesForwardDeclarations.h"
#include "tudat/simulation/estimation_setup/compareEstimationAndCovarianceResultsTestCase.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationTestCaseUtilities.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"


namespace tudat
{
namespace unit_tests
{

using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace estimatable_parameters;
using namespace propagators;
using namespace numerical_integrators;

Eigen::VectorXd getDefaultInitialParameterPerturbation( );

template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::shared_ptr< simulation_setup::EstimationOutput< StateScalarType, TimeType > >, Eigen::VectorXd >
executePlanetaryParameterEstimation(
        const int observableType = 1,
        Eigen::VectorXd parameterPerturbation = getDefaultInitialParameterPerturbation( ),
        Eigen::MatrixXd inverseAPrioriCovariance = Eigen::MatrixXd::Zero( 7, 7 ),
        const double weight = 1.0 );

extern template std::pair< std::shared_ptr< simulation_setup::EstimationOutput< double, double > >, Eigen::VectorXd >
executePlanetaryParameterEstimation< double, double >( const int observableType,
                                                       Eigen::VectorXd parameterPerturbation,
                                                       Eigen::MatrixXd inverseAPrioriCovariance,
                                                       const double weight );

template< typename TimeType, typename StateScalarType >
std::pair< std::shared_ptr< EstimationOutput< StateScalarType, TimeType > >, Eigen::VectorXd >
executePlanetaryParameterEstimation(
        const int observableType,
        Eigen::VectorXd parameterPerturbation,
        Eigen::MatrixXd inverseAPrioriCovariance,
        const double weight )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define setting for total number of bodies and those which need to be integrated numerically.
    // The first numberOfNumericalBodies from the bodyNames vector will be integrated numerically.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Saturn" );

    // Specify initial time.
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = TimeType( 3.0E7 );
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( "Sun" );

    // Create bodies needed in simulation.
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Earth" ] = accelerationsOfEarth;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Earth" );
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    // Define propagator settings.
    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralBodyMap;
    centralBodies.resize( numberOfNumericalBodies );
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodies[ i ] = "SSB";
        centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
    }

    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, centralBodyMap );

    // Set parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
            "Earth",
            propagators::getInitialStateOfBody< TimeType, StateScalarType >(
                    "Earth", centralBodyMap.at( "Earth" ), bodies, initialEphemerisTime ),
            centralBodyMap.at( "Earth" ) ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Moon", gravitational_parameter ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType, TimeType >( parameterNames, bodies );

    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings = std::make_shared< IntegratorSettings< TimeType > >(
            rungeKutta4, TimeType( initialEphemerisTime - 4.0 * maximumTimeStep ), 900.0 );

    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                    centralBodies,
                    accelerationModelMap,
                    bodiesToIntegrate,
                    getInitialStateVectorOfBodiesToEstimate( parametersToEstimate ),
                    TimeType( finalEphemerisTime + 4.0 * maximumTimeStep ),
                    cowell );

    // Define link ends.
    LinkEnds linkEnds;
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    if( observableType == 0 )
    {
        linkEnds[ observed_body ] = LinkEndId( "Earth", "" );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, linkEnds ) );
    }
    else if( observableType == 5 )
    {
        linkEnds[ observed_body ] = LinkEndId( "Earth", "" );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( velocity_observable, linkEnds ) );
    }
    else
    {
        linkEnds[ transmitter ] = LinkEndId( "Earth", "" );
        linkEnds[ receiver ] = LinkEndId( "Mars", "" );

        if( observableType == 1 )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, linkEnds ) );
        }
        else if( observableType == 2 )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( angular_position, linkEnds ) );
        }
        else if( observableType == 3 )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_doppler, linkEnds ) );
        }
        else if( observableType == 4 )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, linkEnds ) );
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_doppler, linkEnds ) );
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( angular_position, linkEnds ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                    bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    // Define observation times.
    double observationTimeStep = 1000.0;
    TimeType observationTime = Time( initialEphemerisTime + 10.0E4 );
    int numberOfObservations = 18000;
    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservations );
    for( int i = 0; i < numberOfObservations; i++ )
    {
        initialObservationTimes[ i ] = observationTime;
        observationTime += observationTimeStep;
    }

    // Define observation simulation settings.
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;
    initialObservationTimes = utilities::addScalarToVector( initialObservationTimes, 30.0 );
    if( observableType == 0 )
    {
        measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                position_observable, linkEnds, initialObservationTimes, observed_body ) );
    }
    else if( observableType == 5 )
    {
        measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                velocity_observable, linkEnds, initialObservationTimes, observed_body ) );
    }
    else
    {
        if( observableType == 1 )
        {
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                    one_way_range, linkEnds, initialObservationTimes, transmitter ) );
        }
        else if( observableType == 2 )
        {
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                    angular_position, linkEnds, initialObservationTimes, transmitter ) );
        }
        else if( observableType == 3 )
        {
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                    one_way_doppler, linkEnds, initialObservationTimes, transmitter ) );
        }
        else if( observableType == 4 )
        {
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                    one_way_range, linkEnds, initialObservationTimes, transmitter ) );
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                    angular_position, linkEnds, initialObservationTimes, transmitter ) );
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                    one_way_doppler, linkEnds, initialObservationTimes, transmitter ) );
        }
    }

    // Simulate observations.
    std::shared_ptr< ObservationCollection< StateScalarType, TimeType > > simulatedObservations =
            simulateObservations< StateScalarType, TimeType >(
                    measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    if( observableType == 4 )
    {
        std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, double > weightsPerObservationParser;
        weightsPerObservationParser[ observationParser( one_way_range ) ] = 1.0 / ( 1.0 * 1.0 );
        weightsPerObservationParser[ observationParser( angular_position ) ] = 1.0 / ( 1.0E-9 * 1.0E-9 );
        weightsPerObservationParser[ observationParser( one_way_doppler ) ] = 1.0 / ( 1.0E-12 * 1.0E-12 );
        simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );
    }
    else
    {
        simulatedObservations->setConstantWeight( weight );
    }

    // Perturb parameter estimate.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    if( parameterPerturbation.rows( ) == 0 )
    {
        parameterPerturbation = Eigen::VectorXd::Zero( 7 );
    }
    for( unsigned int i = 0; i < initialParameterEstimate.rows( ); i++ )
    {
        initialParameterEstimate( i ) += parameterPerturbation( i );
    }
    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    // Define estimation input.
    std::shared_ptr< EstimationInput< StateScalarType, TimeType > > estimationInput =
            std::make_shared< EstimationInput< StateScalarType, TimeType > >( simulatedObservations, inverseAPrioriCovariance );
    std::shared_ptr< CovarianceAnalysisInput< StateScalarType, TimeType > > covarianceInput =
            std::make_shared< EstimationInput< StateScalarType, TimeType > >( simulatedObservations, inverseAPrioriCovariance );
    estimationInput->defineEstimationSettings( true, true, false, true, true );
    covarianceInput->defineCovarianceSettings( true, true, true, false );
    estimationInput->applyFinalParameterCorrection_ = false;

    // Perform estimation.
    std::shared_ptr< EstimationOutput< StateScalarType, TimeType > > estimationOutput =
            orbitDeterminationManager.estimateParameters( estimationInput );

    parametersToEstimate->template resetParameterValues< StateScalarType >(
            estimationOutput->parameterHistory_.at( estimationOutput->bestIteration_ ) );
    std::shared_ptr< CovarianceAnalysisOutput< StateScalarType, TimeType > > covarianceOutput =
            orbitDeterminationManager.computeCovariance( covarianceInput );

    compareEstimationAndCovarianceResults( estimationOutput, covarianceOutput );

    return std::make_pair(
            estimationOutput,
            ( estimationOutput->parameterEstimate_.template cast< double >( ) - truthParameters.template cast< double >( ) ) );
}

}  // namespace unit_tests

}  // namespace tudat

#endif  // EXECUTEPLANETARYPARAMETERESTIMATIONTESTCASE_H
