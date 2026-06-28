/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SIMULATEOBSERVATIONS_H
#define TUDAT_SIMULATEOBSERVATIONS_H

#include <memory>

#include <functional>

#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/math/statistics/randomVariableGenerator.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"
#include "tudat/simulation/estimation_setup/createObservationViability.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tudat
{

namespace simulation_setup
{

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
void addNoiseAndDependentVariableToObservation(
        Eigen::Matrix< ObservationScalarType, ObservationSize, 1 >& calculatedObservation,
        const TimeType& observationTime,
        Eigen::VectorXd& dependentVariables,
        const std::vector< Eigen::Vector6d >& vectorOfStates,
        const std::vector< double >& vectorOfTimes,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings,
        const observation_models::ObservableType observableType,
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr,
        const std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr )
{
    if( dependentVariableCalculator != nullptr )
    {
        dependentVariables = dependentVariableCalculator->calculateDependentVariables(
                vectorOfTimes, vectorOfStates, calculatedObservation.template cast< double >( ), ancillarySettings );
    }

    // Add noise if needed.
    if( noiseFunction != nullptr )
    {
        Eigen::VectorXd noiseToAdd = noiseFunction( observationTime );
        if( noiseToAdd.rows( ) != ObservationSize )
        {
            throw std::runtime_error(
                    "Error wen simulating observation noise, size of noise (" + std::to_string( noiseToAdd.rows( ) ) +
                    ") and size of observable (" + std::to_string( ObservationSize ) +
                    ") are not compatible for observable type: " + observation_models::getObservableName( observableType ) );
        }
        else
        {
            calculatedObservation += noiseToAdd.template cast< ObservationScalarType >( );
        }
    }
}

//! Function to simulate an observable, checking whether it is viable according to settings passed to this function
/*!
 *  Function to simulate an observable, checking whether it is viable according to settings passed to this function
 *  (if viability calculators are passed to this function).
 *  \param observationTime Time at which observable is to be computed
 *  \param observationModel Model used to compute observable
 *  \param referenceLinkEnd Model Reference link end for observable
 *  \param linkViabilityCalculators List of observation viability calculators, which are used to reject simulated
 *  observation if they dont fulfill a given (set of) conditions, e.g. minimum elevation angle (default none).
 *  \return Observation at given time.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::tuple< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, bool, Eigen::VectorXd > simulateObservationWithCheck(
        const TimeType& observationTime,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > linkViabilityCalculators =
                std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > >( ),
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr,
        const std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
{
    // Simulate observable, and retrieve link end times and states
    std::vector< Eigen::Vector6d > vectorOfStates;
    std::vector< double > vectorOfTimes;
    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > calculatedObservation = observationModel->computeObservationsWithLinkEndData(
            observationTime, referenceLinkEnd, vectorOfTimes, vectorOfStates, ancillarySettings );
    Eigen::VectorXd dependentVariables = Eigen::VectorXd::Zero( 0 );

    // Check if observation is feasible
    bool observationFeasible = isObservationViable(
            vectorOfStates, vectorOfTimes, linkViabilityCalculators, calculatedObservation.template cast< double >( ) );

    if( observationFeasible )
    {
        addNoiseAndDependentVariableToObservation< ObservationSize, ObservationScalarType, TimeType >(
                calculatedObservation,
                observationTime,
                dependentVariables,
                vectorOfStates,
                vectorOfTimes,
                ancillarySettings,
                observationModel->getObservableType( ),
                noiseFunction,
                dependentVariableCalculator );
    }

    // Return simulated observable and viability
    return std::make_tuple( calculatedObservation, observationFeasible, dependentVariables );
}

//! Function to simulate observables, checking whether they are viable according to settings passed to this function
/*!
 *  Function to simulate observables, checking whether they are viable according to settings passed to this function
 *  (if viability calculators are passed to this function).
 *  \param observationTimes Times at which observables are to be computed
 *  \param observationModel Model used to compute observables
 *  \param referenceLinkEnd Model Reference link end for observables
 *  \param linkViabilityCalculators List of observation viability calculators, which are used to reject simulated
 *  observation if they dont fulfill a given (set of) conditions, e.g. minimum elevation angle (default none).
 *  \return Observations at given time (concatenated in an Eigen vector) and associated times.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::tuple< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >,
            std::vector< TimeType >,
            std::vector< Eigen::VectorXd > >
simulateObservationsWithCheck(
        const std::vector< TimeType >& observationTimes,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > linkViabilityCalculators =
                std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > >( ),
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr,
        const std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
{
    std::multimap< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
    std::tuple< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, bool, Eigen::VectorXd > simulatedObservation;
    std::vector< Eigen::VectorXd > dependentVariables;

    for( unsigned int i = 0; i < observationTimes.size( ); i++ )
    {
        simulatedObservation =
                simulateObservationWithCheck< ObservationSize, ObservationScalarType, TimeType >( observationTimes.at( i ),
                                                                                                  observationModel,
                                                                                                  referenceLinkEnd,
                                                                                                  linkViabilityCalculators,
                                                                                                  noiseFunction,
                                                                                                  dependentVariableCalculator,
                                                                                                  ancillarySettings );

        // Check if receiving station can view transmitting station.
        if( std::get< 1 >( simulatedObservation ) )
        {
            // If viable, add observable and time to vector of simulated data.
            observations.insert( { observationTimes[ i ], std::get< 0 >( simulatedObservation ) } );
            dependentVariables.push_back( std::get< 2 >( simulatedObservation ) );
        }
    }

    // Return pair of simulated ranges and reception times.
    return std::make_tuple( utilities::createVectorFromMultiMapValues( observations ),
                            utilities::createVectorFromMultiMapKeys( observations ),
                            dependentVariables );
}

//! Function to simulate observables, checking whether they are viable according to settings passed to this function
/*!
 *  Function to simulate observables, checking whether they are viable according to settings passed to this function
 *  (if viability calculators are passed to this function).
 *  \param observationTimes Times at which observables are to be computed
 *  \param observationModel Model used to compute observables
 *  \param referenceLinkEnd Model Reference link end for observables
 *  \param linkViabilityCalculators List of observation viability calculators, which are used to reject simulated
 *  observation if they dont fulfill a given (set of) conditions, e.g. minimum elevation angle (default none).
 *  \return Observations at given times (concatenated in an Eigen vector), with associated times and reference link end.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > >
simulateObservationDatasetWithCheckAndLinkEndIdOutput(
        const std::vector< TimeType >& observationTimes,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > linkViabilityCalculators =
                std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > >( ),
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr,
        const std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
{
    std::tuple< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >,
                std::vector< TimeType >,
                std::vector< Eigen::VectorXd > >
            simulatedObservations =
                    simulateObservationsWithCheck< ObservationSize, ObservationScalarType, TimeType >( observationTimes,
                                                                                                       observationModel,
                                                                                                       referenceLinkEnd,
                                                                                                       linkViabilityCalculators,
                                                                                                       noiseFunction,
                                                                                                       dependentVariableCalculator,
                                                                                                       ancillarySettings );
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset =
            std::make_shared< observation_models::ObservationDataset< ObservationScalarType, TimeType > >( );
    observationDataset->addObservationSet(
            observationModel->getObservableType( ),
            observationModel->getLinkEnds( ),
            std::get< 0 >( simulatedObservations ),
            std::get< 1 >( simulatedObservations ),
            referenceLinkEnd,
            std::get< 2 >( simulatedObservations ),
            ( dependentVariableCalculator == nullptr ) ? nullptr : dependentVariableCalculator->getDependentVariableBookkeeping( ),
            ancillarySettings,
            std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
            true );
    return observationDataset;
}

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > >
simulateObservationsWithCheckAndLinkEndIdOutput(
        const std::vector< TimeType >& observationTimes,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > linkViabilityCalculators =
                std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > >( ),
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr,
        const std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
{
    return observation_models::createSingleObservationSet< ObservationScalarType, TimeType >(
            simulateObservationDatasetWithCheckAndLinkEndIdOutput< ObservationSize, ObservationScalarType, TimeType >(
                    observationTimes,
                    observationModel,
                    referenceLinkEnd,
                    linkViabilityCalculators,
                    noiseFunction,
                    dependentVariableCalculator,
                    ancillarySettings ) );
}

template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > simulatePerArcObservationDataset(
        const std::shared_ptr< PerArcObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const SystemOfBodies& bodies )
{
    using namespace observation_models;

    // Create viability settings for arc-defining constraint
    std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > arcDefiningViabilityCalculators =
            observation_models::createObservationViabilityCalculators( bodies,
                                                                       observationsToSimulate->getLinkEnds( ).linkEnds_,
                                                                       observationsToSimulate->getObservableType( ),
                                                                       { observationsToSimulate->arcDefiningConstraint_ } );

    std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator =
            std::make_shared< ObservationDependentVariableCalculator >(
                    observationsToSimulate->getObservationDependentVariableBookkeeping( ),
                    bodies,
                    observationModel != nullptr
                            ? observationModel->getLegLightTimeCalculators( )
                            : std::map< std::pair< observation_models::LinkEndType, observation_models::LinkEndType >,
                                        std::vector< std::shared_ptr< observation_models::LightTimeCalculatorBase > > >( ) );

    // Define list of arc data
    typedef std::tuple< Eigen::Matrix< ObservationScalarType, ObservationSize, 1 >, std::vector< Eigen::Vector6d >, std::vector< double > >
            SingleObservationData;
    typedef std::map< TimeType, SingleObservationData > SingleArcObservationData;
    std::vector< SingleArcObservationData > simulatedObservations;

    // Declare per-observable variables
    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation;
    std::vector< Eigen::Vector6d > vectorOfStates;
    std::vector< double > vectorOfTimes;
    bool observationFeasible;

    // Declare list of arc observations
    SingleArcObservationData currentObservationArc;

    // Initialize observation simulation
    LinkEndType referenceLinkEnd = observationsToSimulate->getReferenceLinkEndType( );
    TimeType currentObservationTime = observationsToSimulate->startTime_;
    std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings =
            observationsToSimulate->getAncillarySettings( );

    while( currentObservationTime < observationsToSimulate->endTime_ )
    {
        bool addTimeInterval = true;

        // Simulate observation
        currentObservation = observationModel->computeObservationsWithLinkEndData(
                currentObservationTime, referenceLinkEnd, vectorOfTimes, vectorOfStates, ancillarySettings );

        // If observation is feasible, add to arc. If not, check if current arc is to be terminated.
        observationFeasible = isObservationViable(
                vectorOfStates, vectorOfTimes, arcDefiningViabilityCalculators, currentObservation.template cast< double >( ) );
        if( observationFeasible )
        {
            bool isMaximumArcDurationExceeded = false;
            if( currentObservationArc.size( ) > 0 )
            {
                TimeType arcInitialTime = currentObservationArc.begin( )->first;
                if( ( ( currentObservationTime - arcInitialTime ) > observationsToSimulate->maximumArcDuration_ ) &&
                    ( observationsToSimulate->maximumArcDuration_ == observationsToSimulate->maximumArcDuration_ ) )
                {
                    isMaximumArcDurationExceeded = true;
                }
            }

            if( !isMaximumArcDurationExceeded )
            {
                currentObservationArc[ currentObservationTime ] = std::make_tuple( currentObservation, vectorOfStates, vectorOfTimes );
                if( observationsToSimulate->minimumTimeBetweenArcs_ == observationsToSimulate->minimumTimeBetweenArcs_ )
                {
                    currentObservationTime += observationsToSimulate->minimumTimeBetweenArcs_;
                    addTimeInterval = false;
                }
            }
        }
        else if( currentObservationArc.size( ) > 0 )
        {
            TimeType arcInitialTime = currentObservationArc.begin( )->first;
            TimeType arcFinalTime = currentObservationArc.rbegin( )->first;

            if( ( arcFinalTime - arcInitialTime ) > observationsToSimulate->minimumArcDuration_ ||
                !( observationsToSimulate->minimumArcDuration_ == observationsToSimulate->minimumArcDuration_ ) )
            {
                simulatedObservations.push_back( currentObservationArc );
            }

            currentObservationArc.clear( );
        }

        if( addTimeInterval )
        {
            currentObservationTime += observationsToSimulate->intervalBetweenObservations_;
        }
    }

    std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > additionalViabilityCalculators =
            observation_models::createObservationViabilityCalculators( bodies,
                                                                       observationsToSimulate->getLinkEnds( ).linkEnds_,
                                                                       observationsToSimulate->getObservableType( ),
                                                                       observationsToSimulate->additionalViabilitySettingsList_ );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
    std::vector< TimeType > observationTimes;
    std::vector< Eigen::VectorXd > observationsDependentVariables;

    Eigen::VectorXd currentDependentVariable;
    for( unsigned int i = 0; i < simulatedObservations.size( ); i++ )
    {
        for( auto it : simulatedObservations.at( i ) )
        {
            SingleObservationData singleObservation = it.second;
            currentObservation = std::get< 0 >( singleObservation );
            vectorOfStates = std::get< 1 >( singleObservation );
            vectorOfTimes = std::get< 2 >( singleObservation );
            currentDependentVariable = Eigen::VectorXd::Zero( 0 );

            observationFeasible = isObservationViable(
                    vectorOfStates, vectorOfTimes, additionalViabilityCalculators, currentObservation.template cast< double >( ) );
            if( observationFeasible )
            {
                addNoiseAndDependentVariableToObservation< ObservationSize, ObservationScalarType, TimeType >(
                        currentObservation,
                        it.first,
                        currentDependentVariable,
                        vectorOfStates,
                        vectorOfTimes,
                        ancillarySettings,
                        observationModel->getObservableType( ),
                        observationsToSimulate->getObservationNoiseFunction( ),
                        dependentVariableCalculator );
                observations.push_back( currentObservation );
                observationTimes.push_back( it.first );
                observationsDependentVariables.push_back( currentDependentVariable );
            }
        }
    }

    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset =
            std::make_shared< observation_models::ObservationDataset< ObservationScalarType, TimeType > >( );
    observationDataset->addObservationSet( observationModel->getObservableType( ),
                                           observationModel->getLinkEnds( ),
                                           observations,
                                           observationTimes,
                                           referenceLinkEnd,
                                           observationsDependentVariables,
                                           observationsToSimulate->getObservationDependentVariableBookkeeping( ),
                                           ancillarySettings,
                                           std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                           std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                                           true );
    return observationDataset;
}

template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > simulatePerArcSingleObservationSet(
        const std::shared_ptr< PerArcObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const SystemOfBodies& bodies )
{
    return observation_models::createSingleObservationSet< ObservationScalarType, TimeType >(
            simulatePerArcObservationDataset< ObservationScalarType, TimeType, ObservationSize >(
                    observationsToSimulate, observationModel, bodies ) );
}

//! Function to compute observations at times defined by settings object using a given observation model
/*!
 *  Function to compute observations at times defined by settings object using a given observation model
 *  \param observationsToSimulate Object that computes/defines settings for observation times/reference link end
 *  \param observationModel Observation model that is to be used to compute observations
 *  \param currentObservationViabilityCalculators List of observation viability calculators, which are used to reject simulated
 *  observation if they dont fulfill a given (set of) conditions, e.g. minimum elevation angle (default none).
 *  \return Pair of observable values and observation time (with associated reference link end)
 */
template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > simulateSingleObservationDataset(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const SystemOfBodies& bodies )
{
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > simulatedObservations;
    //! Function to create an list of obervation viability conditions for a single set of link ends

    std::function< Eigen::VectorXd( const double ) > noiseFunction = observationsToSimulate->getObservationNoiseFunction( );

    // Simulate observations from tabulated times.
    if( std::dynamic_pointer_cast< TabulatedObservationSimulationSettings< TimeType > >( observationsToSimulate ) != nullptr )
    {
        std::shared_ptr< TabulatedObservationSimulationSettings< TimeType > > tabulatedObservationSettings =
                std::dynamic_pointer_cast< TabulatedObservationSimulationSettings< TimeType > >( observationsToSimulate );

        std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > currentObservationViabilityCalculators =
                observation_models::createObservationViabilityCalculators( bodies,
                                                                           observationsToSimulate->getLinkEnds( ).linkEnds_,
                                                                           observationsToSimulate->getObservableType( ),
                                                                           observationsToSimulate->getViabilitySettingsList( ) );

        std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator =
                std::make_shared< ObservationDependentVariableCalculator >(
                        tabulatedObservationSettings->getObservationDependentVariableBookkeeping( ),
                        bodies,
                        observationModel != nullptr
                                ? observationModel->getLegLightTimeCalculators( )
                                : std::map< std::pair< observation_models::LinkEndType, observation_models::LinkEndType >,
                                            std::vector< std::shared_ptr< observation_models::LightTimeCalculatorBase > > >( ) );

        // Simulate observations at requested pre-defined time.
        simulatedObservations = simulateObservationDatasetWithCheckAndLinkEndIdOutput< ObservationSize, ObservationScalarType, TimeType >(
                tabulatedObservationSettings->simulationTimes_,
                observationModel,
                observationsToSimulate->getReferenceLinkEndType( ),
                currentObservationViabilityCalculators,
                noiseFunction,
                dependentVariableCalculator,
                tabulatedObservationSettings->getAncillarySettings( ) );
    }
    else if( std::dynamic_pointer_cast< PerArcObservationSimulationSettings< TimeType > >( observationsToSimulate ) != nullptr )
    {
        std::shared_ptr< PerArcObservationSimulationSettings< TimeType > > perArcObservationSettings =
                std::dynamic_pointer_cast< PerArcObservationSimulationSettings< TimeType > >( observationsToSimulate );

        simulatedObservations = simulatePerArcObservationDataset< ObservationScalarType, TimeType, ObservationSize >(
                perArcObservationSettings, observationModel, bodies );
    }

    return simulatedObservations;
}

template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > simulateSingleObservationSet(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const SystemOfBodies& bodies )
{
    return observation_models::createSingleObservationSet< ObservationScalarType, TimeType >(
            simulateSingleObservationDataset< ObservationScalarType, TimeType, ObservationSize >(
                    observationsToSimulate, observationModel, bodies ) );
}

//! Function to simulate observations for single observable and single set of link ends.
/*!
 *  Function to simulate observations for single observable and single set of link ends. From the observation time settings and
 *  the observation simulator, the required observations are simulated and returned.
 *  \param observationsToSimulate Object that computes/defines settings for observation times/reference link end
 *  \param observationSimulator Observation simulator for observable for which observations are to be calculated.
 *  \param linkEnds Link end set for which observations are to be calculated.
 *  \return Pair of first: vector of observations; second: vector of times at which observations are taken
 *  (reference to link end defined in observationsToSimulate).
 */
template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > simulateSingleObservationDataset(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > >
                observationSimulator,
        const SystemOfBodies& bodies )
{
    if( observationSimulator == nullptr )
    {
        throw std::runtime_error( "Error when simulating single observation set, Observation simulator is nullptr" );
    }

    return simulateSingleObservationDataset< ObservationScalarType, TimeType, ObservationSize >(
            observationsToSimulate, observationSimulator->getObservationModel( observationsToSimulate->getLinkEnds( ).linkEnds_ ), bodies );
}

template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > simulateSingleObservationSet(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > >
                observationSimulator,
        const SystemOfBodies& bodies )
{
    return observation_models::createSingleObservationSet< ObservationScalarType, TimeType >(
            simulateSingleObservationDataset< ObservationScalarType, TimeType, ObservationSize >(
                    observationsToSimulate, observationSimulator, bodies ) );
}

//! Function to simulate observations from set of observables and link and sets
/*!
 *  Function to simulate observations from set of observables, link ends and observation time settings
 *  Iterates over all observables and link ends and simulates observations.
 *  \param observationsToSimulate List of observation time settings per link end set per observable type.
 *  \param observationSimulators List of Observation simulators per link end set per observable type.
 *  \param viabilityCalculatorList List (per observable type and per link ends) of observation viability calculators, which
 *  are used to reject simulated observation if they dont fulfill a given (set of) conditions, e.g. minimum elevation angle
 *  (default none).
 *  \return Simulated observatoon values and associated times for requested observable types and link end sets.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > simulateObservationDataset(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationsToSimulate,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulators,
        const SystemOfBodies bodies )
{
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset =
            std::make_shared< observation_models::ObservationDataset< ObservationScalarType, TimeType > >( );

    // Iterate over all observables.
    for( unsigned int i = 0; i < observationsToSimulate.size( ); i++ )
    {
        observation_models::ObservableType observableType = observationsToSimulate.at( i )->getObservableType( );
        observation_models::LinkEnds linkEnds = observationsToSimulate.at( i )->getLinkEnds( ).linkEnds_;

        int observationSize = observation_models::getObservableSize( observableType );

        switch( observationSize )
        {
            case 1: {
                std::shared_ptr< observation_models::ObservationSimulator< 1, ObservationScalarType, TimeType > >
                        derivedObservationSimulator =
                                observation_models::getObservationSimulatorOfType< 1 >( observationSimulators, observableType );

                if( derivedObservationSimulator == nullptr )
                {
                    throw std::runtime_error( "Error when simulating observation: dynamic cast to size 1 is nullptr" );
                }

                // Simulate observations for current observable and link ends set.
                observationDataset->addObservationSetFromDataset(
                        *simulateSingleObservationDataset< ObservationScalarType, TimeType, 1 >(
                                observationsToSimulate.at( i ), derivedObservationSimulator, bodies ),
                        0 );
                break;
            }
            case 2: {
                std::shared_ptr< observation_models::ObservationSimulator< 2, ObservationScalarType, TimeType > >
                        derivedObservationSimulator =
                                observation_models::getObservationSimulatorOfType< 2 >( observationSimulators, observableType );

                if( derivedObservationSimulator == nullptr )
                {
                    throw std::runtime_error( "Error when simulating observation: dynamic cast to size 2 is nullptr" );
                }

                // Simulate observations for current observable and link ends set.
                observationDataset->addObservationSetFromDataset(
                        *simulateSingleObservationDataset< ObservationScalarType, TimeType, 2 >(
                                observationsToSimulate.at( i ), derivedObservationSimulator, bodies ),
                        0 );
                break;
            }
            case 3: {
                std::shared_ptr< observation_models::ObservationSimulator< 3, ObservationScalarType, TimeType > >
                        derivedObservationSimulator =
                                observation_models::getObservationSimulatorOfType< 3 >( observationSimulators, observableType );

                if( derivedObservationSimulator == nullptr )
                {
                    throw std::runtime_error( "Error when simulating observation: dynamic cast to size 3 is nullptr" );
                }

                // Simulate observations for current observable and link ends set.
                observationDataset->addObservationSetFromDataset(
                        *simulateSingleObservationDataset< ObservationScalarType, TimeType, 3 >(
                                observationsToSimulate.at( i ), derivedObservationSimulator, bodies ),
                        0 );

                break;
            }
            default:
                throw std::runtime_error( "Error, simulation of observations not yet implemented for size " +
                                          std::to_string( observationSize ) );
        }
    }
    return observationDataset;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > simulateObservations(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationsToSimulate,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulators,
        const SystemOfBodies bodies )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            simulateObservationDataset< ObservationScalarType, TimeType >( observationsToSimulate, observationSimulators, bodies ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > setExistingObservationDataset(
        const std::map< observation_models::ObservableType,
                        std::pair< observation_models::LinkEnds,
                                   std::pair< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >,
                                              std::vector< TimeType > > > > observationsInput,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::map< observation_models::ObservableType, std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > >
                ancillarySettings = std::map< observation_models::ObservableType,
                                              std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > >( ) )
{
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset =
            std::make_shared< observation_models::ObservationDataset< ObservationScalarType, TimeType > >( );

    // Iterate over all observables.
    for( auto itr : observationsInput )
    {
        observation_models::ObservableType observableType = itr.first;
        observation_models::LinkEnds linkEnds = itr.second.first;

        std::pair< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >, std::vector< TimeType > >
                observationsTimesAndValues = itr.second.second;

        if( observationsTimesAndValues.first.size( ) != observationsTimesAndValues.second.size( ) )
        {
            throw std::runtime_error(
                    "Error when setting observation collection from existing observations, size of observation times and observation "
                    "values is inconsistent." );
        }

        std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > currentAncillarySettings = nullptr;
        if( ancillarySettings.count( observableType ) )
        {
            currentAncillarySettings = ancillarySettings.at( observableType );
        }

        observationDataset->addObservationSet( observableType,
                                               linkEnds,
                                               observationsTimesAndValues.first,
                                               observationsTimesAndValues.second,
                                               referenceLinkEnd,
                                               std::vector< Eigen::VectorXd >( ),
                                               nullptr,
                                               currentAncillarySettings,
                                               std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                               std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                                               true );
    }

    return observationDataset;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > setExistingObservations(
        const std::map< observation_models::ObservableType,
                        std::pair< observation_models::LinkEnds,
                                   std::pair< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >,
                                              std::vector< TimeType > > > > observationsInput,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::map< observation_models::ObservableType, std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > >
                ancillarySettings = std::map< observation_models::ObservableType,
                                              std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > >( ) )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            setExistingObservationDataset< ObservationScalarType, TimeType >( observationsInput, referenceLinkEnd, ancillarySettings ) );
}

Eigen::VectorXd getIdenticallyAndIndependentlyDistributedNoise( const std::function< double( const double ) > noiseFunction,
                                                                const int observationSize,
                                                                const double evaluationTime );

//! Function to remove link id from the simulated observations
/*!
 * /param simulatedObservations The simulated observation
 * /return Simulated observations without link end id
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::map< observation_models::ObservableType,
          std::map< observation_models::LinkEnds,
                    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > > >
removeLinkIdFromSimulatedObservations(
        std::map< observation_models::ObservableType,
                  std::map< observation_models::LinkEnds,
                            std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
                                       std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > > simulatedObservations )
{
    std::map< observation_models::ObservableType,
              std::map< observation_models::LinkEnds,
                        std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > > >
            observationsWithoutLinkEndId;

    for( typename std::map< observation_models::ObservableType,
                            std::map< observation_models::LinkEnds,
                                      std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
                                                 std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > >::const_iterator
                 observationIterator = simulatedObservations.begin( );
         observationIterator != simulatedObservations.end( );
         observationIterator++ )
    {
        for( typename std::map< observation_models::LinkEnds,
                                std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
                                           std::pair< std::vector< TimeType >, observation_models::LinkEndType > > >::const_iterator
                     linkIterator = observationIterator->second.begin( );
             linkIterator != observationIterator->second.end( );
             linkIterator++ )
        {
            observationsWithoutLinkEndId[ observationIterator->first ][ linkIterator->first ] =
                    std::make_pair( linkIterator->second.first, linkIterator->second.second.first );
        }
    }
    return observationsWithoutLinkEndId;
}

std::map< double, Eigen::VectorXd > getTargetAnglesAndRange( const simulation_setup::SystemOfBodies& bodies,
                                                             const std::pair< std::string, std::string > groundStationId,
                                                             const std::string& targetBody,
                                                             const std::vector< double > times,
                                                             const bool transmittingToTarget );

/*!
 * Creates observation simulation settings to be used when simulating observations consistent with an observed observation
 * collection.
 *
 * @param observedObservationCollection Observation collection of observed observations (i.e. from ODF file)
 * @return Observation simulation settings for simulated observations.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > >
getObservationSimulationSettingsFromObservationDataset(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observedObservationDataset,
        const SystemOfBodies& bodies )
{
    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > > observationSimulationSettings;

    for( unsigned int setId = 0; setId < observedObservationDataset->getNumberOfObservationSets( ); ++setId )
    {
        const observation_models::ObservationSetMetadata< ObservationScalarType, TimeType >& metadata =
                observedObservationDataset->getObservationSetMetadata( setId );
        const std::shared_ptr< simulation_setup::TabulatedObservationSimulationSettings< TimeType > > singleSetSimulationSettings =
                std::make_shared< simulation_setup::TabulatedObservationSimulationSettings< TimeType > >(
                        metadata.observableType_,
                        observedObservationDataset->getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_,
                        observedObservationDataset->getObservationTimesForSet( setId ),
                        metadata.referenceLinkEnd_,
                        std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
                        nullptr,
                        observedObservationDataset->getAncillarySettings( metadata.ancillarySettingsId_ ) );

        const std::shared_ptr< ObservationDependentVariableBookkeeping > dependentVariableBookkeeping =
                observedObservationDataset->getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
        if( dependentVariableBookkeeping != nullptr )
        {
            const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariablesList =
                    dependentVariableBookkeeping->getDependentVariableSettings( );
            if( dependentVariablesList.size( ) > 0 )
            {
                singleSetSimulationSettings->getObservationDependentVariableBookkeeping( )->addDependentVariables( dependentVariablesList );
            }
        }

        observationSimulationSettings.push_back( singleSetSimulationSettings );
    }

    return observationSimulationSettings;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > >
getObservationSimulationSettingsFromObservations(
        std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observedObservationCollection,
        const SystemOfBodies& bodies )
{
    return getObservationSimulationSettingsFromObservationDataset< ObservationScalarType, TimeType >(
            observedObservationCollection->getObservationDataset( ), bodies );
}

//! Function that simulates observations from a given observation collection (containing real observations), and computes the residuals and relevant observation dependent variables.
//! The derived residuals and dependent variables are then appended to the observed observation collection.
template< typename ObservationScalarType = double, typename TimeType = double >
void computeResidualsAndDependentVariables(
        std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulators,
        const SystemOfBodies& bodies )
{
    // Retrieve observation simulation settings from observed observation collection and simulate observations
    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservationDataset< ObservationScalarType, TimeType >( observationDataset, bodies );
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > computedObservationDataset =
            simulateObservationDataset( observationSimulationSettings, observationSimulators, bodies );

    const observation_models::FlattenedObservationData< ObservationScalarType, TimeType > observationData =
            observationDataset->createComputationFlattenedObservationData( true );
    const observation_models::FlattenedObservationData< ObservationScalarType, TimeType > computedObservationData =
            computedObservationDataset->createComputationFlattenedObservationData( true );

    // Residual/dependent-variable computation includes rejected rows by default so they remain inspectable for restoration decisions.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residuals =
            observationData.getObservationVector( ) - computedObservationData.getObservationVector( );
    observationDataset->setResidualVector( observationData, residuals );

    for( unsigned int setId = 0; setId < observationDataset->getNumberOfObservationSets( ); ++setId )
    {
        std::vector< Eigen::VectorXd > computedDependentVariables = computedObservationDataset->getDependentVariablesForSet( setId );
        if( computedDependentVariables.size( ) > 0 )
        {
            if( computedDependentVariables.at( 0 ).size( ) > 0 )
            {
                observationDataset->setDependentVariablesForSet( setId, computedDependentVariables );
            }
        }
    }
}

template< typename ObservationScalarType = double, typename TimeType = double >
void computeResidualsAndDependentVariables(
        std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulators,
        const SystemOfBodies& bodies )
{
    computeResidualsAndDependentVariables< ObservationScalarType, TimeType >(
            observationCollection->getObservationDataset( ), observationSimulators, bodies );
}

template< typename ObservationScalarType = double, typename TimeType = double >
Eigen::VectorXd getNumericalObservationTimePartial(
        const std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulators,
        const SystemOfBodies& bodies,
        const double timePerturbation )
{
    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > > upPerturbedObservationSimulationSettings;
    for( unsigned int i = 0; i < observationSimulationSettings.size( ); i++ )
    {
        upPerturbedObservationSimulationSettings.push_back(
                perturbObservationTime< TimeType >( observationSimulationSettings.at( i ), timePerturbation ) );
    }
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > computedUpperturbedObservationDataset =
            simulateObservationDataset( upPerturbedObservationSimulationSettings, observationSimulators, bodies );

    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > >
            downPerturbedObservationSimulationSettings;
    for( unsigned int i = 0; i < observationSimulationSettings.size( ); i++ )
    {
        downPerturbedObservationSimulationSettings.push_back(
                perturbObservationTime< TimeType >( observationSimulationSettings.at( i ), -timePerturbation ) );
    }
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > computedDownperturbedObservationDataset =
            simulateObservationDataset( downPerturbedObservationSimulationSettings, observationSimulators, bodies );

    return ( computedUpperturbedObservationDataset->createOrderedFlattenedObservationData( )
                     .getObservationVector( )
                     .template cast< double >( ) -
             computedDownperturbedObservationDataset->createOrderedFlattenedObservationData( )
                     .getObservationVector( )
                     .template cast< double >( ) ) /
            ( 2.0 * timePerturbation );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void estimateTimeBiasPerSet(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const Eigen::VectorXd& timePartials,
        std::vector< double >& timeBiases,
        Eigen::VectorXd& correctedResiduals )
{
    std::vector< std::pair< int, int > > startEndIndices = observationDataset->getObservationSetStartAndSize( );
    Eigen::VectorXd residualVector =
            observationDataset->createOrderedFlattenedObservationData( ).getResidualVector( ).template cast< double >( );
    correctedResiduals.resize( residualVector.rows( ), 1 );

    for( unsigned int i = 0; i < startEndIndices.size( ); i++ )
    {
        Eigen::VectorXd currentPartials = timePartials.segment( startEndIndices.at( i ).first, startEndIndices.at( i ).second );
        double currentTimeBias = ( ( currentPartials.transpose( ) * currentPartials ).inverse( ) *
                                   ( currentPartials.transpose( ) *
                                     residualVector.segment( startEndIndices.at( i ).first, startEndIndices.at( i ).second ) ) )( 0 );
        correctedResiduals.segment( startEndIndices.at( i ).first, startEndIndices.at( i ).second ) =
                residualVector.segment( startEndIndices.at( i ).first, startEndIndices.at( i ).second ) - currentTimeBias * currentPartials;
        timeBiases.push_back( currentTimeBias );
    }
}

template< typename ObservationScalarType = double, typename TimeType = double >
void estimateTimeBiasPerSet(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const Eigen::VectorXd& timePartials,
        std::vector< double >& timeBiases,
        Eigen::VectorXd& correctedResiduals )
{
    estimateTimeBiasPerSet< ObservationScalarType, TimeType >(
            observationCollection->getObservationDataset( ), timePartials, timeBiases, correctedResiduals );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void estimateTimeBiasAndPolynomialFitPerSet(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const Eigen::VectorXd& timePartials,
        std::vector< double >& timeBiases,
        std::vector< Eigen::VectorXd >& polynomialCoefficientsList,
        Eigen::VectorXd& correctedResiduals )
{
    estimateTimeBiasPerSet( observationDataset, timePartials, timeBiases, correctedResiduals );

    std::vector< double > stlTimeVector =
            utilities::staticCastVector< double, TimeType >( observationDataset->createOrderedFlattenedObservationData( ).getTimes( ) );
    Eigen::VectorXd timeVector = utilities::convertStlVectorToEigenVector< double >( stlTimeVector );

    std::vector< std::pair< int, int > > startEndIndices = observationDataset->getObservationSetStartAndSize( );

    //    for( unsigned int i = 0; i < startEndIndices.size( ); i++ )
    //    {
    //        Eigen::VectorXd currentTimes =
    //            timeVector.segment( startEndIndices.at( i ).first, startEndIndices.at( i ).second ).array( ) - timeVector(
    //            startEndIndices.at( i ).first );
    //        Eigen::VectorXd currentResiduals = correctedResiduals.segment( startEndIndices.at( i ).first, startEndIndices.at( i ).second
    //        ); Eigen::VectorXd polynomialCoefficients = linear_algebra::getLeastSquaresPolynomialFit(
    //            currentTimes, currentResiduals, { 0, 1 } );
    //        polynomialCoefficientsList.push_back( polynomialCoefficients );
    //        Eigen::VectorXd polynomialValues = linear_algebra::evaluatePolynomial( currentTimes, polynomialCoefficients, { 0, 1, 2, 3 } );
    //        correctedResiduals.segment( startEndIndices.at( i ).first, startEndIndices.at( i ).second ) -= polynomialValues;
    //    }
}

template< typename ObservationScalarType = double, typename TimeType = double >
void estimateTimeBiasAndPolynomialFitPerSet(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const Eigen::VectorXd& timePartials,
        std::vector< double >& timeBiases,
        std::vector< Eigen::VectorXd >& polynomialCoefficientsList,
        Eigen::VectorXd& correctedResiduals )
{
    estimateTimeBiasAndPolynomialFitPerSet< ObservationScalarType, TimeType >(
            observationCollection->getObservationDataset( ), timePartials, timeBiases, polynomialCoefficientsList, correctedResiduals );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void getResidualStatistics(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        Eigen::VectorXd& startTimes,
        Eigen::VectorXd& durations,
        Eigen::VectorXd& meanValues,
        Eigen::VectorXd& rmsValues )
{
    const observation_models::FlattenedObservationData< ObservationScalarType, TimeType > flattenedObservationData =
            observationDataset->createOrderedFlattenedObservationData( );
    std::vector< double > stlTimeVector = utilities::staticCastVector< double, TimeType >( flattenedObservationData.getTimes( ) );
    Eigen::VectorXd timeVector = utilities::convertStlVectorToEigenVector< double >( stlTimeVector );

    Eigen::VectorXd residuals = flattenedObservationData.getResidualVector( ).template cast< double >( );

    std::vector< std::pair< int, int > > startEndIndices = observationDataset->getObservationSetStartAndSize( );
    startTimes = Eigen::VectorXd::Zero( startEndIndices.size( ) );
    durations = Eigen::VectorXd::Zero( startEndIndices.size( ) );
    meanValues = Eigen::VectorXd::Zero( startEndIndices.size( ) );
    rmsValues = Eigen::VectorXd::Zero( startEndIndices.size( ) );

    for( unsigned int i = 0; i < startEndIndices.size( ); i++ )
    {
        Eigen::VectorXd currentTimes = timeVector.segment( startEndIndices.at( i ).first, startEndIndices.at( i ).second ).array( ) -
                timeVector( startEndIndices.at( i ).first );
        Eigen::VectorXd currentResiduals = residuals.segment( startEndIndices.at( i ).first, startEndIndices.at( i ).second );

        startTimes( i ) = timeVector( startEndIndices.at( i ).first );
        durations( i ) = currentTimes( currentTimes.rows( ) - 1 ) - currentTimes( 0 );
        meanValues( i ) = linear_algebra::getVectorEntryMean( currentResiduals );
        rmsValues( i ) = linear_algebra::getVectorEntryRootMeanSquare( currentResiduals );
    }
}

template< typename ObservationScalarType = double, typename TimeType = double >
void getResidualStatistics(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        Eigen::VectorXd& startTimes,
        Eigen::VectorXd& durations,
        Eigen::VectorXd& meanValues,
        Eigen::VectorXd& rmsValues )
{
    getResidualStatistics< ObservationScalarType, TimeType >(
            observationCollection->getObservationDataset( ), startTimes, durations, meanValues, rmsValues );
}

}  // namespace simulation_setup

}  // namespace tudat
#endif  // TUDAT_SIMULATEOBSERVATIONS_H
