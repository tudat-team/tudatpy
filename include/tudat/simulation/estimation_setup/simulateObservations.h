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
#include <boost/bind.hpp>
#include <functional>

#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/simulation/estimation_setup/observations.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/statistics/randomVariableGenerator.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tudat
{

namespace simulation_setup
{

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
        const std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr )
{
    // Simulate observable, and retrieve link end times and states
    std::vector< Eigen::Vector6d > vectorOfStates;
    std::vector< double > vectorOfTimes;
    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > calculatedObservation =
            observationModel->computeObservationsWithLinkEndData(
                observationTime, referenceLinkEnd, vectorOfTimes, vectorOfStates );
    Eigen::VectorXd dependentVariables = Eigen::VectorXd::Zero( 0 );
    if( dependentVariableCalculator != nullptr )
    {
        dependentVariables = dependentVariableCalculator->calculateDependentVariables(
                    vectorOfTimes, vectorOfStates, calculatedObservation.template cast< double >( ) );
    }
    // Check if observation is feasible
    bool observationFeasible = isObservationViable( vectorOfStates, vectorOfTimes, linkViabilityCalculators );

    // Add noise if needed.
    if( observationFeasible && ( noiseFunction != nullptr ) )
    {
        Eigen::VectorXd noiseToAdd = noiseFunction( observationTime );
        if( noiseToAdd.rows( ) != ObservationSize )
        {
            std::cout<<noiseToAdd.rows( )<<" "<<ObservationSize<<" "<<observationModel->getObservableType( )<<std::endl;
            throw std::runtime_error(
                        "Error wen simulating observation noise, size of noise (" + std::to_string( noiseToAdd.rows( ) ) +
                        ") and size of observable (" + std::to_string( ObservationSize ) +
                        ") are not compatible for observable type: " + observation_models::getObservableName( observationModel->getObservableType( ) ) );
        }
        else
        {
            calculatedObservation += noiseToAdd;
        }
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
        const std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr  )
{
    std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
    std::tuple< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, bool, Eigen::VectorXd > simulatedObservation;
    std::vector< Eigen::VectorXd > dependentVariables;

    for( unsigned int i = 0; i < observationTimes.size( ); i++ )
    {
        simulatedObservation = simulateObservationWithCheck< ObservationSize, ObservationScalarType, TimeType >(
                    observationTimes.at( i ), observationModel, referenceLinkEnd, linkViabilityCalculators, noiseFunction, dependentVariableCalculator );

        // Check if receiving station can view transmitting station.
        if( std::get< 1 >( simulatedObservation ) )
        {
            // If viable, add observable and time to vector of simulated data.
            observations[ observationTimes[ i ] ] = std::get< 0 >( simulatedObservation );
            dependentVariables.push_back( std::get< 2 >( simulatedObservation ) );
        }
    }

    // Return pair of simulated ranges and reception times.
    return std::make_tuple( utilities::createVectorFromMapValues( observations ),
                            utilities::createVectorFromMapKeys( observations ),
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
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > >
simulateObservationsWithCheckAndLinkEndIdOutput(
        const std::vector< TimeType >& observationTimes,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > linkViabilityCalculators =
        std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > >( ),
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr,
        const std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr )
{
    std::tuple< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >,
    std::vector< TimeType >,
    std::vector< Eigen::VectorXd > >
            simulatedObservations =
            simulateObservationsWithCheck( observationTimes, observationModel, referenceLinkEnd, linkViabilityCalculators, noiseFunction, dependentVariableCalculator );

    return std::make_shared< observation_models::SingleObservationSet< ObservationScalarType, TimeType > >(
                observationModel->getObservableType( ), observationModel->getLinkEnds( ),
                std::get< 0 >( simulatedObservations ), std::get< 1 >( simulatedObservations ), referenceLinkEnd, std::get< 2 >( simulatedObservations ),
                dependentVariableCalculator );
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
template< typename ObservationScalarType = double, typename TimeType = double,
          int ObservationSize = 1 >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > >
simulateSingleObservationSet(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const SystemOfBodies& bodies )
{
    // Delcare return type.
    std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > simulatedObservations;
    //! Function to create an list of obervation viability conditions for a single set of link ends
    std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > currentObservationViabilityCalculators =
            observation_models::createObservationViabilityCalculators(
                bodies,
                observationsToSimulate->getLinkEnds( ),
                observationsToSimulate->getObservableType( ),
                observationsToSimulate->getViabilitySettingsList( ) );

    std::function< Eigen::VectorXd( const double ) > noiseFunction = observationsToSimulate->getObservationNoiseFunction( );

    // Simulate observations from tabulated times.
    if( std::dynamic_pointer_cast< TabulatedObservationSimulationSettings< TimeType > >( observationsToSimulate ) != nullptr )
    {
        std::shared_ptr< TabulatedObservationSimulationSettings< TimeType > > tabulatedObservationSettings =
                std::dynamic_pointer_cast< TabulatedObservationSimulationSettings< TimeType > >( observationsToSimulate );

        // Simulate observations at requested pre-defined time.
        simulatedObservations = simulateObservationsWithCheckAndLinkEndIdOutput<
                ObservationSize, ObservationScalarType, TimeType >(
                    tabulatedObservationSettings->simulationTimes_, observationModel,
                    observationsToSimulate->getReferenceLinkEndType( ),
                    currentObservationViabilityCalculators, noiseFunction,
                    observationsToSimulate->getDependentVariableCalculator( ) );

    }
    //    // Simulate observations per arc from settings
    //    else if( std::dynamic_pointer_cast< ArcLimitedObservationSimulationSettings< TimeType > >( observationsToSimulate ) != NULL )
    //    {
    //        std::shared_ptr< ArcLimitedObservationSimulationSettings< TimeType > > arcLimitedObservationSettings =
    //                std::dynamic_pointer_cast< ArcLimitedObservationSimulationSettings< TimeType > >( observationsToSimulate );

    //        // Define constituents of return pair.
    //        std::vector< TimeType > times;
    //        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observations;
    //        int observableSize = observationModel->getObservationSize( );

    //        // Define vector to be used for storing single arc observations.
    //        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleArcObservations =
    //                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
    //                    observableSize * arcLimitedObservationSettings->observationLimitPerArc_ );

    //        // Set start and end times.
    //        TimeType currentTime = arcLimitedObservationSettings->startTime_;
    //        TimeType currentArcStartTime = arcLimitedObservationSettings->startTime_;

    //        int totalNumberOfObservations = 0;
    //        int numberOfObservationsInCurrentArc = 0;

    //        // Simulate observations arcs until provided end time.
    //        while( currentTime < arcLimitedObservationSettings->endTime_ - 0.1 )
    //        {
    //            // Reset variables for start of new arc.
    //            numberOfObservationsInCurrentArc = 0;
    //            singleArcObservations.setZero( );

    //            // Simulate observations for sinlge arc
    //            while( ( currentTime < currentArcStartTime + arcLimitedObservationSettings->arcDuration_ ) &&
    //                   ( numberOfObservationsInCurrentArc < arcLimitedObservationSettings->observationLimitPerArc_ ) &&
    //                   ( currentTime < arcLimitedObservationSettings->endTime_ ) )
    //            {
    //                // Simulate single observation
    //                std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, bool > currentObservation =
    //                        simulateObservationWithCheck<
    //                        ObservationSize, ObservationScalarType, TimeType >(
    //                            currentTime, observationModel, observationsToSimulate->getReferenceLinkEndType( ),
    //                            currentObservationViabilityCalculators, noiseFunction );

    //                // If observation is possible, set it in current arc observations.
    //                if( currentObservation.second )
    //                {
    //                    times.push_back( currentTime );
    //                    singleArcObservations.segment( numberOfObservationsInCurrentArc * observableSize, observableSize ) = currentObservation.first;
    //                    numberOfObservationsInCurrentArc += observableSize;
    //                }
    //                currentTime += arcLimitedObservationSettings->observationInterval_;
    //            }

    //            // Add single arc observations to total observations.
    //            observations.conservativeResize( totalNumberOfObservations + numberOfObservationsInCurrentArc );
    //            observations.segment( totalNumberOfObservations, numberOfObservationsInCurrentArc ) = singleArcObservations.segment( 0, numberOfObservationsInCurrentArc );
    //            totalNumberOfObservations += numberOfObservationsInCurrentArc;

    //            // Update times to next arc
    //            currentArcStartTime += arcLimitedObservationSettings->arcDuration_;
    //            currentTime = currentArcStartTime;

    //        }

    //        simulatedObservations = std::make_pair( observations, std::make_pair( times, arcLimitedObservationSettings->getReferenceLinkEndType( ) ) );
    //    }

    return simulatedObservations;
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
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > >
simulateSingleObservationSet(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > observationSimulator,
        const SystemOfBodies& bodies )
{
    if( observationSimulator == nullptr )
    {
        throw std::runtime_error( "Error when simulating single observation set, Observation simulator is nullptr" );
    }

    return simulateSingleObservationSet< ObservationScalarType, TimeType, ObservationSize >(
                observationsToSimulate, observationSimulator->getObservationModel( observationsToSimulate->getLinkEnds( ) ),
                bodies );
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
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > simulateObservations(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationsToSimulate,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >& observationSimulators,
        const SystemOfBodies bodies )
{
    // Declare return map.
    typename observation_models::ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets sortedObservations;

    // Iterate over all observables.
    for( unsigned int i = 0; i < observationsToSimulate.size( ); i++ )
    {
        observation_models::ObservableType observableType = observationsToSimulate.at( i )->getObservableType( );
        observation_models::LinkEnds linkEnds = observationsToSimulate.at( i )->getLinkEnds( );

        int observationSize = observation_models::getObservableSize( observableType );

        switch( observationSize )
        {
        case 1:
        {
            std::shared_ptr< observation_models::ObservationSimulator< 1, ObservationScalarType, TimeType > > derivedObservationSimulator =
                    observation_models::getObservationSimulatorOfType< 1 >( observationSimulators, observableType );

            if( derivedObservationSimulator == nullptr )
            {
                throw std::runtime_error( "Error when simulating observation: dynamic cast to size 1 is nullptr" );
            }

            // Simulate observations for current observable and link ends set.
            sortedObservations[ observableType ][ linkEnds ].push_back(
                        simulateSingleObservationSet< ObservationScalarType, TimeType, 1 >(
                            observationsToSimulate.at( i ), derivedObservationSimulator, bodies ) );
            break;
        }
        case 2:
        {
            std::shared_ptr< observation_models::ObservationSimulator< 2, ObservationScalarType, TimeType > > derivedObservationSimulator =
                    observation_models::getObservationSimulatorOfType< 2 >( observationSimulators, observableType );

            if( derivedObservationSimulator == nullptr )
            {
                throw std::runtime_error( "Error when simulating observation: dynamic cast to size 2 is nullptr" );
            }

            // Simulate observations for current observable and link ends set.
            sortedObservations[ observableType ][ linkEnds ].push_back(
                        simulateSingleObservationSet< ObservationScalarType, TimeType, 2 >(
                            observationsToSimulate.at( i ), derivedObservationSimulator, bodies ) );
            break;
        }
        case 3:
        {
            std::shared_ptr< observation_models::ObservationSimulator< 3, ObservationScalarType, TimeType > > derivedObservationSimulator =
                    observation_models::getObservationSimulatorOfType< 3 >( observationSimulators, observableType );

            if( derivedObservationSimulator == nullptr )
            {
                throw std::runtime_error( "Error when simulating observation: dynamic cast to size 3 is nullptr" );
            }

            // Simulate observations for current observable and link ends set.
            sortedObservations[ observableType ][ linkEnds ].push_back(
                        simulateSingleObservationSet< ObservationScalarType, TimeType, 3 >(
                            observationsToSimulate.at( i ), derivedObservationSimulator, bodies ) );

            break;
        }
        default:
            throw std::runtime_error( "Error, simulation of observations not yet implemented for size " +
                                      std::to_string( observationSize ) );

        }
    }
    std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection =
            std::make_shared< observation_models::ObservationCollection< ObservationScalarType, TimeType > >( sortedObservations );

    return observationCollection;
}

Eigen::VectorXd getIdenticallyAndIndependentlyDistributedNoise(
        const std::function< double( const double ) > noiseFunction,
        const int observationSize,
        const double evaluationTime );


//! Function to remove link id from the simulated observations
/*!
 * /param simulatedObservations The simulated observation
 * /return Simulated observations without link end id
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
std::vector< TimeType > > > > removeLinkIdFromSimulatedObservations(
        std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
        std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > > simulatedObservations )
{
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
            std::vector< TimeType > > > > observationsWithoutLinkEndId;

    for( typename std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
         std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > >::const_iterator observationIterator = simulatedObservations.begin( );
         observationIterator != simulatedObservations.end( ); observationIterator++ )
    {
        for( typename std::map< observation_models::LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
             std::pair< std::vector< TimeType >, observation_models::LinkEndType > > >::const_iterator linkIterator = observationIterator->second.begin( );
             linkIterator != observationIterator->second.end( ); linkIterator++ )
        {
            observationsWithoutLinkEndId[ observationIterator->first ][ linkIterator->first ] =
                    std::make_pair( linkIterator->second.first, linkIterator->second.second.first );
        }
    }
    return observationsWithoutLinkEndId;
}


std::map< double, Eigen::VectorXd > getTargetAnglesAndRange(
        const simulation_setup::SystemOfBodies& bodies,
        const std::pair< std::string, std::string > groundStationId,
        const std::string& targetBody,
        const std::vector< double > times,
        const bool transmittingToTarget );


}

}
#endif // TUDAT_SIMULATEOBSERVATIONS_H
