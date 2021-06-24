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

#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/statistics/randomVariableGenerator.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Base struct for defining times at which observations are to be simulated.
/*!
 *  Base struct for defining times at which observations are to be simulated. Here, only the link end from which the
 *  observation is to be calculated is defined. Derived classes are used for defining the times themselves
 *  (either directly or through some algorithm).
 */
template< typename TimeType = double >
struct ObservationSimulationSettings
{
    //! Constructor, defines link end type.
    /*!
     *  Constructor, defines link end type from which observations are to be simulated.
     *  \param linkEndType Link end type from which observations are to be simulated.
     */
    ObservationSimulationSettings(
            const observation_models::ObservableType observableType,
            const observation_models::LinkEnds linkEnds,
            const observation_models::LinkEndType linkEndType = observation_models::receiver,
            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
            std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ) ):
        observableType_( observableType ), linkEnds_( linkEnds ), linkEndType_( linkEndType ),
        viabilitySettingsList_( viabilitySettingsList ){ }

    //! Destructor.
    virtual ~ObservationSimulationSettings( ){ }

    observation_models::ObservableType observableType_;

    observation_models::LinkEnds linkEnds_;

    //! Link end type from which observations are to be simulated.
    observation_models::LinkEndType linkEndType_;

    std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > viabilitySettingsList_;

    std::function< Eigen::VectorXd( const double ) > observationNoiseFunction_;
};


//! Struct to define a list of observation times, fully defined before simulating the observations
/*!
 *  Struct to define a list of observation times, fully defined before simulating the observations. Simulations are simulated
 *  at the times stored in this struct. Some may be discarded due to the use of vaibility settins
 */
template< typename TimeType = double >
struct TabulatedObservationSimulationSettings: public ObservationSimulationSettings< TimeType >
{
    //! Constructor
    /*!
     * Constructor
     * \param linkEndType Link end type from which observations are to be simulated.
     * \param simulationTimes List of times at which to perform the observation simulation
     */
    TabulatedObservationSimulationSettings(
            const observation_models::ObservableType observableType,
            const observation_models::LinkEnds linkEnds,
            const std::vector< TimeType >& simulationTimes,
            const observation_models::LinkEndType linkEndType = observation_models::receiver,
            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
            std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ) ):
        ObservationSimulationSettings< TimeType >( observableType, linkEnds, linkEndType, viabilitySettingsList ),
        simulationTimes_( simulationTimes ){ }

    //! Destructor
    ~TabulatedObservationSimulationSettings( ){ }

    //! List of times at which to perform the observation simulation
    std::vector< TimeType > simulationTimes_;
};


//! Function to simulate a fixed number of simulations, in an arcwise manner, taking into account viability settings
/*!
 *  Function to simulate a fixed number of simulations, in an arcwise manner, taking into account viability settings. This class
 *  defines an observation arc starting every X seconds, with observations simulated every M seconds from the start of this
 *  interval until N observations have been simulated, taking into account that some are discarded due to vaibiloty settings.
 */
template< typename TimeType = double >
struct ArcLimitedObservationSimulationSettings: public ObservationSimulationSettings< TimeType >
{
    //! Constructor
    /*!
     * Constructor
     * \param linkEndType Reference link end type
     * \param startTime Time at which to start simulating observations
     * \param endTime Time at which to end simulating observations
     * \param observationInterval Time between two subsequent observations (e.g. integration time)
     * \param arcDuration Duration of an arc over which observationLimitPerArc observations are to be simulated
     * \param observationLimitPerArc Number of observations that are to be simulated per arc
     */
    ArcLimitedObservationSimulationSettings(
            const observation_models::ObservableType observableType,
            const observation_models::LinkEnds linkEnds,
            const TimeType startTime, const TimeType endTime, const TimeType observationInterval,
            const TimeType arcDuration, const int observationLimitPerArc,
            const observation_models::LinkEndType linkEndType = observation_models::receiver,
            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
            std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ) ):
        ObservationSimulationSettings< TimeType >(
            observableType, linkEnds, linkEndType, viabilitySettingsList ),
        startTime_( startTime ), endTime_( endTime ), observationInterval_( observationInterval ),
        arcDuration_( arcDuration ), observationLimitPerArc_( observationLimitPerArc ){ }

    ~ArcLimitedObservationSimulationSettings( ){ }

    //! Time at which to start simulating observations
    TimeType startTime_;

    //! Time at which to end simulating observations
    TimeType endTime_;

    //! Time between two subsequent observations (e.g. integration time)
    TimeType observationInterval_;

    //! Duration of an arc over which observationLimitPerArc observations are to be simulated
    TimeType arcDuration_;

    //! Number of observations that are to be simulated per arc
    int observationLimitPerArc_;
};

template< typename TimeType >
std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >  getObservationSimulationSettings(
    const std::map< observation_models::ObservableType, std::vector< observation_models::LinkEnds > >& linkEndsPerObservable,
        const std::vector< TimeType >& observationTimes,
        const observation_models::LinkEndType referenceLinkEnd = observation_models::receiver )
{
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;
    for( auto it : linkEndsPerObservable )
    {
        observation_models::ObservableType currentObservable = it.first;
        std::vector< observation_models::LinkEnds > currentLinkEndsList = it.second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                            currentObservable, currentLinkEndsList.at( i ), observationTimes, referenceLinkEnd ) );
        }
    }
    return measurementSimulationInput;
}

//! Function to simulate an observable, checking whether it is viable according to settings passed to this function
/*!
 *  Function to simulate an observable, checking whether it is viable according to settings passed to this function
 *  (if viability calculators are passed to this function).
 *  \param observationTime Time at which observable is to be computed
 *  \param observationModel Model used to compute observable
 *  \param linkEndAssociatedWithTime Model Reference link end for observable
 *  \param linkViabilityCalculators List of observation viability calculators, which are used to reject simulated
 *  observation if they dont fulfill a given (set of) conditions, e.g. minimum elevation angle (default none).
 *  \return Observation at given time.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, bool > simulateObservationWithCheck(
        const TimeType& observationTime,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const observation_models::LinkEndType linkEndAssociatedWithTime,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > linkViabilityCalculators =
        std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > >( ),
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr )
{
    // Initialize vector with reception times.
    bool observationFeasible = 1;

    std::vector< Eigen::Vector6d > vectorOfStates;
    std::vector< double > vectorOfTimes;

    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > calculatedObservation =
            observationModel->computeObservationsWithLinkEndData(
                observationTime, linkEndAssociatedWithTime, vectorOfTimes, vectorOfStates );

    observationFeasible = isObservationViable( vectorOfStates, vectorOfTimes, linkViabilityCalculators );

    if( observationFeasible && ( noiseFunction != nullptr ) )
    {
        Eigen::VectorXd noiseToAdd = noiseFunction( observationTime );
        if( noiseToAdd.rows( ) != ObservationSize )
        {
            throw std::runtime_error( "Error wen simulating observation noise, size of noise and observable are not compatible" );
        }
        else
        {
            calculatedObservation += noiseToAdd;
        }
    }
    // Return pair of simulated ranges and reception times.
    return std::make_pair( calculatedObservation, observationFeasible );
}

//! Function to simulate observables, checking whether they are viable according to settings passed to this function
/*!
 *  Function to simulate observables, checking whether they are viable according to settings passed to this function
 *  (if viability calculators are passed to this function).
 *  \param observationTimes Times at which observables are to be computed
 *  \param observationModel Model used to compute observables
 *  \param linkEndAssociatedWithTime Model Reference link end for observables
 *  \param linkViabilityCalculators List of observation viability calculators, which are used to reject simulated
 *  observation if they dont fulfill a given (set of) conditions, e.g. minimum elevation angle (default none).
 *  \return Observations at given time (concatenated in an Eigen vector) and associated times.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > >
simulateObservationsWithCheck(
        const std::vector< TimeType >& observationTimes,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const observation_models::LinkEndType linkEndAssociatedWithTime,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > linkViabilityCalculators =
        std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > >( ),
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr )
{
    std::map< TimeType, Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > > observations;
    std::pair< Eigen::Matrix< ObservationScalarType, ObservationSize, 1 >, bool > simulatedObservation;

    for( unsigned int i = 0; i < observationTimes.size( ); i++ )
    {
        simulatedObservation = simulateObservationWithCheck< ObservationSize, ObservationScalarType, TimeType >(
                    observationTimes.at( i ), observationModel, linkEndAssociatedWithTime, linkViabilityCalculators, noiseFunction );

        // Check if receiving station can view transmitting station.
        if( simulatedObservation.second )
        {
            // If viable, add observable and time to vector of simulated data.
            observations[ observationTimes[ i ] ] = simulatedObservation.first;
        }
    }

    // Return pair of simulated ranges and reception times.
    return std::make_pair( utilities::createConcatenatedEigenMatrixFromMapValues( observations ),
                           utilities::createVectorFromMapKeys( observations ) );
}

//! Function to simulate observables, checking whether they are viable according to settings passed to this function
/*!
 *  Function to simulate observables, checking whether they are viable according to settings passed to this function
 *  (if viability calculators are passed to this function).
 *  \param observationTimes Times at which observables are to be computed
 *  \param observationModel Model used to compute observables
 *  \param linkEndAssociatedWithTime Model Reference link end for observables
 *  \param linkViabilityCalculators List of observation viability calculators, which are used to reject simulated
 *  observation if they dont fulfill a given (set of) conditions, e.g. minimum elevation angle (default none).
 *  \return Observations at given times (concatenated in an Eigen vector), with associated times and reference link end.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, observation_models::LinkEndType > >
simulateObservationsWithCheckAndLinkEndIdOutput(
        const std::vector< TimeType >& observationTimes,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const observation_models::LinkEndType linkEndAssociatedWithTime,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > linkViabilityCalculators =
        std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > >( ),
        const std::function< Eigen::VectorXd( const double ) > noiseFunction = nullptr )
{
    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > simulatedObservations =
            simulateObservationsWithCheck( observationTimes, observationModel, linkEndAssociatedWithTime, linkViabilityCalculators, noiseFunction );

    return std::make_pair( simulatedObservations.first, std::make_pair( simulatedObservations.second, linkEndAssociatedWithTime ) );
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
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,std::pair< std::vector< TimeType >, observation_models::LinkEndType > >
simulateSingleObservationSet(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationsToSimulate,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const SystemOfBodies& bodies )
{
    // Delcare return type.
    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, observation_models::LinkEndType > >
            simulatedObservations;
    //! Function to create an list of obervation viability conditions for a single set of link ends
    std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > currentObservationViabilityCalculators =
            observation_models::createObservationViabilityCalculators(
                bodies,
                observationsToSimulate->linkEnds_,
                observationsToSimulate->observableType_,
                observationsToSimulate->viabilitySettingsList_ );

    std::function< Eigen::VectorXd( const double ) > noiseFunction = observationsToSimulate->observationNoiseFunction_;

    // Simulate observations from tabulated times.
    if( std::dynamic_pointer_cast< TabulatedObservationSimulationSettings< TimeType > >( observationsToSimulate ) != nullptr )
    {
        std::shared_ptr< TabulatedObservationSimulationSettings< TimeType > > tabulatedObservationSettings =
                std::dynamic_pointer_cast< TabulatedObservationSimulationSettings< TimeType > >( observationsToSimulate );

        // Simulate observations at requested pre-defined time.
        simulatedObservations = simulateObservationsWithCheckAndLinkEndIdOutput<
                ObservationSize, ObservationScalarType, TimeType >(
                    tabulatedObservationSettings->simulationTimes_, observationModel, observationsToSimulate->linkEndType_,
                    currentObservationViabilityCalculators, noiseFunction );

    }
    // Simulate observations per arc from settings
    else if( std::dynamic_pointer_cast< ArcLimitedObservationSimulationSettings< TimeType > >( observationsToSimulate ) != NULL )
    {
        std::shared_ptr< ArcLimitedObservationSimulationSettings< TimeType > > arcLimitedObservationSettings =
                std::dynamic_pointer_cast< ArcLimitedObservationSimulationSettings< TimeType > >( observationsToSimulate );

        // Define constituents of return pair.
        std::vector< TimeType > times;
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observations;
        int observableSize = observationModel->getObservationSize( );

        // Define vector to be used for storing single arc observations.
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleArcObservations =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
                    observableSize * arcLimitedObservationSettings->observationLimitPerArc_ );

        // Set start and end times.
        TimeType currentTime = arcLimitedObservationSettings->startTime_;
        TimeType currentArcStartTime = arcLimitedObservationSettings->startTime_;

        int totalNumberOfObservations = 0;
        int numberOfObservationsInCurrentArc = 0;

        // Simulate observations arcs until provided end time.
        while( currentTime < arcLimitedObservationSettings->endTime_ - 0.1 )
        {
            // Reset variables for start of new arc.
            numberOfObservationsInCurrentArc = 0;
            singleArcObservations.setZero( );

            // Simulate observations for sinlge arc
            while( ( currentTime < currentArcStartTime + arcLimitedObservationSettings->arcDuration_ ) &&
                   ( numberOfObservationsInCurrentArc < arcLimitedObservationSettings->observationLimitPerArc_ ) &&
                   ( currentTime < arcLimitedObservationSettings->endTime_ ) )
            {
                // Simulate single observation
                std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, bool > currentObservation =
                        simulateObservationWithCheck<
                        ObservationSize, ObservationScalarType, TimeType >(
                            currentTime, observationModel, observationsToSimulate->linkEndType_,
                            currentObservationViabilityCalculators, noiseFunction );

                // If observation is possible, set it in current arc observations.
                if( currentObservation.second )
                {
                    times.push_back( currentTime );
                    singleArcObservations.segment( numberOfObservationsInCurrentArc * observableSize, observableSize ) = currentObservation.first;
                    numberOfObservationsInCurrentArc += observableSize;
                }
                currentTime += arcLimitedObservationSettings->observationInterval_;
            }

            // Add single arc observations to total observations.
            observations.conservativeResize( totalNumberOfObservations + numberOfObservationsInCurrentArc );
            observations.segment( totalNumberOfObservations, numberOfObservationsInCurrentArc ) = singleArcObservations.segment( 0, numberOfObservationsInCurrentArc );
            totalNumberOfObservations += numberOfObservationsInCurrentArc;

            // Update times to next arc
            currentArcStartTime += arcLimitedObservationSettings->arcDuration_;
            currentTime = currentArcStartTime;

        }

        simulatedObservations = std::make_pair( observations, std::make_pair( times, arcLimitedObservationSettings->linkEndType_ ) );
    }

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
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,std::pair< std::vector< TimeType >, observation_models::LinkEndType > >
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
                observationsToSimulate, observationSimulator->getObservationModel( observationsToSimulate->linkEnds_ ),
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
std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > >
simulateObservations(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationsToSimulate,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >& observationSimulators,
        const SystemOfBodies bodies )
{
    // Declare return map.
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
            std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > > observations;

    // Iterate over all observables.
    for( unsigned int i = 0; i < observationsToSimulate.size( ); i++ )
    {
        observation_models::ObservableType observableType = observationsToSimulate.at( i )->observableType_;
        observation_models::LinkEnds linkEnds = observationsToSimulate.at( i )->linkEnds_;

        int observationSize = observation_models::getObservableSize( observableType );

        switch( observationSize )
        {
        case 1:
        {
            std::shared_ptr< observation_models::ObservationSimulator< 1, ObservationScalarType, TimeType > > derivedObservationSimulator =
                    observation_models::getObservationSimulatorOfType< 1 >( observationSimulators, observableType );

            if( derivedObservationSimulator == nullptr )
            {
                throw std::runtime_error( "Error when simulating observation: dynamic case to size 1 is nullptr" );
            }

            // Simulate observations for current observable and link ends set.
            observations[ observableType ][ linkEnds ] =
                    simulateSingleObservationSet< ObservationScalarType, TimeType, 1 >(
                        observationsToSimulate.at( i ), derivedObservationSimulator, bodies );
            break;
        }
        case 2:
        {
            std::shared_ptr< observation_models::ObservationSimulator< 2, ObservationScalarType, TimeType > > derivedObservationSimulator =
                    observation_models::getObservationSimulatorOfType< 2 >( observationSimulators, observableType );

            if( derivedObservationSimulator == nullptr )
            {
                throw std::runtime_error( "Error when simulating observation: dynamic case to size 2 is nullptr" );
            }

            // Simulate observations for current observable and link ends set.
            observations[ observableType ][ linkEnds ] =
                    simulateSingleObservationSet< ObservationScalarType, TimeType, 2 >(
                        observationsToSimulate.at( i ), derivedObservationSimulator, bodies );
            break;
        }
        case 3:
        {
            std::shared_ptr< observation_models::ObservationSimulator< 3, ObservationScalarType, TimeType > > derivedObservationSimulator =
                    observation_models::getObservationSimulatorOfType< 3 >( observationSimulators, observableType );

            if( derivedObservationSimulator == nullptr )
            {
                throw std::runtime_error( "Error when simulating observation: dynamic case to size 3 is nullptr" );
            }

            // Simulate observations for current observable and link ends set.
            observations[ observableType ][ linkEnds ] =
                    simulateSingleObservationSet< ObservationScalarType, TimeType, 3 >(
                        observationsToSimulate.at( i ), derivedObservationSimulator, bodies );

            break;
        }
        default:
            throw std::runtime_error( "Error, simulation of observations not yet implemented for size " +
                                      std::to_string( observationSize ) );

        }
    }
    return observations;
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

inline std::function< double( const double ) > getGaussianDistributionNoiseFunction(
        const double standardDeviation,
        const double mean = 0.0,
        const double seed = 0.0 )
{
    std::function< double( ) > inputFreeNoiseFunction = statistics::createBoostContinuousRandomVariableGeneratorFunction(
                statistics::normal_boost_distribution, { mean, standardDeviation }, seed );
    return std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                      inputFreeNoiseFunction, std::placeholders::_1 );
}

}

}
#endif // TUDAT_SIMULATEOBSERVATIONS_H
