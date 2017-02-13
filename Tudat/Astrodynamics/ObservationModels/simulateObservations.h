/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/ObservationModels/observationManager.h"

namespace tudat
{

namespace observation_models
{

//! Base struct for defining times at which observations are to be simulated.
/*!
 *  Base struct for defining times at which observations are to be simulated. Here, only the link end from which the
 *  observation is to be calculated is defined. Derived classes are used for defining the times themselves
 *  (either directly or through some algorithm).
 */
template< typename TimeType >
struct ObservationSimulationTimeSettings
{
    //! Constructor, defines link end type.
    /*!
     *  Constructor, defines link end type from which observations are to be simulated.
     *  \param linkEndType Link end type from which observations are to be simulated.
     */
    ObservationSimulationTimeSettings( const LinkEndType linkEndType ):linkEndType_( linkEndType ){ }

    //! Destructor.
    virtual ~ObservationSimulationTimeSettings( ){ }

    //! Link end type from which observations are to be simulated.
    LinkEndType linkEndType_;
};

template< typename TimeType >
struct TabulatedObservationSimulationTimeSettings: public ObservationSimulationTimeSettings< TimeType >
{
    TabulatedObservationSimulationTimeSettings(
            const LinkEndType linkEndType, const std::vector< TimeType >& simulationTimes ): ObservationSimulationTimeSettings< TimeType >( linkEndType ),
        simulationTimes_( simulationTimes ){ }

    ~TabulatedObservationSimulationTimeSettings( ){ }

    std::vector< TimeType > simulationTimes_;
};

//! Function to compute observations at times defined by settings object using a given observation model
/*!
 *  Function to compute observations at times defined by settings object using a given observation model
 *  \param observationsToSimulate Object that computes/defines settings for observation times/reference link end
 *  \param observationModel Observation model that is to be used to compute observations
 *  \return Pair of observable values and observation time (with associated reference link end)
 */
template< typename ObservationScalarType = double, typename TimeType = double,
          int ObservationSize = 1 >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,std::pair< std::vector< TimeType >, LinkEndType > >
simulateSingleObservationSet(
        const boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > observationsToSimulate,
        const boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel )
{
    // Delcare return type.
    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > >
            simulatedObservations;

    // Simulate observations from tabulated times.
    if( boost::dynamic_pointer_cast< TabulatedObservationSimulationTimeSettings< TimeType > >( observationsToSimulate ) != NULL )
    {
        boost::shared_ptr< TabulatedObservationSimulationTimeSettings< TimeType > > tabulatedObservationSettings =
                boost::dynamic_pointer_cast< TabulatedObservationSimulationTimeSettings< TimeType > >( observationsToSimulate );

        // Simulate observations at requested pre-defined time.
        simulatedObservations = simulateObservationsWithCheckAndLinkEndIdOutput<
                ObservationSize, ObservationScalarType, TimeType >(
                    tabulatedObservationSettings->simulationTimes_, observationModel, observationsToSimulate->linkEndType_ );

    }

    return simulatedObservations;
}

//! Function to simulate observations for single observable and single set of link ends.
/*!
 *  Function to simulate observations for single observable and single set of link ends. From the observation time settings and
 *  the observation manager (model), the required observations are simulated and returned.
 *  \param observationsToSimulate Object that computes/defines settings for observation times/reference link end
 *  \param observationManager Observation manager for observable for which observations are to be calculated.
 *  \param linkEnds Link end set for which observations are to be calculated.
 *  \return Pair of first: vector of observations; second: vector of times at which observations are taken
 *  (reference to link end defined in observationsToSimulate).
 */
template< typename ObservationScalarType = double, typename TimeType = double, int ObservationSize = 1 >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,std::pair< std::vector< TimeType >, LinkEndType > >
simulateSingleObservationSet(
        const boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > observationsToSimulate,
        const boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > observationManager,
        const LinkEnds& linkEnds )
{
    if( observationManager == NULL )
    {
        throw std::runtime_error( "Error when simulating single observation set, observation manager is NULL" );
    }

    return simulateSingleObservationSet< ObservationScalarType, TimeType, ObservationSize >(
                observationsToSimulate, observationManager->getObservationModel( linkEnds ) );
}

//! Function to generate ObservationSimulationTimeSettings objects from simple time list input.
/*!
 *  Function to generate ObservationSimulationTimeSettings objects, as required for observation simulation from
 *  simulateSingleObservationSet from simple time vectors input.
 *  \param originalMap List of observation times per link end set per observable type.
 *  \return TabulatedObservationSimulationTimeSettings objects from time list input (originalMap)
 */
template< typename TimeType = double >
std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > > >
createObservationSimulationTimeSettingsMap(
        const std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > >&
        originalMap )
{
    // Declare return map.
    std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > > >
            newMap;

    // Iterate over all observables.
    for( typename std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > >::
         const_iterator it = originalMap.begin( ); it != originalMap.end( ); it++ )
    {
        // Iterate over all link end sets
        for( typename std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > >::const_iterator
             linkEndIterator = it->second.begin( ); linkEndIterator != it->second.end( ); linkEndIterator++ )
        {
            // Create ObservationSimulationTimeSettings from vector of times.
            newMap[ it->first ][ linkEndIterator->first ] = boost::make_shared<
                    TabulatedObservationSimulationTimeSettings< TimeType > >(
                        linkEndIterator->second.second, linkEndIterator->second.first );
        }
    }
    return newMap;
}

//! Function to simulate observations from set of observables and link and sets and simple vectors of requested times.
/*!
 *  Function to simulate observations from set of observables and link and sets and simple vectors of requested times.
 *  Function calls function for creation of ObservationSimulationTimeSettings and subsequently uses these for function call
 *  to observation simulation function.
 *  \param observationsToSimulate List of observation times per link end set per observable type.
 *  \param observationManagers List of observation managers per link end set per observable type.
 *  \return Simulated observatoon values and associated times for requested observable types and link end sets.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
std::pair< std::vector< TimeType >, LinkEndType > > > >
simulateObservations(
        const std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > >&
        observationsToSimulate,
        const std::map< ObservableType, boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > >&
        observationManagers )
{
    return simulateObservations< ObservationScalarType, TimeType >(
                createObservationSimulationTimeSettingsMap( observationsToSimulate ), observationManagers );
}

//! Function to simulate observations from set of observables and link and sets
/*!
 *  Function to simulate observations from set of observables and link and setss and simple vectors of requested times.
 *  Iterates over all observables and link ends and simulates observations.
 *  \param observationsToSimulate List of observation time settings per link end set per observable type.
 *  \param observationManagers List of observation managers per link end set per observable type.
 *  \return Simulated observatoon values and associated times for requested observable types and link end sets.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
std::pair< std::vector< TimeType >, LinkEndType > > > >
simulateObservations(
        const std::map< ObservableType, std::map< LinkEnds,
        boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > > >& observationsToSimulate,
        const std::map< ObservableType,
        boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > >& observationManagers )
{
    // Declare return map.
    std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
            std::pair< std::vector< TimeType >, LinkEndType > > > > observations;

    // Iterate over all observables.
    for( typename std::map< ObservableType, std::map< LinkEnds,
         boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > >  > >::const_iterator observationIterator =
         observationsToSimulate.begin( ); observationIterator != observationsToSimulate.end( ); observationIterator++ )
    {

        // Iterate over all link ends for current observable.
        for( typename std::map< LinkEnds,
             boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > >::const_iterator linkEndIterator =
             observationIterator->second.begin( ); linkEndIterator != observationIterator->second.end( ); linkEndIterator++ )
        {
            int observationSize = observationManagers.at( observationIterator->first )->getObservationSize(
                        linkEndIterator->first );

            switch( observationSize )
            {
            case 1:
            {
                boost::shared_ptr< ObservationManager< 1, ObservationScalarType, TimeType > > derivedObservationManager =
                        boost::dynamic_pointer_cast< ObservationManager< 1, ObservationScalarType, TimeType > >(
                            observationManagers.at( observationIterator->first ) );

                if( derivedObservationManager == NULL )
                {
                    throw std::runtime_error( "Error when simulating observation: dynamic case to size 1 is NULL" );
                }

                // Simulate observations for current observable and link ends set.
                observations[ observationIterator->first ][ linkEndIterator->first ] =
                        simulateSingleObservationSet< ObservationScalarType, TimeType, 1 >(
                            linkEndIterator->second, derivedObservationManager->getObservationSimulator( ),
                            linkEndIterator->first );
                break;
            }
            case 2:
            {
                boost::shared_ptr< ObservationManager< 2, ObservationScalarType, TimeType > > derivedObservationManager =
                        boost::dynamic_pointer_cast< ObservationManager< 2, ObservationScalarType, TimeType > >(
                            observationManagers.at( observationIterator->first ) );

                if( derivedObservationManager == NULL )
                {
                    throw std::runtime_error( "Error when simulating observation: dynamic case to size 2 is NULL" );
                }

                // Simulate observations for current observable and link ends set.
                observations[ observationIterator->first ][ linkEndIterator->first ] =
                        simulateSingleObservationSet< ObservationScalarType, TimeType, 2 >(
                            linkEndIterator->second, derivedObservationManager->getObservationSimulator( ),
                            linkEndIterator->first );
                break;
            }
            case 3:
            {
                boost::shared_ptr< ObservationManager< 3, ObservationScalarType, TimeType > > derivedObservationManager =
                        boost::dynamic_pointer_cast< ObservationManager< 3, ObservationScalarType, TimeType > >(
                            observationManagers.at( observationIterator->first ) );

                if( derivedObservationManager == NULL )
                {
                    throw std::runtime_error( "Error when simulating observation: dynamic case to size 3 is NULL" );
                }

                // Simulate observations for current observable and link ends set.
                observations[ observationIterator->first ][ linkEndIterator->first ] =
                        simulateSingleObservationSet< ObservationScalarType, TimeType, 3 >(
                            linkEndIterator->second, derivedObservationManager->getObservationSimulator( ), linkEndIterator->first );
                break;
            }
            default:
                throw std::runtime_error( "Error, simulation of observations not yet implemented for size " +
                                          boost::lexical_cast< std::string >( observationSize ) );

            }
        }
    }
    return observations;
}


}

}
#endif // TUDAT_SIMULATEOBSERVATIONS_H
