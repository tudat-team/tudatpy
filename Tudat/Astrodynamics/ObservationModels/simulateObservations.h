#ifndef SIMULATEOBSERVATIONS_H
#define SIMULATEOBSERVATIONS_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/ObservationModels/observationManager.h"

namespace tudat
{

namespace observation_models
{

//! Base struct for defining times at which observations are to be simulated.
/*!
 *  Base struct for defining times at which observations are to be simulated. Here, only the link end from which the observation
 *  is to be calculated is defined. Derived classes are used for defining the times themselves (either directly or through some algorithm).
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

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType,
          int ObservationSize = 1 >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,std::pair< std::vector< TimeType >, LinkEndType > > simulateSingleObservationSet(
        const boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > observationsToSimulate,
        const boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > observationModel )
{
    // Delcare return type.
    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > > simulatedObservations;

    // Simulate observations from tabulated times.
    if( boost::dynamic_pointer_cast< TabulatedObservationSimulationTimeSettings< TimeType > >( observationsToSimulate ) != NULL )
    {
        boost::shared_ptr< TabulatedObservationSimulationTimeSettings< TimeType > > tabulatedObservationSettings =
                boost::dynamic_pointer_cast< TabulatedObservationSimulationTimeSettings< TimeType > >( observationsToSimulate );

        // Simulate observations at requested pre-defined time.
        simulatedObservations = simulateObservationsWithCheckAndLinkEndIdOutput<
                ObservationSize, ObservationScalarType, TimeType, StateScalarType >(
                    tabulatedObservationSettings->simulationTimes_, observationModel, observationsToSimulate->linkEndType_ );
    }

    return simulatedObservations;
}

//! Function to simulate observations for single observable and single set of link ends.
/*!
 *  Function to simulate observations for single observable and single set of link ends. From the observation time settings and
 *  the observation manager (model), the required observations are simulated and returned.
 *  \param observationsToSimulate Settings for times (and link end types) at which observations are to be calculated.
 *  \param observationManager Observation manager for observable for which observations are to be calculated.
 *  \param linkEnds Link end set for which observations are to be calculated.
 *  \return Pair of first: vector of observations; second: vector of times at which observations are taken (reference to link end defined in
 *  observationsToSimulate).
 */
template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType, int ObservationSize = 1 >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,std::pair< std::vector< TimeType >, LinkEndType > > simulateSingleObservationSet(
        const boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > observationsToSimulate,
        const boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > observationManager,
        const LinkEnds& linkEnds )
{
    if( observationManager == NULL )
    {
        std::cerr<<"Error when simulating single observation set, observable size is not 1 when making dynamic cast"<<std::endl;
    }

    return simulateSingleObservationSet< ObservationScalarType, TimeType, StateScalarType, ObservationSize >(
                observationsToSimulate, observationManager->getObservationModel( linkEnds ) );
}

//! Function to generate ObservationSimulationTimeSettings objects from simple time list input.
/*!
 *  Function to generate ObservationSimulationTimeSettings objects, as required for observation simulation from simulateSingleObservationSet,
 *  from simple time vectors input.
 *  \param originalMap List of observation times per link end set per observable type.
 *  \return TabulatedObservationSimulationTimeSettings objects from simple time list input (originalMap)
 */
template< typename TimeType = double >
std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > > >
createObservationSimulationTimeSettingsMap(
        const std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > >& originalMap )
{
    // Declare return map.
    std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > > > newMap;

    // Iterate over all observables.
    for( typename std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > >::const_iterator it =
         originalMap.begin( ); it != originalMap.end( ); it++ )
    {
        // Iterate over all link end sets
        for( typename std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > >::const_iterator linkEndIterator = it->second.begin( );
             linkEndIterator != it->second.end( ); linkEndIterator++ )
        {
            // Create ObservationSimulationTimeSettings from vector of times.
            newMap[ it->first ][ linkEndIterator->first ] = boost::make_shared< TabulatedObservationSimulationTimeSettings< TimeType > >(
                        linkEndIterator->second.second, linkEndIterator->second.first );
        }
    }
    return newMap;
}

//! Function to simulate observations from set of observables and link and sets and simple vectors of requested times.
/*!
 *  Function to simulate observations from set of observables and link and sets and simple vectors of requested times. Function calls function
 *  for creation of ObservationSimulationTimeSettings and subsequently uses these for function call to observation simulation function.
 *  \param observationsToSimulate List of observation times per link end set per observable type.
 *  \param observationManagers List of observation managers per link end set per observable type.
 *  \return Simulated observatoon values and associated times for requested observable types and link end sets.
 */
template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
std::pair< std::vector< TimeType >, LinkEndType > > > >
simulateObservations(
        const std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > >& observationsToSimulate,
        const std::map< ObservableType, boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > >& observationManagers )
{
    return simulateObservations< ObservationScalarType, TimeType, StateScalarType >(
                createObservationSimulationTimeSettingsMap( observationsToSimulate ), observationManagers );
}

//! Function to simulate observations from set of observables and link and sets
/*!
 *  Function to simulate observations from set of observables and link and setss and simple vectors of requested times. Iterates over all
 *  observables and link ends and simulates observations.
 *  \param observationsToSimulate List of observation time settings per link end set per observable type.
 *  \param observationManagers List of observation managers per link end set per observable type.
 *  \return Simulated observatoon values and associated times for requested observable types and link end sets.
 */
template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
std::pair< std::vector< TimeType >, LinkEndType > > > >
simulateObservations(
        const std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > > >& observationsToSimulate,
        const std::map< ObservableType, boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > >& observationManagers )
{
    // Declare return map.
    std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
            std::pair< std::vector< TimeType >, LinkEndType > > > > observations;

    // Iterate over all observables.
    for( typename std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > >  > >::const_iterator observationIterator =
         observationsToSimulate.begin( ); observationIterator != observationsToSimulate.end( ); observationIterator++ )
    {
        //std::cout<<"Simulating observable of type: ";
        //printObservableType( observationIterator->first );
        std::cout<<std::endl;

        // Iterate over all link ends for current observable.
        for( typename std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > >::const_iterator linkEndIterator =
             observationIterator->second.begin( ); linkEndIterator != observationIterator->second.end( ); linkEndIterator++ )
        {
            //std::cout<<"Link ends: ";
            //printLinkEnds( linkEndIterator->first );

            //std::cout<<"Obs: "<<observationManagers.count( observationIterator->first )<<std::endl;

            int observationSize = observationManagers.at( observationIterator->first )->getObservationSize( linkEndIterator->first );

            //std::cout<<"Obs. size: "<<observationSize<<std::endl;
            switch( observationSize )
            {
            case 1:
            {
                boost::shared_ptr< ObservationManager< 1, ObservationScalarType, TimeType, StateScalarType > > derivedObservationManager =
                        boost::dynamic_pointer_cast< ObservationManager< 1, ObservationScalarType, TimeType, StateScalarType > >(
                            observationManagers.at( observationIterator->first ) );

                if( derivedObservationManager == NULL )
                {
                    std::cerr<<"Error when simulating observation: dynamic case to size 1 is NULL"<<std::endl;
                }

                // Simulate observations for current observable and link ends set.
                observations[ observationIterator->first ][ linkEndIterator->first ] =
                        simulateSingleObservationSet< ObservationScalarType, TimeType, StateScalarType, 1 >(
                            linkEndIterator->second, derivedObservationManager->getObservationSimulator( ), linkEndIterator->first );
                break;
            }
            case 2:
            {
                std::cout<<"Derived obs. set"<<std::endl;
                boost::shared_ptr< ObservationManager< 2, ObservationScalarType, TimeType, StateScalarType > > derivedObservationManager =
                        boost::dynamic_pointer_cast< ObservationManager< 2, ObservationScalarType, TimeType, StateScalarType > >(
                            observationManagers.at( observationIterator->first ) );
                std::cout<<"Derived obs. set"<<( derivedObservationManager == NULL )<<
                           std::endl;

                if( derivedObservationManager == NULL )
                {
                    std::cerr<<"Error when simulating observation: dynamic case to size 2 is NULL"<<std::endl;
                }

                // Simulate observations for current observable and link ends set.
                observations[ observationIterator->first ][ linkEndIterator->first ] =
                        simulateSingleObservationSet< ObservationScalarType, TimeType, StateScalarType, 2 >(
                            linkEndIterator->second, derivedObservationManager->getObservationSimulator( ), linkEndIterator->first );
                break;
            }
            case 3:
            {
                boost::shared_ptr< ObservationManager< 3, ObservationScalarType, TimeType, StateScalarType > > derivedObservationManager =
                        boost::dynamic_pointer_cast< ObservationManager< 3, ObservationScalarType, TimeType, StateScalarType > >(
                            observationManagers.at( observationIterator->first ) );

                if( derivedObservationManager == NULL )
                {
                    std::cerr<<"Error when simulating observation: dynamic case to size 3 is NULL"<<std::endl;
                }

                // Simulate observations for current observable and link ends set.
                observations[ observationIterator->first ][ linkEndIterator->first ] =
                        simulateSingleObservationSet< ObservationScalarType, TimeType, StateScalarType, 3 >(
                            linkEndIterator->second, derivedObservationManager->getObservationSimulator( ), linkEndIterator->first );
                break;
            }
            default:
                std::cerr<<"Error, simulation of observations not yet implemented for size "<<observationSize<<std::endl;
            }


            std::cout<<". Number of observations: "<<observations[ observationIterator->first ][ linkEndIterator->first ].second.first.size( )<<std::endl;
        }
    }
    return observations;
}


}

}
#endif // SIMULATEOBSERVATIONS_H
