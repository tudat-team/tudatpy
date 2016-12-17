#ifndef OBSERVATIONSIMULATOR_H
#define OBSERVATIONSIMULATOR_H

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_models
{


template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double,
          typename StateScalarType = ObservationScalarType >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, bool > simulateObservationWithCheck(
        const TimeType& observationTime,
        const boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > observationModel,
        const LinkEndType linkEndAssociatedWithTime )
{
    // Initialize vector with reception times.
    bool isObservationFeasible = 1;

    std::vector< basic_mathematics::Vector6d > vectorOfStates;
    std::vector< double > vectorOfTimes;

    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > calculatedObservation =
            observationModel->computeObservationsWithLinkEndData(
                observationTime, linkEndAssociatedWithTime, vectorOfTimes, vectorOfStates );

    // Return pair of simulated ranges and reception times.
    return std::make_pair( calculatedObservation, isObservationFeasible );
}

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double,
          typename StateScalarType = ObservationScalarType >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > >
simulateObservationsWithCheck(
        const std::vector< TimeType >& observationTimes,
        const boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > observationModel,
        const LinkEndType linkEndAssociatedWithTime )
{
    std::map< TimeType, Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > > observations;
    std::pair< Eigen::Matrix< ObservationScalarType, ObservationSize, 1 >, bool > simulatedObservation;

    for( unsigned int i = 0; i < observationTimes.size( ); i++ )
    {
        simulatedObservation = simulateObservationWithCheck< ObservationSize, ObservationScalarType, TimeType, StateScalarType >(
                    observationTimes.at( i ), observationModel, linkEndAssociatedWithTime );

        // Check if receiving station can view transmitting station.
        if( simulatedObservation.second )
        {
            // If visible, add range and time to vecots of simulated data.
            observations[ observationTimes[ i ]  ] = simulatedObservation.first;
        }
    }

    // Return pair of simulated ranges and reception times.
    return std::make_pair( utilities::createConcatenatedEigenMatrixFromMapValues( observations ),
                           utilities::createVectorFromMapKeys( observations ) );
}

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double,
          typename StateScalarType = ObservationScalarType >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > >
simulateObservationsWithCheckAndLinkEndIdOutput(
        const std::vector< TimeType >& observationTimes,
        const boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > observationModel,
        const LinkEndType linkEndAssociatedWithTime )
{
    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > simulatedObservations =
            simulateObservationsWithCheck( observationTimes, observationModel, linkEndAssociatedWithTime );

    return std::make_pair( simulatedObservations.first, std::make_pair( simulatedObservations.second, linkEndAssociatedWithTime ) );
}



template< int ObservationSize, typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class ObservationSimulator
{
public:

    ObservationSimulator(
            const ObservableType observableType,
            const std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > >& observationModels ):
        observableType_( observableType ), observationModels_( observationModels ){ }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    virtual ~ObservationSimulator( ){ }

    ObservableType getObservableType( )
    {
        return observableType_;
    }

    int getObservationSize( const LinkEnds& linkEnds )
    {
        return observationModels_.at( linkEnds )->getObservationSize( );
    }

    boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > getObservationModel(
            const LinkEnds linkEnds )
    {
        if( observationModels_.count( linkEnds ) == 0 )
        {
            std::cerr<<"Error in observation manager when getting observation model, did not find model for given link ends "<<observationModels_.size( )<<std::endl;
        }
        return observationModels_.at( linkEnds );
    }

    std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > > getObservationModels( )
    {
        return observationModels_;
    }


    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, bool >
    simulateObservation( const TimeType& observationTime,
                         const LinkEnds linkEnds,
                         const LinkEndType linkEndAssociatedWithTime,
                         const bool checkTimes = true )
    {
        boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > selectedObservationModel =
                observationModels_.at( linkEnds );

        return simulateObservationWithCheck< ObservationSize, ObservationScalarType, TimeType, StateScalarType >(
                    observationTime, selectedObservationModel, linkEndAssociatedWithTime );
    }

    //! Function to simulate observations between specified link ends.
    /*!
     *  Function to simulate observations between specified link ends. Users can specify whether to check for availability of
     *  link at given reception time.
     *  \param recetionTimes Vector of times at which observations taked place (i.e. reception time)
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     *  \param checkReceptionTimes Boolean whether the availability of the link should be checked.
     *  \return Pair of observable values and associated observation times.
     */
    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > >
    simulateObservations( const std::vector< TimeType >& observationTimes,
                          const LinkEnds linkEnds,
                          const LinkEndType linkEndAssociatedWithTime,
                          const bool checkTimes = true )

    {
        if( observationModels_.count( linkEnds ) == 0 )
        {
            std::cerr<<"Error when simulating observtions, could not find observation model for given linke ends"<<std::endl;
        }

        boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > selectedObservationModel =
                observationModels_.at( linkEnds );

        return simulateObservationsWithCheck< ObservationSize, ObservationScalarType, TimeType, StateScalarType >(
                    observationTimes, selectedObservationModel, linkEndAssociatedWithTime );
    }

    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > >
    simulateObservationsWithLinkEndId( const std::vector< TimeType >& times,
                                       const LinkEnds linkEnds,
                                       const LinkEndType linkEndAssociatedWithTime,
                                       const bool checkTimes = true )
    {
        std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > simulatedObservations =
                simulateObservations( times, linkEnds, linkEndAssociatedWithTime, checkTimes );

        return std::make_pair( simulatedObservations.first, std::make_pair( simulatedObservations.second, linkEndAssociatedWithTime ) );
    }

protected:

    ObservableType observableType_;

    std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > > observationModels_;
};

}

}
#endif // OBSERVATIONSIMULATOR_H
