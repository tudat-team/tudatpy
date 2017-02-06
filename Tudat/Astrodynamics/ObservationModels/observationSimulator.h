/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONSIMULATOR_H
#define TUDAT_OBSERVATIONSIMULATOR_H

#include <boost/lexical_cast.hpp>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_models
{


//! Function to simulate an observable, checking whether it is viable according to settings passed to this function
/*!
 *  Function to simulate an observable, checking whether it is viable according to settings passed to this function
 *  NOTE: Viability check is turned of at present, will be incorporated in subsequent pull request
 *  \param observationTime Time at which observable is to be computed
 *  \param observationModel Model used to compute observable
 *  \param linkEndAssociatedWithTime Model Reference link end for observable
 *  \return Observation at given time.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, bool > simulateObservationWithCheck(
        const TimeType& observationTime,
        const boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const LinkEndType linkEndAssociatedWithTime )
{
    // Initialize vector with reception times.
    bool isObservationFeasible = 1;

    std::vector< Eigen::Vector6d > vectorOfStates;
    std::vector< double > vectorOfTimes;

    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > calculatedObservation =
            observationModel->computeObservationsWithLinkEndData(
                observationTime, linkEndAssociatedWithTime, vectorOfTimes, vectorOfStates );

    // Return pair of simulated ranges and reception times.
    return std::make_pair( calculatedObservation, isObservationFeasible );
}

//! Function to simulate observables, checking whether they are viable according to settings passed to this function
/*!
 *  Function to simulate observables, checking whether they are viable according to settings passed to this function
 *  NOTE: Viability check is turned of at present, will be incorporated in subsequent pull request
 *  \param observationTimes Times at which observables are to be computed
 *  \param observationModel Model used to compute observables
 *  \param linkEndAssociatedWithTime Model Reference link end for observables
 *  \return Observations at given time (concatenated in an Eigen vector) and associated times.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > >
simulateObservationsWithCheck(
        const std::vector< TimeType >& observationTimes,
        const boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const LinkEndType linkEndAssociatedWithTime )
{
    std::map< TimeType, Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > > observations;
    std::pair< Eigen::Matrix< ObservationScalarType, ObservationSize, 1 >, bool > simulatedObservation;

    for( unsigned int i = 0; i < observationTimes.size( ); i++ )
    {
        simulatedObservation = simulateObservationWithCheck< ObservationSize, ObservationScalarType, TimeType >(
                    observationTimes.at( i ), observationModel, linkEndAssociatedWithTime );

        // Check if receiving station can view transmitting station.
        if( simulatedObservation.second )
        {
            // If viable, add observable and time to vector of simulated data.
            observations[ observationTimes[ i ]  ] = simulatedObservation.first;
        }
    }

    // Return pair of simulated ranges and reception times.
    return std::make_pair( utilities::createConcatenatedEigenMatrixFromMapValues( observations ),
                           utilities::createVectorFromMapKeys( observations ) );
}

//! Function to simulate observables, checking whether they are viable according to settings passed to this function
/*!
 *  Function to simulate observables, checking whether they are viable according to settings passed to this function
 *  NOTE: Viability check is turned of at present, will be incorporated in subsequent pull request
 *  \param observationTimes Times at which observables are to be computed
 *  \param observationModel Model used to compute observables
 *  \param linkEndAssociatedWithTime Model Reference link end for observables
 *  \return Observations at given times (concatenated in an Eigen vector), with associated times and reference link end.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > >
simulateObservationsWithCheckAndLinkEndIdOutput(
        const std::vector< TimeType >& observationTimes,
        const boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
        const LinkEndType linkEndAssociatedWithTime )
{
    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > simulatedObservations =
            simulateObservationsWithCheck( observationTimes, observationModel, linkEndAssociatedWithTime );

    return std::make_pair( simulatedObservations.first, std::make_pair( simulatedObservations.second, linkEndAssociatedWithTime ) );
}


//! Objects used to simulate a set of observations of a given kind
/*!
 *  Objects used to simulate a set of observations of a given kind.
 *  NOTE: In the current application, this class does not add much to the existing interfaces. It is included in
 *  anticipation of the inclusion of observation viability checking, measurement noise simulation, etc., which
 *  will be added in an upcoming pull request
 */
template< int ObservationSize, typename ObservationScalarType = double, typename TimeType = double >
class ObservationSimulator
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Type of observable for which this object computes observations
     * \param observationModels List of observation models of type observableType
     */
    ObservationSimulator(
            const ObservableType observableType,
            const std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize,
            ObservationScalarType, TimeType > > >& observationModels ):
        observableType_( observableType ), observationModels_( observationModels ){ }

    //! Virtual destructor
    virtual ~ObservationSimulator( ){ }

    //! Function to get the type of observable for which this object computes observations
    /*!
     * Function to get the type of observable for which this object computes observations
     * \return Type of observable for which this object computes observations
     */
    ObservableType getObservableType( )
    {
        return observableType_;
    }

    //! Function to get the size of the observable for a given set of link ends
    /*!
     * Function to get the size of the observable for a given set of link ends
     * \return Size of the observable for a given set of link ends
     */
    int getObservationSize( const LinkEnds& linkEnds )
    {
        if( observationModels_.count( linkEnds ) == 0 )
        {
            throw std::runtime_error( "Error, could not find observation models for requested link ends" );
        }
        return observationModels_.at( linkEnds )->getObservationSize( );
    }

    //! Function to get the observation model for a given set of link ends
    /*!
     * Function to get the observation model for a given set of link ends
     * \return Observation model for a given set of link ends
     */
    boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > getObservationModel(
            const LinkEnds linkEnds )
    {
        if( observationModels_.count( linkEnds ) == 0 )
        {
            throw std::runtime_error(
                        "Error in observation manager when getting observation model, did not find model for given link ends " +
                        boost::lexical_cast< std::string >( observationModels_.size( ) ) );
        }
        return observationModels_.at( linkEnds );
    }

    //! Function to get the full list of observation models
    /*!
     * Function to get the full list of observation models
     * \return Full list of observation models
     */
    std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > >
    getObservationModels( )
    {
        return observationModels_;
    }

    //! Function to simulate observations between specified link ends.
    /*!
     *  Function to simulate observations between specified link ends. Users can specify whether to check for availability of
     *  link at given reception time.
     *  \param observationTimes Vector of times at which observations taked place (i.e. reception time)
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     *  \param linkEndAssociatedWithTime Reference link end for observable
     *  \param checkTimes Boolean denoting whether the observation times are to be checked for viability
     *  \return Observations at given times (concatenated in an Eigen vector), with associated times
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

        boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > selectedObservationModel =
                observationModels_.at( linkEnds );

        return simulateObservationsWithCheck< ObservationSize, ObservationScalarType, TimeType >(
                    observationTimes, selectedObservationModel, linkEndAssociatedWithTime );
    }

    //! Function to simulate observations between specified link ends.
    /*!
     *  Function to simulate observations between specified link ends. Users can specify whether to check for availability of
     *  link at given reception time.
     *  \param observationTimes Vector of times at which observations taked place (i.e. reception time)
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     *  \param linkEndAssociatedWithTime Reference link end for observable
     *  \param checkTimes Boolean denoting whether the observation times are to be checked for viability
     *  \return Observations at given times (concatenated in an Eigen vector), with associated times and reference link end.
     */
    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > >
    simulateObservationsWithLinkEndId( const std::vector< TimeType >& observationTimes,
                                       const LinkEnds linkEnds,
                                       const LinkEndType linkEndAssociatedWithTime,
                                       const bool checkTimes = true )
    {
        std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > simulatedObservations =
                simulateObservations( observationTimes, linkEnds, linkEndAssociatedWithTime, checkTimes );

        return std::make_pair( simulatedObservations.first, std::make_pair(
                                   simulatedObservations.second, linkEndAssociatedWithTime ) );
    }

protected:

    //! Type of observable for which this object computes observations
    ObservableType observableType_;

    //! List of observation models of type observableType
    std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > >
    observationModels_;
};

}

}
#endif // TUDAT_OBSERVATIONSIMULATOR_H
