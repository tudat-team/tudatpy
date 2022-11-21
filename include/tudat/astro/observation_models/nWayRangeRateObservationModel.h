/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NWAYRANGERATEOBSERVATIONMODEL_H
#define TUDAT_NWAYRANGERATEOBSERVATIONMODEL_H

#include "tudat/astro/observation_models/nWayRangeObservationModel.h"

namespace tudat
{

namespace observation_models
{

inline double getDifferencedNWayRangeScalingFactor(
        const std::vector< double >& linkEndTimes,
        const observation_models::LinkEndType referenceLinkEnd )
{
    int numberOfEntries = linkEndTimes.size( ) / 2;
    double arcDuration = TUDAT_NAN;
    if ( referenceLinkEnd == observation_models::transmitter )
    {
        arcDuration = linkEndTimes[ numberOfEntries ] - linkEndTimes[ 0 ];
    }
    else if ( referenceLinkEnd == observation_models::receiver )
    {
        arcDuration = linkEndTimes[ 2 * numberOfEntries - 1 ] - linkEndTimes[ numberOfEntries - 1 ];
    }
    else
    {
        throw std::runtime_error( "Error when getting differenced n-way range scaling factor; link end " +
                                  getLinkEndTypeString( referenceLinkEnd ) + " not recognized." );
    }
    return 1.0 / arcDuration;
}


template< typename ObservationScalarType = double,
          typename TimeType = double >
class NWayDifferencedRangeObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:    
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;


    NWayDifferencedRangeObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( n_way_differenced_range, linkEnds, observationBiasCalculator ),
        arcStartObservationModel_( arcStartObservationModel ),
        arcEndObservationModel_( arcEndObservationModel ),
        numberOfLinkEnds_( linkEnds.size( ) ){ }

    //! Destructor
    ~NWayDifferencedRangeObservationModel( ){ }

    //! Function to compute one-way range rate observation at given time.
    /*!
     *  This function computes the one-way observation at a given time.
     *  \param time Time at which observation is to be simulated.
     *  \param isTimeAtReception True if given time is to be the reception time, false if it is transmission time.
     *  \return Calculated observed one-way range rate value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            const std::shared_ptr< ObservationAncilliarySimulationSettings< TimeType > > ancilliarySetings = nullptr ) const
    {
        if( ancilliarySetings == nullptr )
        {
            throw std::runtime_error( "Error when simulating one-way averaged Doppler observable; no ancilliary settings found. Ancilliary settings are requiured for integration time" );
        }

        TimeType integrationTime;
        try
        {
            integrationTime = ancilliarySetings->getAncilliaryDoubleData( averaged_doppler_integration_time, true );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error( "Error when retrieving integration time for one-way averaged Doppler observable: " +
                            std::string( caughtException.what( ) ) );
        }

        return ( arcEndObservationModel_->computeObservations( time, linkEndAssociatedWithTime ) -
                arcStartObservationModel_->computeObservations( time - integrationTime, linkEndAssociatedWithTime ) ) /
                static_cast< ObservationScalarType >( integrationTime );
    }


    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings< TimeType > > ancilliarySetings = nullptr )
    {
        std::vector< double > arcStartLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcStartLinkEndStates;
        std::vector< double > arcEndLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcEndLinkEndStates;

        TimeType integrationTime;
        try
        {
            integrationTime = ancilliarySetings->getAncilliaryDoubleData( averaged_doppler_integration_time, true );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error( "Error when retrieving integration time for one-way averaged Doppler observable: " +
                            std::string( caughtException.what( ) ) );
        }

        Eigen::Matrix< ObservationScalarType, 1, 1 > observation =
                ( arcEndObservationModel_->computeIdealObservationsWithLinkEndData(
                    time + integrationTime / 2.0, linkEndAssociatedWithTime, arcEndLinkEndTimes, arcEndLinkEndStates, ancilliarySetings ) -
                arcStartObservationModel_->computeIdealObservationsWithLinkEndData(
                    time - integrationTime / 2.0, linkEndAssociatedWithTime, arcStartLinkEndTimes, arcStartLinkEndStates, ancilliarySetings ) ) /
                static_cast< ObservationScalarType >( integrationTime );

        linkEndTimes.clear( );
        linkEndStates.clear( );
        linkEndTimes.resize( 4 * ( numberOfLinkEnds_ - 1 ) );
        linkEndStates.resize( 4 * ( numberOfLinkEnds_ - 1 ) );

        for( int i = 0; i < 2 * ( numberOfLinkEnds_ - 1 ) ; i++ )
        {
            linkEndTimes[ i ] = arcStartLinkEndTimes[ i ];
            linkEndTimes[ i + 2 * ( numberOfLinkEnds_ - 1 ) ] = arcEndLinkEndTimes[ i ];

            linkEndStates[ i ] = arcStartLinkEndStates[ i ];
            linkEndStates[ i + 2 * ( numberOfLinkEnds_ - 1 ) ] = arcEndLinkEndStates[ i ];
        }

        return observation;
    }

    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > getArcEndObservationModel( )
    {
        return arcEndObservationModel_;
    }

    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > getArcStartObservationModel( )
    {
        return arcStartObservationModel_;
    }
private:

    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel_;

    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel_;

    int numberOfLinkEnds_;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_NWAYRANGERATEOBSERVATIONMODEL_H
