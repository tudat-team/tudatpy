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

template< typename ObservationScalarType = double,
          typename TimeType = double >
class NWayDifferencedRangeObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:    
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;


    NWayDifferencedRangeObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > arcStartObservationModel,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > arcEndObservationModel,
            const std::function< double( const double ) > integrationTimeFunction,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( n_way_differenced_range, linkEnds, observationBiasCalculator ),
        arcStartObservationModel_( arcStartObservationModel ),
        arcEndObservationModel_( arcEndObservationModel ),
        integrationTimeFunction_( integrationTimeFunction )
    {
        numberOfLinkEnds_ = linkEnds.size( );
    }

    //! Destructor
    ~NWayDifferencedRangeObservationModel( ){ }

    //! Function to compute one-way range rate observation at given time.
    /*!
     *  This function computes the one-way observation at a given time.
     *  \param time Time at which observation is to be simulated.
     *  \param isTimeAtReception True if given time is to be the reception time, false if it is transmission time.
     *  \return Calculated observed one-way range rate value.
     */
    ObservationScalarType computeObservation(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime ) const
    {
        ObservationScalarType countIntervalDuration = countIntervalDurationFunction_( );
        return ( arcEndObservationModel_->computeObservation( time, linkEndAssociatedWithTime ) -
                arcStartObservationModel_->computeObservation( time - countIntervalDuration, linkEndAssociatedWithTime ) ) /
                static_cast< ObservationScalarType >( countIntervalDuration );
    }


    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< TimeType >& linkEndTimes,
            std::vector< Eigen::Matrix< ObservationScalarType, 6, 1 > >& linkEndStates )  const
    {
        std::vector< TimeType > arcStartLinkEndTimes;
        std::vector< Eigen::Matrix< ObservationScalarType, 6, 1 > > arcStartLinkEndStates;
        std::vector< TimeType > arcEndLinkEndTimes;
        std::vector< Eigen::Matrix< ObservationScalarType, 6, 1 > > arcEndLinkEndStates;

        ObservationScalarType countIntervalDuration = countIntervalDurationFunction_( );

        ObservationScalarType observation =
                ( arcEndObservationModel_->computeObservationAndFullPrecisionLinkEndData(
                    time, linkEndAssociatedWithTime, arcEndLinkEndTimes, arcEndLinkEndStates ) -
                arcStartObservationModel_->computeObservationAndFullPrecisionLinkEndData(
                    time - countIntervalDuration, linkEndAssociatedWithTime, arcStartLinkEndTimes, arcStartLinkEndStates ) ) /
                static_cast< ObservationScalarType >( countIntervalDuration );

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

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << observation ).finished( );
    }

//    void resetRetransmissionDelaysFunctions(
//            const boost::function< std::vector< ObservationScalarType >( ) > arcStartRetransmissionDelaysFunction,
//            const boost::function< std::vector< ObservationScalarType >( ) > arcEndRetransmissionDelaysFunction )
//    {
//        arcStartObservationModel_->resetRetransmissionDelaysFunction( arcStartRetransmissionDelaysFunction );
//        arcEndObservationModel_->resetRetransmissionDelaysFunction( arcEndRetransmissionDelaysFunction );
//    }

    boost::function< ObservationScalarType( ) > getIntegrationTimeFunction( )
    {
        return integrationTimeFunction_;
    }

    boost::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType, ObservationScalarType > > getArcEndObservationModel( )
    {
        return arcEndObservationModel_;
    }

    boost::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType, ObservationScalarType > > getArcStartObservationModel( )
    {
        return arcStartObservationModel_;
    }
private:

    boost::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel_;

    boost::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel_;

    std::function< double( const double ) > integrationTimeFunction_;

    int numberOfLinkEnds_;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_NWAYRANGERATEOBSERVATIONMODEL_H
