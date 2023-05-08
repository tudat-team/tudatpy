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
        const simulation_setup::SystemOfBodies&,
        const LinkEnds& linkEnds,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings,
        const bool isFirstPartial )
{
    double integrationTime;
    try
    {
        integrationTime = ancillarySettings->getAncilliaryDoubleData( doppler_integration_time, true );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error when retrieving integration time for one-way averaged Doppler observable: " +
                        std::string( caughtException.what( ) ) );
    }
    return 1.0 / integrationTime;
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
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetings = nullptr ) const
    {
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;

        return computeIdealObservationsWithLinkEndData(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancilliarySetings );
    }


    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetings = nullptr )
    {
        std::vector< double > arcStartLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcStartLinkEndStates;
        std::vector< double > arcEndLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcEndLinkEndStates;

        TimeType integrationTime;
        try
        {
            integrationTime = ancilliarySetings->getAncilliaryDoubleData( doppler_integration_time, true );
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
