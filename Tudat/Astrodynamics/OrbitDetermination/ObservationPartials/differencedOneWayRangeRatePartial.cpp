/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/differencedOneWayRangeRatePartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to calculate the observation partial(s) at required time and state
std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > >
DifferencedOneWayRangeRatePartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
{

    using namespace observation_partials;

    // Split input times/states for arc start and end ranges
    std::vector< Eigen::Vector6d > arcStartStates;
    arcStartStates.push_back( states.at( 0 ) );
    arcStartStates.push_back( states.at( 1 ) );
    std::vector< Eigen::Vector6d > arcEndStates;
    arcEndStates.push_back( states.at( 2 ) );
    arcEndStates.push_back( states.at( 3 ) );

    std::vector< double > arcStartTimes;
    arcStartTimes.push_back( times.at( 0 ) );
    arcStartTimes.push_back( times.at( 1 ) );
    std::vector< double > arcEndTimes;
    arcEndTimes.push_back( times.at( 2 ) );
    arcEndTimes.push_back( times.at( 3 ) );


    // Obtain arc start and end range partials
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > arcStartPartials =
            arcStartRangePartial_->calculatePartial( arcStartStates, arcStartTimes, linkEndOfFixedTime );

    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > arcEndPartials =
            arcEndRangePartial_->calculatePartial( arcEndStates, arcEndTimes, linkEndOfFixedTime );

    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > differencedRangeRatePartials;

    if( arcStartPartials.size( ) != arcEndPartials.size( ) )
    {
        std::cerr << "Error when making differenced one way range rate partials, arc start and end partials inconsistent" << std::endl;
    }

    // Retrieve arc length
    double arcDuration = TUDAT_NAN;
    if ( linkEndOfFixedTime == observation_models::transmitter )
    {
        arcDuration = times[ 2 ] - times[ 0 ];
    }
    else if ( linkEndOfFixedTime == observation_models::receiver )
    {
        arcDuration = times[ 3 ] - times[ 1 ];
    }


    // Scale partials by arc duration
    for( unsigned int i = 0; i < arcStartPartials.size( ); i++ )
    {
        differencedRangeRatePartials.push_back(
                    std::make_pair( -arcStartPartials[ i ].first / arcDuration, arcStartPartials[ i ].second ) );
    }

    for( unsigned int i = 0; i < arcEndPartials.size( ); i++ )
    {
        differencedRangeRatePartials.push_back(
                    std::make_pair( arcEndPartials[ i ].first / arcDuration, arcEndPartials[ i ].second ) );
    }

    return differencedRangeRatePartials;
}

}

}

