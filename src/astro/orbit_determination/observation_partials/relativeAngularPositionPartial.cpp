/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/relativeAngularPositionPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/angularPositionPartial.h"

namespace tudat
{

namespace observation_partials
{


//! Update the scaling object to the current times and states
void RelativeAngularPositionScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                             const std::vector< double >& times,
                                             const observation_models::LinkEndType fixedLinkEnd,
                                             const Eigen::VectorXd currentObservation )
{
    if ( fixedLinkEnd != observation_models::receiver )
    {
        throw std::runtime_error( "Error when updating a relative angular position scaling object, fixed link end time different from receiver." );
    }

//    Eigen::Vector6d firstTransmitterState = linkEndStates[ 0 ];
//    Eigen::Vector6d secondTransmitterState = linkEndStates[ 1 ];
//    Eigen::Vector6d receiverState = linkEndStates[ 2 ];

    Eigen::Vector3d relativeRangeVectorFirstTransmitter = ( linkEndStates[ 2 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d relativeRangeVectorSecondTransmitter = ( linkEndStates[ 2 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

    Eigen::Vector3d normalizedRelativeRangeVectorFirstTransmitter = relativeRangeVectorFirstTransmitter.normalized( );
    Eigen::Vector3d normalizedRelativeRangeVectorSecondTransmitter = relativeRangeVectorSecondTransmitter.normalized( );

    // Compute common scaling factor
    scalingFactorFirstTransmitter_ = calculatePartialOfAngularPositionWrtLinkEndPosition( relativeRangeVectorFirstTransmitter, true );
    scalingFactorSecondTransmitter_ = calculatePartialOfAngularPositionWrtLinkEndPosition( relativeRangeVectorSecondTransmitter, true );

    // Compute scaling for receiver reference
//    if( fixedLinkEnd == observation_models::receiver )
//    {
    referenceLightTimeCorrectionScalingFirstTransmitter_ = scalingFactorFirstTransmitter_ * linkEndStates[ 0 ].segment( 3, 3 ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVectorFirstTransmitter ) );
    referenceScalingFactorFirstTransmitter_ =
            scalingFactorFirstTransmitter_ *
            ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 0 ].segment( 3, 3 ) * normalizedRelativeRangeVectorFirstTransmitter.transpose( ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVectorFirstTransmitter ) ) );

    referenceLightTimeCorrectionScalingSecondTransmitter_ = scalingFactorSecondTransmitter_ * linkEndStates[ 1 ].segment( 3, 3 ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVectorSecondTransmitter ) );
    referenceScalingFactorSecondTransmitter_ = scalingFactorSecondTransmitter_ *
            ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 1 ].segment( 3, 3 ) * normalizedRelativeRangeVectorSecondTransmitter.transpose( ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVectorSecondTransmitter ) ) );
//    }


    currentLinkEndType_ = fixedLinkEnd;

}


}

}
