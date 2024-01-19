/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/oneWayRangePartial.h"

namespace tudat
{

namespace observation_partials
{


//! Update the scaling object to the current times and states
void OneWayRangeScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                 const std::vector< double >& times,
                                 const observation_models::LinkEndType fixedLinkEnd,
                                 const Eigen::VectorXd currentObservation )
{
    // Compute Euclidean distance vector
    Eigen::Vector3d rangeVector = linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 0 ].segment( 0, 3 );
    fixedLinkEndScalingFactor_ = rangeVector.transpose( ) / rangeVector.norm( );

    // Compute scaling for receiver reference
    if( fixedLinkEnd == observation_models::receiver )
    {
        referenceLightTimeCorrectionScaling_ = 1.0 / ( 1.0 - fixedLinkEndScalingFactor_.transpose( ).dot( linkEndStates[ 0 ].segment( 3, 3 ) ) /
                physical_constants::SPEED_OF_LIGHT );
        referenceScalingFactor_ =  fixedLinkEndScalingFactor_ * referenceLightTimeCorrectionScaling_;
    }

    // Compute scaling for transmitter reference
    else if( fixedLinkEnd == observation_models::transmitter )
    {
        referenceLightTimeCorrectionScaling_ =
                1.0 / ( 1.0 - fixedLinkEndScalingFactor_.transpose( ).dot( linkEndStates[ 1 ].segment( 3, 3 ) ) /
                physical_constants::SPEED_OF_LIGHT );
        referenceScalingFactor_ =  fixedLinkEndScalingFactor_ * referenceLightTimeCorrectionScaling_;
    }


    currentLinkEndType_ = fixedLinkEnd;
}


}

}
