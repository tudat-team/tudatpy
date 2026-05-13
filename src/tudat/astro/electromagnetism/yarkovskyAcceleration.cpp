/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/electromagnetism/yarkovskyAcceleration.h"

#include <cmath>
#include <limits>

namespace tudat
{
namespace electromagnetism
{

//! Compute Yarkovsky acceleration using the transverse direction of the RTN frame.
Eigen::Vector3d computeYarkovskyAcceleration( double yarkovskyParameter, const Eigen::Vector6d& stateVector )
{
    const Eigen::Vector3d currentPosition = stateVector.segment( 0, 3 );
    const Eigen::Vector3d currentVelocity = stateVector.segment( 3, 3 );

    const double currentDistance = currentPosition.norm( );
    if( currentDistance <= std::numeric_limits< double >::epsilon( ) )
    {
        return Eigen::Vector3d::Zero( );
    }

    const double auOverDistance = physical_constants::ASTRONOMICAL_UNIT / currentDistance;
    const double yarkovskyMagnitude = yarkovskyParameter * auOverDistance * auOverDistance;

    const Eigen::Vector3d radialUnitVector = currentPosition / currentDistance;
    const Eigen::Vector3d transverseVelocity = currentVelocity - currentVelocity.dot( radialUnitVector ) * radialUnitVector;
    const double transverseSpeed = transverseVelocity.norm( );
    const double transverseSpeedTolerance =
            std::sqrt( std::numeric_limits< double >::epsilon( ) ) * currentVelocity.norm( );
    if( transverseSpeed <= transverseSpeedTolerance )
    {
        return Eigen::Vector3d::Zero( );
    }

    const Eigen::Vector3d yarkovskyDirection = transverseVelocity / transverseSpeed;

    return yarkovskyMagnitude * yarkovskyDirection;
}

}  // namespace electromagnetism
}  // namespace tudat
