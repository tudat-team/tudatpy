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

#include <iostream>

#include "tudat/astro/electromagnetism/yarkovskyAcceleration.h"

namespace tudat
{
namespace electromagnetism
{

//! Compute Yarkovsky acceleration using the transverse direction of the RTN frame.
Eigen::Vector3d computeYarkovskyAcceleration( double yarkovskyParameter, const Eigen::Vector6d& stateVector )
{
    const Eigen::Vector3d currentPosition = stateVector.segment( 0, 3 );
    const Eigen::Vector3d currentVelocity = stateVector.segment( 3, 3 );

    const double minAngle = 1.0e-12;
    if( currentPosition.norm( ) == 0.0 || currentVelocity.norm( ) == 0.0 ||
        std::abs( currentPosition.dot( currentVelocity ) ) / ( currentPosition.norm( ) * currentVelocity.norm( ) ) > std::cos( minAngle ) )
    {
        std::cerr << "Warning when computing yarkovsky acceleration, position and velocity are close to parallel, or one has norm zero; "
                     "returning zero acceleration"
                  << std::endl;
        return Eigen::Vector3d::Zero( );
    }

    const double currentDistance = currentPosition.norm( );
    const double auOverDistance = physical_constants::ASTRONOMICAL_UNIT / currentDistance;
    const double yarkovskyMagnitude = yarkovskyParameter * auOverDistance * auOverDistance;

    const Eigen::Vector3d radialUnitVector = currentPosition / currentDistance;
    const Eigen::Vector3d transverseVelocity = currentVelocity - currentVelocity.dot( radialUnitVector ) * radialUnitVector;
    const Eigen::Vector3d yarkovskyDirection = transverseVelocity.normalized( );

    return yarkovskyMagnitude * yarkovskyDirection;
}

}  // namespace electromagnetism
}  // namespace tudat
