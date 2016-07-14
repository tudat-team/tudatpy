/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Propulsion/thrustGuidance.h"

namespace tudat
{

namespace propulsion
{

Eigen::Vector3d getThrustDirectionColinearWithVelocity(
        const basic_mathematics::Vector6d& currentState, const double currentTime, const bool putThrustInOppositeDirection )
{
    return ( ( putThrustInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 3, 3 ) ).normalized( );
}

Eigen::Vector3d getThrustDirectionColinearWithPosition(
        const basic_mathematics::Vector6d& currentState, const double currentTime, const bool putThrustInOppositeDirection )
{
    return ( ( putThrustInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 0, 3 ) ).normalized( );
}

Eigen::Vector3d getThrustDirectionFromTimeOnlyFunction(
        const basic_mathematics::Vector6d& currentState, const double currentTime,
                 const boost::function< Eigen::Vector3d( const double ) > timeOnlyFunction )
{
    return timeOnlyFunction( currentTime ).normalized( );
}

} // namespace propulsion

} // namespace tudat

