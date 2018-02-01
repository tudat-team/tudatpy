/*    Copyright (c) 2010-2018, Delft University of Technology
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

//! Function to get the unit vector colinear with velocity segment of a translational state.
Eigen::Vector3d getForceDirectionColinearWithVelocity(
        const boost::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putForceInOppositeDirection )
{
    static Eigen::Vector6d currentState;
    currentStateFunction( currentState );
    return ( ( putForceInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 3, 3 ) ).normalized( );
}

//! Function to get the unit vector colinear with position segment of a translational state.
Eigen::Vector3d getForceDirectionColinearWithPosition(
        const boost::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putForceInOppositeDirection )
{
    static Eigen::Vector6d currentState;
    currentStateFunction( currentState );
    return ( ( putForceInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 0, 3 ) ).normalized( );
}

//! Function to get the force direction from a time-only function.
Eigen::Vector3d getForceDirectionFromTimeOnlyFunction(
        const double currentTime,
        const boost::function< Eigen::Vector3d( const double ) > timeOnlyFunction )
{
    return timeOnlyFunction( currentTime ).normalized( );
}

} // namespace propulsion

} // namespace tudat

