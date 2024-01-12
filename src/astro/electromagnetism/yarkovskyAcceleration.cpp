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


namespace tudat
{
namespace electromagnetism
{

//! Compute Yarkovsky Acceleration using a simplified tangential model.
Eigen::Vector3d computeYarkovskyAcceleration( double yarkovskyParameter, const Eigen::Vector6d& stateVector )
{
    const double currentDistance = stateVector.segment( 0, 3 ).norm( );
    const double auOverDistance = physical_constants::ASTRONOMICAL_UNIT / currentDistance;
    const double yarkovskyMagnitude = yarkovskyParameter * auOverDistance * auOverDistance;

    // Find the direction (tangential to the orbit and thus parallel to the velocity)
    const Eigen::Vector3d yarkovskyDirection = stateVector.segment( 3, 3 ).normalized( );

    return yarkovskyMagnitude * yarkovskyDirection;
}

} // namespace electromagnetism
} // namespace tudat
