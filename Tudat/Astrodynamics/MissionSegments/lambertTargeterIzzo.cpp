/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Battin, R.H. An Introduction to the Mathematics and Methods of Astrodynamics,
 *          AIAA Education Series, 1999.
 *      Izzo, D. lambert_problem.h, keptoolbox.
 *
 *    Notes
 *      This code is an implementation of the method developed by Dario Izzo from ESA/ACT and
 *      publicly available at: http://keptoolbox.sourceforge.net/.
 *      After verification and validation, it was proven that this algorithm is faster and more
 *      robust than the implemented Lancaster & Blanchard and Gooding method. Notably, this method
 *      does not suffer from the near-pi singularity (pi-transfers are by nature singular).
 *
 */

#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterIzzo.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"

#include <Eigen/Geometry>

namespace tudat
{
namespace mission_segments
{

//! Execute Lambert targeting solver.
void LambertTargeterIzzo::execute( )
{
    // Call Izzo's Lambert targeting routine.
    solveLambertProblemIzzo( cartesianPositionAtDeparture, cartesianPositionAtArrival,
                             timeOfFlight, gravitationalParameter, cartesianVelocityAtDeparture,
                             cartesianVelocityAtArrival, isRetrograde_, convergenceTolerance_,
                             maximumNumberOfIterations_ );
}

//! Get radial velocity at departure.
double LambertTargeterIzzo::getRadialVelocityAtDeparture( )
{
    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtDeparture
            = cartesianPositionAtDeparture.normalized( );

    // Compute radial velocity at departure.
    return cartesianVelocityAtDeparture.dot( radialUnitVectorAtDeparture );
}

//! Get radial velocity at arrival.
double LambertTargeterIzzo::getRadialVelocityAtArrival( )
{
    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival.normalized( );

    // Compute radial velocity at arrival.
    return cartesianVelocityAtArrival.dot( radialUnitVectorAtArrival );
}

//! Get transverse velocity at departure.
double LambertTargeterIzzo::getTransverseVelocityAtDeparture( )
{
    // Compute angular momemtum vector.
    const Eigen::Vector3d angularMomentumVector =
            cartesianPositionAtDeparture.cross( cartesianVelocityAtDeparture );

    // Compute normalized angular momentum vector.
    const Eigen::Vector3d angularMomentumUnitVector = angularMomentumVector.normalized( );

    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtDeparture
            = cartesianPositionAtDeparture.normalized( );

    // Compute tangential unit vector.
    Eigen::Vector3d tangentialUnitVectorAtDeparture =
                angularMomentumUnitVector.cross( radialUnitVectorAtDeparture );

    // Compute tangential velocity at departure.
    return cartesianVelocityAtDeparture.dot( tangentialUnitVectorAtDeparture );
}

//! Get transverse velocity at arrival.
double LambertTargeterIzzo::getTransverseVelocityAtArrival( )
{
    // Compute angular momemtum vector.
    const Eigen::Vector3d angularMomentumVector =
            cartesianPositionAtArrival.cross( cartesianVelocityAtArrival );

    // Compute normalized angular momentum vector.
    const Eigen::Vector3d angularMomentumUnitVector = angularMomentumVector.normalized( );

    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival.normalized( );

    // Compute tangential unit vector.
    Eigen::Vector3d tangentialUnitVectorAtArrival
            = angularMomentumUnitVector.cross( radialUnitVectorAtArrival );

    // Compute tangential velocity at departure.
    return cartesianVelocityAtArrival.dot( tangentialUnitVectorAtArrival );
}

//! Get semi-major axis.
double LambertTargeterIzzo::getSemiMajorAxis( )
{
    // Compute specific orbital energy: eps = v^2/ - mu/r.
    const double specificOrbitalEnergy = cartesianVelocityAtDeparture.squaredNorm( ) / 2.0
            - gravitationalParameter / cartesianPositionAtDeparture.norm( );

    // Compute semi-major axis: a = -mu / 2*eps.
    return -gravitationalParameter / ( 2.0 * specificOrbitalEnergy );
}

} // namespace mission_segments
} // namespace tudat
