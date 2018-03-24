/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <boost/make_shared.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterGooding.h"

//! Tudat library namespace.
namespace tudat
{
namespace mission_segments
{

using namespace root_finders;

//! Constructor with immediate definition of parameters and execution of the algorithm.
LambertTargeterGooding::LambertTargeterGooding( 
        const Eigen::Vector3d& aCartesianPositionAtDeparture,
        const Eigen::Vector3d& aCartesianPositionAtArrival,
        const double aTimeOfFlight,
        const double aGravitationalParameter,
        RootFinderPointer aRootFinder )
    : LambertTargeter( aCartesianPositionAtDeparture, aCartesianPositionAtArrival,
                       aTimeOfFlight, aGravitationalParameter ),
      rootFinder( aRootFinder )
{
    // Required because the make_shared in the function definition gives problems for MSVC.
    if ( !rootFinder.get( ) )
    {
        rootFinder = boost::make_shared< NewtonRaphson >( 1.0e-12, 1000 );
    }

    // Execute algorithm.
    execute( );
}

//! Execute Lambert targeting solver.
void LambertTargeterGooding::execute( )
{
    // Call Gooding's Lambert targeting routine.
    solveLambertProblemGooding( cartesianPositionAtDeparture, cartesianPositionAtArrival,
                                timeOfFlight, gravitationalParameter,
                                cartesianVelocityAtDeparture, cartesianVelocityAtArrival,
                                rootFinder );
}

//! Get radial velocity at departure.
double LambertTargeterGooding::getRadialVelocityAtDeparture( )
{
    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtDeparture = cartesianPositionAtDeparture.normalized( );

    // Compute radial velocity at departure.
    return cartesianVelocityAtDeparture.dot( radialUnitVectorAtDeparture );
}

//! Get radial velocity at arrival.
double LambertTargeterGooding::getRadialVelocityAtArrival( )
{
    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival.normalized( );

    // Compute radial velocity at arrival.
    return cartesianVelocityAtArrival.dot( radialUnitVectorAtArrival );
}

//! Get transverse velocity at departure.
double LambertTargeterGooding::getTransverseVelocityAtDeparture( )
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
double LambertTargeterGooding::getTransverseVelocityAtArrival( )
{
    // Compute angular momemtum vector.
    const Eigen::Vector3d angularMomentumVector =
            cartesianPositionAtArrival.cross( cartesianVelocityAtArrival );

    // Compute normalized angular momentum vector.
    const Eigen::Vector3d angularMomentumUnitVector = angularMomentumVector.normalized( );

    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival.normalized( );

    // Compute tangential unit vector.
    Eigen::Vector3d tangentialUnitVectorAtArrival =
                angularMomentumUnitVector.cross( radialUnitVectorAtArrival );

    // Compute tangential velocity at departure.
    return cartesianVelocityAtArrival.dot( tangentialUnitVectorAtArrival );
}

//! Get semi-major axis.
double LambertTargeterGooding::getSemiMajorAxis( )
{
    // Compute specific orbital energy: eps = v^2/ - mu/r.
    const double specificOrbitalEnergy = cartesianVelocityAtDeparture.squaredNorm( ) / 2.0
            - gravitationalParameter / cartesianPositionAtDeparture.norm( );

    // Compute semi-major axis: a = -mu / 2*eps.
    return -gravitationalParameter / ( 2.0 * specificOrbitalEnergy );
}

} // namespace mission_segments
} // namespace tudat
