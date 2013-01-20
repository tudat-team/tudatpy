/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120210    T. Secretin       First creation of code.
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

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterIzzo.h"

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
