/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      101111    E. Iorfida        Creation of code.
 *      101111    E. Iorfida        Implementation of all the equations up to the Newton method.
 *      101117    E. Iorfida        Velocities computations added.
 *      101126    E. Iorfida        Get/set codes deleted.
 *      101206    E. Iorfida        LambertTargetingElements class deleted,
 *                                  added setInitialState, modified punctuation. Set single
 *                                  variables, change variables names in more understandable ones.
 *      101209    E. Iorfida        Corrected some coding errors.
 *      101213    E. Iorfida        Deleted lambertAngle, added numberOfRevolution, modified
 *                                  implementation.
 *      101214    E. Iorfida        Implementation only for the case with numberOfRevolution = 0.
 *      110113    E. Iorfida        Added necessary elements to build pointer-to-member-function
 *                                  to RootFinderAlgorithms and NewtonRaphsonMethod classes.
 *      110124    E. Iorfida        Added necessary piece of code to be able to use the last
 *                                  version of Newton-Raphson code.
 *      110126    E. Iorfida        Initialized member functions.
 *      110130    J. Melman         Simplified variable names, e.g., 'normOfdeletVelocityVector'
 *                                  became 'speed'. Requested references to specific formulas. Also
 *                                  corrected 'tangential' to 'transverse'. Simplified computation
 *                                  of radial unit vector. Corrected computation of transverse
 *                                  heliocentric velocity.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified variable names (from
 *                                  heliocentric, to inertial). Added patch for negative case of
 *                                  initialLambertGuess_. Added equations references.
 *      110206    E. Iorfida        Added unique function for Newton-Raphson method. Added
 *                                  computeAbsoluteValue to the initialLambertGuess_ for
 *                                  non-converging cases.
 *      110208    E. Iorfida        Added CartesianPositionElements objects as input and
 *                                  CartesianVelocityElements objects as output.
 *      110418    E. Iorfida        Added a new normal plane that take into account the case of two
 *                                  parallel position vector (with a relative angle of 180
 *                                  degrees). Better defined the pointers to the output
 *                                  CartesianVelocityElements.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120620    T. Secretin       Adapted and moved code from LambertTargeter.cpp.
 *
 *    References
 *
 *    Notes
 *
 */

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterGooding.h"

//! Tudat library namespace.
namespace tudat
{
namespace mission_segments
{

//! Execute Lambert targeting solver.
void LambertTargeterGooding::execute( )
{
    // Call Gooding's Lambert targeting routine.
    solveLambertProblemGooding( cartesianPositionAtDeparture_, cartesianPositionAtArrival_,
                                timeOfFlight_, gravitationalParameter_,
                                cartesianVelocityAtDeparture_, cartesianVelocityAtArrival_,
                                newtonRaphson_, convergenceTolerance_,
                                maximumNumberOfIterations_ );
}

//! Get radial velocity at departure.
double LambertTargeterGooding::getRadialVelocityAtDeparture( )
{
    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtDeparture = cartesianPositionAtDeparture_.normalized( );

    // Compute radial velocity at departure.
    return cartesianVelocityAtDeparture_.dot( radialUnitVectorAtDeparture );
}

//! Get radial velocity at arrival.
double LambertTargeterGooding::getRadialVelocityAtArrival( )
{
    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival_.normalized( );

    // Compute radial velocity at arrival.
    return cartesianVelocityAtArrival_.dot( radialUnitVectorAtArrival );
}

//! Get transverse velocity at departure.
double LambertTargeterGooding::getTransverseVelocityAtDeparture( )
{
    // Compute angular momemtum vector.
    const Eigen::Vector3d angularMomentumVector =
            cartesianPositionAtDeparture_.cross( cartesianVelocityAtDeparture_ );

    // Compute normalized angular momentum vector.
    const Eigen::Vector3d angularMomentumUnitVector = angularMomentumVector.normalized( );

    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtDeparture
            = cartesianPositionAtDeparture_.normalized( );

    // Compute tangential unit vector.
    Eigen::Vector3d tangentialUnitVectorAtDeparture =
                angularMomentumUnitVector.cross( radialUnitVectorAtDeparture );

    // Compute tangential velocity at departure.
    return cartesianVelocityAtDeparture_.dot( tangentialUnitVectorAtDeparture );
}

//! Get transverse velocity at arrival.
double LambertTargeterGooding::getTransverseVelocityAtArrival( )
{
    // Compute angular momemtum vector.
    const Eigen::Vector3d angularMomentumVector =
            cartesianPositionAtArrival_.cross( cartesianVelocityAtArrival_ );

    // Compute normalized angular momentum vector.
    const Eigen::Vector3d angularMomentumUnitVector = angularMomentumVector.normalized( );

    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival_.normalized( );

    // Compute tangential unit vector.
    Eigen::Vector3d tangentialUnitVectorAtArrival =
                angularMomentumUnitVector.cross( radialUnitVectorAtArrival );

    // Compute tangential velocity at departure.
    return cartesianVelocityAtArrival_.dot( tangentialUnitVectorAtArrival );
}

//! Get semi-major axis.
double LambertTargeterGooding::getSemiMajorAxis( )
{
    // Compute specific orbital energy: eps = v^2/ - mu/r.
    const double specificOrbitalEnergy = cartesianVelocityAtDeparture_.squaredNorm( ) / 2.0
            - gravitationalParameter_ / cartesianPositionAtDeparture_.norm( );

    // Compute semi-major axis: a = -mu / 2*eps.
    return -gravitationalParameter_ / ( 2.0 * specificOrbitalEnergy );
}

} // namespace mission_segments
} // namespace tudat
