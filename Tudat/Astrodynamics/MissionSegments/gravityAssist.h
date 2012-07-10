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
 *      110117    E. Iorfida        Creation of code.
 *      110128    E. Iorfida        Added member variables necessary for the unit tests.
 *      110202    J. Melman         Changed several names and suggested minor restructuring of
 *                                  class.
 *      110203    E. Iorfida        Changed some variable names and modified punctuation.
 *      110205    J. Melman         Removed the trailing underscores in some public variables.
 *                                  Changed and added some notes.
 *      110208    E. Iorfida        Added CartesianVelocityElements objects as input. Deleted
 *                                  inheritance from TrajectoryDesignMethod.
 *      110212    J. Melman         Made delta-V private. getDeltaV changed into computeDeltaV.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius, replaced by an element
 *                                  of GeometricShapes.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120417    T. Secretin       Moved set functions to constructor.
 *      120508    P. Musegaas       The gravitational parameter is now passed as double, made some
 *                                  variables constant.
 *      120530    P. Musegaas       Complete revision. Removed class structure, made it a free
 *                                  function. Added two functions for propagating a gravity assist
 *                                  (powered and unpowered).
 *      120625    P. Musegaas       Minor changes.
 *      120703    T. Secretin       Minor layout changes. Changed constructor.
 *
 *    Notes
 *      Note that the exact implementation of Newton-Raphson as root finder should be updated if
 *      someone would want to use a different root-finding technique.
 *
 */

#ifndef TUDAT_GRAVITY_ASSIST_H
#define TUDAT_GRAVITY_ASSIST_H

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"

namespace tudat
{
namespace mission_segments
{

//! Calculate deltaV of a gravity assist.
/*!
 * Calculates the deltaV required to perform a certain gravity assist. This function essentially
 * tries to patch the incoming and outgoing velocity using an unpowered gravity assist. If however
 * the required bending angle cannot be obtained, the deltaV required to patch this is calculated.
 * Likewise the deltaV required to patch the hyperbolic excess velocities is calculated. The sum
 * of the two is the total deltaV required to perform the maneuver and is returned.
 * \param centralBodyGravitationalParameter Gravitational parameter of the swing-by body.[m^3 s^-2]
 * \param centralBodyVelocity Heliocentric velocity of the swing-by body.                  [m s^-1]
 * \param incomingVelocity Heliocentric velocity of the spacecraft before the swing-by.    [m s^-1]
 * \param outgoingVelocity Heliocentric velocity of the spacecraft after the swing-by.     [m s^-1]
 * \param smallestPeriapsisDistance Closest allowable distance to the swing-by body.            [m]
 * \param newtonRaphson Pointer to the Newton Raphson that the user wants to use.               [-]
 * \return deltaV The deltaV required for the gravity assist maneuver.                     [m s^-1]
 */
double gravityAssist( const double centralBodyGravitationalParameter,
                      const Eigen::Vector3d& centralBodyVelocity,
                      const Eigen::Vector3d& incomingVelocity,
                      const Eigen::Vector3d& outgoingVelocity,
                      const double smallestPeriapsisDistance,
                      boost::shared_ptr< NewtonRaphson > newtonRaphson =
                                                        boost::make_shared< NewtonRaphson >( ) );

//! Propagate an unpowered gravity assist.
/*!
 * Calculates the outgoing velocity of an unpowered gravity assist. The gravity assist is defined
 * by a 3D rotation angle and the pericenter radius of the swing-by maneuver.
 * \param centralBodyGravitationalParameter Gravitational parameter of the swing-by body.[m^3 s^-2]
 * \param centralBodyVelocity Heliocentric velocity of the swing-by body.                  [m s^-1]
 * \param incomingVelocity Heliocentric velocity of the spacecraft before the swing-by.    [m s^-1]
 * \param rotationAngle Angle defining the rotation due to the swing-by in the 3D plane.      [rad]
 * \param pericenterRadius Pericenter radius of the swing-by maneuver.                          [m]
 * \return outgoingVelocity Heliocentric velocity of the spacecraft after the swing-by.    [m s^-1]
 */
Eigen::Vector3d gravityAssist( const double centralBodyGravitationalParameter,
                               const Eigen::Vector3d& centralBodyVelocity,
                               const Eigen::Vector3d& incomingVelocity,
                               const double rotationAngle,
                               const double pericenterRadius );

//! Propagate a powered gravity assist.
/*!
 * Calculates the outgoing velocity of a powered gravity assist. The gravity assist is defined by
 * a 3D rotation angle, the pericenter radius of the swing-by maneuver and the deltaV magnitude
 * that is applied at the pericenter passage of the gravity assist.
 * \param centralBodyGravitationalParameter Gravitational parameter of the swing-by body.[m^3 s^-2]
 * \param centralBodyVelocity Heliocentric velocity of the swing-by body.                  [m s^-1]
 * \param incomingVelocity Heliocentric velocity of the spacecraft before the swing-by.    [m s^-1]
 * \param rotationAngle Angle defining the rotation due to the swing-by in the 3D plane.      [rad]
 * \param pericenterRadius Pericenter radius of the swing-by maneuver.                          [m]
 * \param deltaV DeltaV magnitude of the gravity assist that is applied at pericenter      [m s^-1]
 * \return outgoingVelocity Heliocentric velocity of the spacecraft after the swing-by.    [m s^-1]
 */
Eigen::Vector3d gravityAssist( const double centralBodyGravitationalParameter,
                               const Eigen::Vector3d& centralBodyVelocity,
                               const Eigen::Vector3d& incomingVelocity,
                               const double rotationAngle,
                               const double pericenterRadius,
                               const double deltaV );

//! Gravity assist functions class.
/*!
 * This class contains the functions required by the root-finders in the gravity assist Delta-V
 * computations to find the eccentricity.
 */
class GravityAssistFunctions
{
public:

    //! Constructor with immediate definition of parameters.
    /*!
     * Constructor that sets all the parameters in the velocity-effect functions for use in the
     * Newton Raphson rootfinder.
     */
    GravityAssistFunctions ( const double incomingSemiMajorAxis,
                             const double outgoingSemiMajorAxis,
                             const double bendingAngle )
        : incomingSemiMajorAxis_( incomingSemiMajorAxis),
          outgoingSemiMajorAxis_( outgoingSemiMajorAxis ),
          bendingAngle_ ( bendingAngle )
    { }

    //! Compute velocity-effect.
    /*!
     * Computes velocity-effect delta-V. This function is used by the Newton-Raphson root-finder.
     * \param incomingEccentricity Incoming eccentricity.
     * \return Velocity-effect at defined eccentricity.
     * \sa NewtonRaphson().
     */
    double computeVelocityEffectFunction( double& incomingEccentricity );

    //! Compute first-derivative of velocity-effect.
    /*!
     * Computes first-derivative of velocity-effect delta-V. This function is used by the
     * Newton-Raphson root-finder
     * \param incomingEccentricity Incoming eccentricity.
     * \return Value of first derivative of root-finder function at defined eccentricity.
     * \sa NewtonRapshon().
     */
    double computeFirstDerivativeVelocityEffect( double& incomingEccentricity );

protected:

private:

    //! Semi-major axis of the incoming hyperbolic leg.
    /*!
     * Semi-major axis of the incoming hyperbolic leg.
     */
    const double incomingSemiMajorAxis_;

    //! Semi-major axis of the outgoing hyperbolic leg.
    /*!
     * Semi-major axis of the outgoing hyperbolic leg.
     */
    const double outgoingSemiMajorAxis_;

    //! Bending angle between the excess velocities.
    /*!
     * Bending angle between the excess velocities.
     */
    const double bendingAngle_;
};

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_GRAVITY_ASSIST_H
