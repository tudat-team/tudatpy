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
 *      110117    E. Iorfida        File created.
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
 *      120713    P. Musegaas       Added option to iterate on pericenter radius instead of
 *                                  eccentricity. Added separate class for this and a flag in free
 *                                  function.
 *      120813    P. Musegaas       Changed code to new root finding structure. Added option to
 *                                  specify which rootfinder and termination conditions to use.
 *      130121    K. Kumar          Added shared-ptr typedefs.
 *      140117    E. Brandon        Corrected doxygen documentation.
 *
 *    References
 *
 *    Notes
 *      Note that the exact implementation of Newton-Raphson as root finder should be updated if
 *      someone would want to use a different root-finding technique.
 *
 *      By default the eccentricity is used as the iteration procedure. This is because in
 *      optimizing a Cassini-like trajectory, the pericenter radius had about 2-4 NaN values in
 *      100000 times the gravity assist calculation. The eccentricity iteration had no NaN values
 *      for a similar run in which 100000 gravity assist calculations were done. Also the
 *      eccentricity seemed to require slightly less iterations (does not necessarily mean it is
 *      faster or more accurate).
 *
 */

#ifndef TUDAT_GRAVITY_ASSIST_H
#define TUDAT_GRAVITY_ASSIST_H

#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"

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
 * \param useEccentricityInsteadOfPericenter Flag to indicate the iteration procedure for matching
 *                                           the bending angle.                                 [-]
 * \param speedTolerance Tolerance at which the velocity effect deltaV is deemed 0.0.           [-]
 * \param rootFinder Shared-pointer to the rootfinder that is to be used. Default is Newton-Raphson
 *          using 1000 iterations as maximum and 1.0e-12 relative X-tolerance.
 * \return deltaV The deltaV required for the gravity assist maneuver.                     [m s^-1]
 */
double gravityAssist( const double centralBodyGravitationalParameter,
                      const Eigen::Vector3d& centralBodyVelocity,
                      const Eigen::Vector3d& incomingVelocity,
                      const Eigen::Vector3d& outgoingVelocity,
                      const double smallestPeriapsisDistance,
                      const bool useEccentricityInsteadOfPericenter = true,
                      const double speedTolerance = 1.0e-6,
                      root_finders::RootFinderPointer rootFinder
                        = boost::make_shared< root_finders::NewtonRaphson >( 1.0e-12, 1000 ) );

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

//! Pericenter finding functions class.
/*!
 * This class contains the functions required by the root-finders to find the pericenter radius in
 * the gravity assist function to find the deltaV.
 */
class PericenterFindingFunctions
{
public:

    //! Constructor with immediate definition of parameters.
    /*!
     * Constructor that sets all the parameters in the pericenter finding functions for use in the
     * Newton Raphson rootfinder. Note that absolute values of the semi-major axis are required due
     * to the rootfinding process.
     *
     * \param absoluteIncomingSemiMajorAxis The absolute semi-major axis of the incoming hyperbolic
     *          leg.                                                                            [m]
     * \param absoluteOutgoingSemiMajorAxis The absolute semi-major axis of the outgoing hyperbolic
     *          leg.                                                                            [m]
     * \param bendingAngle The bending angle between the excess velocities.                   [rad]
     */
    PericenterFindingFunctions ( const double absoluteIncomingSemiMajorAxis,
                                 const double absoluteOutgoingSemiMajorAxis,
                                 const double bendingAngle )
        : absoluteIncomingSemiMajorAxis_( absoluteIncomingSemiMajorAxis ),
          absoluteOutgoingSemiMajorAxis_( absoluteOutgoingSemiMajorAxis ),
          bendingAngle_( bendingAngle )
    { }

    //! Compute pericenter radius function.
    /*!
     * Computes pericenter radius function. This function is used by the Newton-Raphson root-finder
     * to find the pericenter radius that matches the bending angle required in the gravity assist.
     * \param pericenterRadius Pericenter radius.
     * \return Pericenter radius root finding function value.
     * \sa NewtonRaphson().
     */
    double computePericenterRadiusFunction( const double pericenterRadius );

    //! Compute first-derivative of the pericenter radius function.
    /*!
     * Computes the first-derivative of the pericenter radius function. This function is used by
     * the Newton-Raphson root-finder to find the pericenter radius that matches the bending angle
     * required in the gravity assist.
     * \param pericenterRadius Pericenter radius.
     * \return Pericenter radius root finding function first-derivative value.
     * \sa NewtonRapshon().
     */
    double computeFirstDerivativePericenterRadiusFunction( const double pericenterRadius );

protected:

private:

    //! The absolute semi-major axis of the incoming hyperbolic leg.
    /*!
     * The absolute semi-major axis of the incoming hyperbolic leg. The absolute value is required
     * because otherwisely the first derivative of the pericenter radius finding function will
     * compute the root of a negative value.
     */
    const double absoluteIncomingSemiMajorAxis_;

    //! The absolute semi-major axis of the outgoing hyperbolic leg.
    /*!
     * The absolute semi-major axis of the outgoing hyperbolic leg. The absolute value is required
     * because otherwisely the first derivative of the pericenter radius finding function will
     * compute the root of a negative value.
     */
    const double absoluteOutgoingSemiMajorAxis_;

    //! Bending angle between the excess velocities.
    /*!
     * Bending angle between the excess velocities.
     */
    const double bendingAngle_;
};

//! Typedef for shared-pointer to PericenterFindingFunctions object.
typedef boost::shared_ptr< PericenterFindingFunctions > PericenterFindingFunctionsPointer;

//! Eccentricity finding functions class.
/*!
 * This class contains the functions required by the root-finders to find the incoming eccentricity
 * in the gravity assist function to find the deltaV.
 */
class EccentricityFindingFunctions
{
public:

    //! Constructor with immediate definition of parameters.
    /*!
     * Constructor that sets all the parameters in the eccentricity finding functions for use in
     * the Newton Raphson rootfinder.
     * \param incomingSemiMajorAxis The semi-major axis of the incomming hyperbolic leg.        [m]
     * \param outgoingSemiMajorAxis The semi-major axis of the outgoing hyperbolic leg.         [m]
     * \param bendingAngle The bending angle between the excess velocities.                   [rad]
     */
    EccentricityFindingFunctions( const double incomingSemiMajorAxis,
                                  const double outgoingSemiMajorAxis,
                                  const double bendingAngle )
        : incomingSemiMajorAxis_( incomingSemiMajorAxis),
          outgoingSemiMajorAxis_( outgoingSemiMajorAxis ),
          bendingAngle_( bendingAngle )
    { }

    //! Compute incoming eccentricity function.
    /*!
     * Computes incoming eccentricity function. This function is used by the Newton-Raphson root-
     * finder to find the incoming eccentricity that matches the bending angle required in the
     * gravity assist.
     * \param incomingEccentricity Incoming eccentricity.
     * \return Incoming eccentricity root finding function value.
     * \sa NewtonRaphson().
     */
    double computeIncomingEccentricityFunction( const double incomingEccentricity );

    //! Compute first-derivative of the incoming eccentricity function.
    /*!
     * Computes the first-derivative of the incoming eccentricity function. This function is used
     * by the Newton-Raphson root-finder to find the incoming eccentricity that matches the bending
     * angle required in the gravity assist.
     * \param incomingEccentricity Incoming eccentricity.
     * \return Incoming eccentricity root finding function first-derivative value.
     * \sa NewtonRapshon().
     */
    double computeFirstDerivativeIncomingEccentricityFunction( const double incomingEccentricity );

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

//! Typedef for shared-pointer to EccentricityFindingFunctions object.
typedef boost::shared_ptr< EccentricityFindingFunctions > EccentricityFindingFunctionsPointer;

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_GRAVITY_ASSIST_H
