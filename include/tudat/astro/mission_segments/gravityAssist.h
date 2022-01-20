/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

#include "tudat/math/root_finders/newtonRaphson.h"
#include "tudat/math/root_finders/rootFinder.h"
#include "tudat/math/root_finders/terminationConditions.h"

namespace tudat
{
namespace mission_segments
{

double calculateUnpoweredGravityAssistPericenter(
        const double absoluteIncomingSemiMajorAxis,
        const double absoluteOutgoingSemiMajorAxis,
        const double bendingAngle,
        const double initialGuess,
        const root_finders::RootFinderPointer rootFinder =
        std::make_shared< root_finders::NewtonRaphson< > >( 1.0e-12, 1000 ) );

double calculateGravityAssistDeltaVThroughPericenter(
        const double centralBodyGravitationalParameter,
        const double absoluteIncomingExcessVelocity,
        const double absoluteOutgoingExcessVelocity,
        const double bendingAngle,
        const double initialGuess,
        const root_finders::RootFinderPointer rootFinder =
        std::make_shared< root_finders::NewtonRaphson< > >( 1.0e-12, 1000 ) );

double calculateGravityAssistDeltaVThroughEccentricity(
        const double centralBodyGravitationalParameter,
        const double absoluteIncomingExcessVelocity,
        const double absoluteOutgoingExcessVelocity,
        const double bendingAngle,
        const root_finders::RootFinderPointer rootFinder =
        std::make_shared< root_finders::NewtonRaphson< > >( 1.0e-12, 1000 ) );

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
//<<<<<<< HEAD
//double gravityAssist( const double centralBodyGravitationalParameter,
//                      const Eigen::Vector3d& centralBodyVelocity,
//                      const Eigen::Vector3d& incomingVelocity,
//                      const Eigen::Vector3d& outgoingVelocity,
//                      const double smallestPeriapsisDistance,
//                      const bool useEccentricityInsteadOfPericenter = true,
//                      const double speedTolerance = 1.0e-6,
//                      root_finders::RootFinderPointer rootFinder
//                        = std::make_shared< root_finders::NewtonRaphson< > >( 1.0e-12, 1000 ) );
//=======
double calculateGravityAssistDeltaV(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& centralBodyVelocity,
        const Eigen::Vector3d& incomingVelocity,
        const Eigen::Vector3d& outgoingVelocity,
        const double smallestPeriapsisDistance,
        const bool useEccentricityInsteadOfPericenter = true,
        const double speedTolerance = 1.0e-6,
        root_finders::RootFinderPointer rootFinder
        = std::make_shared< root_finders::NewtonRaphson< > >( 1.0e-12, 1000 ) );
//>>>>>>> feature/mga_dsm_refactor

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
Eigen::Vector3d calculateUnpoweredGravityAssistOutgoingVelocity(
        const double centralBodyGravitationalParameter,
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
Eigen::Vector3d calculatePoweredGravityAssistOutgoingVelocity(
        const double centralBodyGravitationalParameter,
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
typedef std::shared_ptr< PericenterFindingFunctions > PericenterFindingFunctionsPointer;

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
typedef std::shared_ptr< EccentricityFindingFunctions > EccentricityFindingFunctionsPointer;

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_GRAVITY_ASSIST_H
