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
 *      110128    E. Iorfida        Added boolean variable that sets the necessity of
 *                                  Newton-Raphson method, particularly in unit test.
 *      110202    J. Melman         Renamed certain parameters and added comments to clarify the
 *                                  code more.
 *      110203    E. Iorfida        Changed some variables names and modified punctuation.
 *      110205    J. Melman         Removed the trailing underscores in some public variables. Some
 *                                  comment rephrasing. Changed and added some notes.
 *      110212    J. Melman         Added a reference to my own thesis.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius, replaced by an element of
 *                                  GeometricShapes.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120508    P. Musegaas       The gravitational parameter is now passed as double.
 *      120530    P. Musegaas       Complete revision. Removed class structure, made it a free
 *                                  function. Added two functions for propagating a gravity assist
 *                                  (powered and unpowered).
 *      120625    P. Musegaas       Minor changes.
 *      120703    T. Secretin       Minor layout changes.
 *      120704    P. Musegaas       Minor change. Reduced negligence of velocity effect from 1 cm/s
 *                                  to 1 micrometer/second.
 *
 *    References
 *      Reference for deltaV computation function:
 *          Melman J. Trajectory optimization for a mission to Neptune and Triton, MSc thesis
 *              report, Delft University of Technology, 2007.
 *      Reference for unpowered gravity assist propagation function:
 *          Conway, B.A., Spacecraft Trajectory Optimization, Chapter 7, Cambridge University
 *              Press, 2010.
 *      Reference for powered gravity assist propagation function:
 *          Musegaas, P., Optimization of Space Trajectories Including Multiple Gravity Assists and
 *              Deep Space Maneuvers, MSc thesis report, Delft University of Technology, 2012.
 *              [unpublished so far].
 *
 *    Notes
 *      Gravity assist and swing-by are two different words for the same thing. The delta-V that is
 *      computed for a powered swing-by has not been proven to be the optimum (lowest) to achieve
 *      the desired geometry of incoming and outgoing hyperbolic legs. Some literature research will
 *      have to be done to look at the alternatives.
 *
 *      Note that the exact implementation of Newton Raphson as root finder should be updated if
 *      someone would want to use a different root finding technique.
 *
 *      Note that a velocity effect deltaV of less than 1 micrometer/second is deemed negligable in
 *      this code.
 *
 */

#include <cmath>

#include <Eigen/Dense>

#include <TudatCore/Mathematics/BasicMathematics/linearAlgebra.h>

#include "Tudat/Astrodynamics/MissionSegments/gravityAssist.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

namespace tudat
{
namespace mission_segments
{

//! Calculate deltaV of a gravity assist.
double gravityAssist( const double centralBodyGravitationalParameter,
                      const Eigen::Vector3d& centralBodyVelocity,
                      const Eigen::Vector3d& incomingVelocity,
                      const Eigen::Vector3d& outgoingVelocity,
                      const double smallestPeriapsisDistance,
                      boost::shared_ptr< NewtonRaphson > newtonRaphson )
{
    // Compute incoming and outgoing hyperbolic excess velocity.
    const Eigen::Vector3d incomingHyperbolicExcessVelocity
            = incomingVelocity - centralBodyVelocity;
    const Eigen::Vector3d outgoingHyperbolicExcessVelocity
            = outgoingVelocity - centralBodyVelocity;

    // Compute absolute values of the hyperbolic excess velocities.
    const double absoluteIncomingExcessVelocity = incomingHyperbolicExcessVelocity.norm( );
    const double absoluteOutgoingExcessVelocity = outgoingHyperbolicExcessVelocity.norm( );

    // Compute bending angle.
    const double bendingAngle = mathematics::linear_algebra::computeAngleBetweenVectors(
                            incomingHyperbolicExcessVelocity, outgoingHyperbolicExcessVelocity );

    // Compute maximum achievable bending angle.
    const double maximumBendingAngle =
            std::asin( 1.0 / ( 1.0 + ( smallestPeriapsisDistance *
                                       absoluteIncomingExcessVelocity *
                                       absoluteIncomingExcessVelocity /
                                       centralBodyGravitationalParameter ) ) ) +
            std::asin( 1.0 / ( 1.0 + ( smallestPeriapsisDistance *
                                       absoluteOutgoingExcessVelocity *
                                       absoluteOutgoingExcessVelocity /
                                       centralBodyGravitationalParameter ) ) );

    // Initialize bending effect deltaV, which is zero, unless extra bending angle is required.
    double bendingEffectDeltaV = 0.0;

    // Check if an additional bending angle is required.
    if ( bendingAngle > maximumBendingAngle )
    {
        // Compute required extra bending angle that cannot be delivered by an unpowered swing-by.
        const double extraBendingAngle = bendingAngle - maximumBendingAngle;

        // Compute necessary delta-V due to bending-effect.
        bendingEffectDeltaV = 2.0 * std::min( absoluteIncomingExcessVelocity,
                                                     absoluteOutgoingExcessVelocity ) *
                                     std::sin( extraBendingAngle / 2.0 );
    }

    // Initialize velocity effect delta V parameter.
    double velocityEffectDeltaV = 0.0;

    // An excess speed difference of less than 1 micrometer/second is deemed to be negligible.
    const double speedTolerance = 1.0e-6;

    // Check if it is necessary to apply a swing-by delta-V due to the effect of a difference in
    // absolute excess velocities.
    if ( std::fabs( absoluteIncomingExcessVelocity - absoluteOutgoingExcessVelocity )
         >= speedTolerance )
    {
        // Compute semi-major axis of hyperbolic legs.
        const double incomingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter /
                                             absoluteIncomingExcessVelocity /
                                             absoluteIncomingExcessVelocity;
        const double outgoingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter /
                                             absoluteOutgoingExcessVelocity /
                                             absoluteOutgoingExcessVelocity;

        // Set the gravity assist function with the variables to perform root finder calculations.
        GravityAssistFunctions gravityAssistFunctions( incomingSemiMajorAxis,
                                                       outgoingSemiMajorAxis, bendingAngle);

        // Set initial guess of the variable computed in Newton-Rapshon method.
        newtonRaphson->setInitialGuessOfRoot( 1.01 );

        // Set the class that contains the functions needed for Newton-Raphson.
        NewtonRaphsonAdaptor< GravityAssistFunctions > newtonRaphsonAdaptorForGravityAssist;
        newtonRaphsonAdaptorForGravityAssist.setClass( &gravityAssistFunctions );

        // Set the functions needed for the Newton Raphson method.
        newtonRaphson->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptorForGravityAssist );
        newtonRaphsonAdaptorForGravityAssist.setPointerToFunction(
                    &GravityAssistFunctions::computeVelocityEffectFunction );
        newtonRaphsonAdaptorForGravityAssist.setPointerToFirstDerivativeFunction(
                    &GravityAssistFunctions::computeFirstDerivativeVelocityEffect );

        // Execute Newton-Raphson method.
        newtonRaphson->execute( );

        // Set incoming hyperbolic leg eccentricity as the output value of Newton-Raphson method.
        const double incomingEccentricity = newtonRaphson->getComputedRootOfFunction( );

        // Compute outgoing hyperbolic leg eccentricity.
        const double outgoingEccentricity = 1.0 - ( incomingSemiMajorAxis /
                                                    outgoingSemiMajorAxis ) *
                                            ( 1.0 - incomingEccentricity );

        // Compute incoming and outgoing velocities at periapsis.
        const double incomingVelocityAtPeriapsis = absoluteIncomingExcessVelocity *
                    std::sqrt( ( incomingEccentricity + 1.0 ) / ( incomingEccentricity - 1.0 ) );
        const double outgoingVelocityAtPeriapsis = absoluteOutgoingExcessVelocity *
                    std::sqrt( ( outgoingEccentricity + 1.0 ) / ( outgoingEccentricity - 1.0 ) );

        // Compute necessary delta-V due to velocity-effect.
        velocityEffectDeltaV = std::fabs( incomingVelocityAtPeriapsis -
                                          outgoingVelocityAtPeriapsis );
    }

    // Compute and return the total delta-V.
    return bendingEffectDeltaV + velocityEffectDeltaV;
}

//! Propagate an unpowered gravity assist.
Eigen::Vector3d gravityAssist( const double centralBodyGravitationalParameter,
                               const Eigen::Vector3d& centralBodyVelocity,
                               const Eigen::Vector3d& incomingVelocity,
                               const double rotationAngle,
                               const double pericenterRadius )
{
    // Calculate the incoming velocity.
    const Eigen::Vector3d relativeIncomingVelocity = incomingVelocity - centralBodyVelocity;
    const double absoluteRelativeIncomingVelocity = relativeIncomingVelocity.norm( );

    // Calculate the eccentricity and bending angle.
    const double eccentricity = 1.0 + pericenterRadius / centralBodyGravitationalParameter *
                            absoluteRelativeIncomingVelocity * absoluteRelativeIncomingVelocity;
    const double bendingAngle = 2.0 * std::asin ( 1.0 / eccentricity );

    // Calculate the unit vectors.
    const Eigen::Vector3d unitVector1 = relativeIncomingVelocity /
                                        absoluteRelativeIncomingVelocity;
    const Eigen::Vector3d unitVector2 = unitVector1.cross( centralBodyVelocity ).normalized( );
    const Eigen::Vector3d unitVector3 = unitVector1.cross( unitVector2 );

    // Calculate the relative outgoing velocity.
    const Eigen::Vector3d relativeOutgoingVelocity = absoluteRelativeIncomingVelocity *
            ( std::cos( bendingAngle ) * unitVector1 + std::sin( bendingAngle ) *
              std::cos( rotationAngle ) * unitVector2 + std::sin( bendingAngle ) *
              std::sin( rotationAngle ) * unitVector3 );

    // Add the relative outgoing velocity to the swing-by body velocity and return it.
    return centralBodyVelocity + relativeOutgoingVelocity;
}

//! Propagate a powered gravity assist.
Eigen::Vector3d gravityAssist( const double centralBodyGravitationalParameter,
                               const Eigen::Vector3d& centralBodyVelocity,
                               const Eigen::Vector3d& incomingVelocity,
                               const double rotationAngle,
                               const double pericenterRadius,
                               const double deltaV )
{
    // Calculate the incoming velocity.
    const Eigen::Vector3d relativeIncomingVelocity = incomingVelocity - centralBodyVelocity;
    const double absoluteRelativeIncomingVelocity = relativeIncomingVelocity.norm( );

    // Calculate the incoming eccentricity and bending angle.
    const double incomingEccentricity = 1.0 + pericenterRadius /
                                        centralBodyGravitationalParameter *
                                        absoluteRelativeIncomingVelocity *
                                        absoluteRelativeIncomingVelocity;
    const double incomingBendingAngle = std::asin ( 1.0 / incomingEccentricity );

    // Calculate the pericenter velocities.
    const double incomingPericenterVelocity = std::sqrt( absoluteRelativeIncomingVelocity *
                                                         absoluteRelativeIncomingVelocity *
                                                         ( incomingEccentricity + 1.0 ) /
                                                         ( incomingEccentricity - 1.0 ) );
    const double outgoingPericenterVelocity = incomingPericenterVelocity + deltaV;

    // Calculate magnitude of the absolute relative outgoing velocity.
    const double absoluteRelativeOutgoingVelocity =
            std::sqrt( outgoingPericenterVelocity * outgoingPericenterVelocity -
                       2.0 * centralBodyGravitationalParameter / pericenterRadius );

    // Calculate the remaining bending angles.
    const double outgoingBendingAngle =
            std::asin ( 1.0 / ( 1.0 + absoluteRelativeOutgoingVelocity *
                                absoluteRelativeOutgoingVelocity * pericenterRadius /
                                centralBodyGravitationalParameter ) );
    const double bendingAngle = incomingBendingAngle + outgoingBendingAngle;

    // Calculate the unit vectors.
    const Eigen::Vector3d unitVector1 = relativeIncomingVelocity /
                                        absoluteRelativeIncomingVelocity;
    const Eigen::Vector3d unitVector2 = unitVector1.cross( centralBodyVelocity ).normalized( );
    const Eigen::Vector3d unitVector3 = unitVector1.cross( unitVector2 );

    // Calculate the relative outgoing velocity.
    const Eigen::Vector3d relativeOutgoingVelocity = absoluteRelativeOutgoingVelocity *
            ( std::cos( bendingAngle ) * unitVector1 + std::sin( bendingAngle ) *
              std::cos( rotationAngle ) * unitVector2 + std::sin( bendingAngle ) *
              std::sin( rotationAngle ) * unitVector3 );

    // Add the relative outgoing velocity to the swing-by body velocity and return it.
    return centralBodyVelocity + relativeOutgoingVelocity;
}

//! Compute velocity-effect.
double GravityAssistFunctions::computeVelocityEffectFunction( double& incomingEccentricity )
{
    return std::asin( 1.0 /incomingEccentricity )
            + std::asin( 1.0 / ( 1.0 - incomingSemiMajorAxis_ / outgoingSemiMajorAxis_ *
                                 ( 1.0 -incomingEccentricity ) ) ) - bendingAngle_;
}

//! Compute first-derivative of velocity-effect.
double GravityAssistFunctions::computeFirstDerivativeVelocityEffect(
        double& incomingEccentricity )
{
    const double eccentricitySquareMinusOne_ = incomingEccentricity * incomingEccentricity - 1.0;
    const double semiMajorAxisRatio_ = incomingSemiMajorAxis_ / outgoingSemiMajorAxis_ ;
    const double bParameter_ = 1.0 - semiMajorAxisRatio_ * ( 1.0 - incomingEccentricity );

    return -1.0 / ( incomingEccentricity * std::sqrt( eccentricitySquareMinusOne_ ) ) -
           semiMajorAxisRatio_ / ( bParameter_ * std::sqrt( bParameter_ * bParameter_ - 1.0 ) );
}

} // namespace mission_segments
} // namespace tudat
