/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      References for deltaV computation function:
 *          Melman J. Trajectory optimization for a mission to Neptune and Triton, MSc thesis
 *              report, Delft University of Technology, 2007.
 *          Musegaas, P., Optimization of Space Trajectories Including Multiple Gravity Assists and
 *              Deep Space Maneuvers, MSc thesis report, Delft University of Technology, 2012.
 *              [unpublished so far].
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
 *      Note that by default a velocity effect deltaV of less than 1 micrometer/second is deemed
 *      negligable in this code. This value can be set though.
 *
 */

#include <cmath>

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;


#include <Eigen/Dense>

#include "tudat/math/basic/linearAlgebra.h"

#include "tudat/astro/mission_segments/gravityAssist.h"
#include "tudat/math/basic/functionProxy.h"
#include "tudat/math/root_finders/createRootFinder.h"


namespace tudat
{
namespace mission_segments
{

using namespace root_finders;

double calculateUnpoweredGravityAssistMaximumBendingAngle(
        const double smallestPeriapsisDistance,
        const double absoluteIncomingExcessVelocity,
        const double absoluteOutgoingExcessVelocity,
        const double centralBodyGravitationalParameter )
{
    // Compute maximum achievable bending angle.
    return std::asin( 1.0 / ( 1.0 + ( smallestPeriapsisDistance *
                                      absoluteIncomingExcessVelocity *
                                      absoluteIncomingExcessVelocity /
                                      centralBodyGravitationalParameter ) ) ) +
            std::asin( 1.0 / ( 1.0 + ( smallestPeriapsisDistance *
                                       absoluteOutgoingExcessVelocity *
                                       absoluteOutgoingExcessVelocity /
                                       centralBodyGravitationalParameter ) ) );
}

double calculateUnpoweredGravityAssistPericenter(
        const double absoluteIncomingSemiMajorAxis,
        const double absoluteOutgoingSemiMajorAxis,
        const double bendingAngle,
        const double initialGuess,
        const RootFinderPointer rootFinder )
{

    // Set the gravity assist function with the variables to perform root finder calculations.
    PericenterFindingFunctions pericenterFindingFunctions( absoluteIncomingSemiMajorAxis,
                                                           absoluteOutgoingSemiMajorAxis,
                                                           bendingAngle);

    // Create an object containing the function of which we whish to obtain the root from.
    basic_mathematics::UnivariateProxyPointer rootFunction = std::make_shared<
            basic_mathematics::UnivariateProxy >(
                std::bind( &PericenterFindingFunctions::computePericenterRadiusFunction,
                           pericenterFindingFunctions, std::placeholders::_1 ) );

    // Add the first derivative of the root function.
    rootFunction->addBinding( -1, std::bind( &PericenterFindingFunctions::
                                             computeFirstDerivativePericenterRadiusFunction,
                                             pericenterFindingFunctions, std::placeholders::_1 ) );

    // Set pericenter radius based on result of Newton-Raphson root-finding algorithm.
    return rootFinder->execute( rootFunction, initialGuess );
}

double calculateGravityAssistDeltaVThroughPericenter(
        const double centralBodyGravitationalParameter,
        const double absoluteIncomingExcessVelocity,
        const double absoluteOutgoingExcessVelocity,
        const double bendingAngle,
        const double initialGuess,
        const RootFinderPointer rootFinder )
{
    // Compute semi-major axis of hyperbolic legs. This is the absolute semi-major axis, because
    // it will otherwisely result in the root of a negative function for various cases during
    // the rootfinding process.
    const double absoluteIncomingSemiMajorAxis = 1.0 * centralBodyGravitationalParameter /
            absoluteIncomingExcessVelocity /
            absoluteIncomingExcessVelocity;
    const double absoluteOutgoingSemiMajorAxis = 1.0 * centralBodyGravitationalParameter /
            absoluteOutgoingExcessVelocity /
            absoluteOutgoingExcessVelocity;

    // Calculate pericenter radius
    const double pericenterRadius = calculateUnpoweredGravityAssistPericenter(
                absoluteIncomingSemiMajorAxis,
                absoluteOutgoingSemiMajorAxis,
                bendingAngle,
                initialGuess,
                rootFinder );

    // Compute incoming hyperbolic leg eccentricity.
    const double incomingEccentricity = 1.0 + pericenterRadius / absoluteIncomingSemiMajorAxis;

    // Compute outgoing hyperbolic leg eccentricity.
    const double outgoingEccentricity = 1.0 + pericenterRadius / absoluteOutgoingSemiMajorAxis;

    // Compute incoming and outgoing velocities at periapsis.
    const double incomingVelocityAtPeriapsis = absoluteIncomingExcessVelocity *
            std::sqrt( ( incomingEccentricity + 1.0 ) / ( incomingEccentricity - 1.0 ) );
    const double outgoingVelocityAtPeriapsis = absoluteOutgoingExcessVelocity *
            std::sqrt( ( outgoingEccentricity + 1.0 ) / ( outgoingEccentricity - 1.0 ) );

    // Compute necessary delta-V due to velocity-effect.
    return std::fabs( incomingVelocityAtPeriapsis -
                      outgoingVelocityAtPeriapsis );
}

double calculateGravityAssistDeltaVThroughEccentricity(
        const double centralBodyGravitationalParameter,
        const double absoluteIncomingExcessVelocity,
        const double absoluteOutgoingExcessVelocity,
        const double bendingAngle,
        const RootFinderPointer rootFinder )
{
    // Compute semi-major axis of hyperbolic legs.
    const double incomingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter /
            absoluteIncomingExcessVelocity /
            absoluteIncomingExcessVelocity;
    const double outgoingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter /
            absoluteOutgoingExcessVelocity /
            absoluteOutgoingExcessVelocity;

    // Set the gravity assist function with the variables to perform root finder calculations.
    EccentricityFindingFunctions eccentricityFindingFunctions( incomingSemiMajorAxis,
                                                               outgoingSemiMajorAxis,
                                                               bendingAngle );

    // Create an object containing the function of which we whish to obtain the root from.
    basic_mathematics::UnivariateProxyPointer rootFunction = std::make_shared<
            basic_mathematics::UnivariateProxy >(
                std::bind( &EccentricityFindingFunctions::
                           computeIncomingEccentricityFunction,
                           eccentricityFindingFunctions, std::placeholders::_1 ) );

    // Add the first derivative of the root function.
    rootFunction->addBinding( -1, std::bind(
                                  &EccentricityFindingFunctions::
                                  computeFirstDerivativeIncomingEccentricityFunction,
                                  eccentricityFindingFunctions, std::placeholders::_1 ) );

    // Initialize incoming eccentricity.
    double incomingEccentricity = TUDAT_NAN;

    // Set initial guess of the variable computed in Newton-Rapshon method.
    if ( ( absoluteOutgoingExcessVelocity / absoluteIncomingExcessVelocity ) < 100.0 )
    {
        // In these cases the very low estimate (which is given under else) may in some cases
        // result in no convergence. Hence a higher value of 1.01 is necessary. This will not
        // result in 'going through' 1.0 as mentioned below, because the eccentricity in these
        // cases is always high!
        try
        {
            incomingEccentricity = rootFinder->execute( rootFunction, 1.0 + 1.0e-2 );
        }
        catch(std::runtime_error&)
        {
            root_finders::RootFinderPointer rootFinder_temp
                    = std::make_shared< root_finders::Bisection< > >( 1.0e-12, 1000 ) ;
            incomingEccentricity = rootFinder_temp->execute( rootFunction, 1.0 + 1.0e-2 );

        }
    }

    else
    {
        // This is set to a value that is close to 1.0. This is more robust than higher values,
        // because for those higher values Newton Raphson sometimes 'goes through' 1.0. This
        // results in NaN values for the derivative of the eccentricity finding function.
        try
        {
            incomingEccentricity = rootFinder->execute( rootFunction, 1.0 + 1.0e-10 );
        }
        catch(std::runtime_error&)
        {
            root_finders::RootFinderPointer rootFinder_temp
                    = std::make_shared< root_finders::Bisection< > >( 1.0e-12, 1000 ) ;
            incomingEccentricity = rootFinder_temp->execute( rootFunction, 1.0 + 1.0e-10 );

        }
    }

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
    return std::fabs( incomingVelocityAtPeriapsis -
                      outgoingVelocityAtPeriapsis );
}
//! Calculate deltaV of a gravity assist.
double calculateGravityAssistDeltaV( const double centralBodyGravitationalParameter,
                                     const Eigen::Vector3d& centralBodyVelocity,
                                     const Eigen::Vector3d& incomingVelocity,
                                     const Eigen::Vector3d& outgoingVelocity,
                                     const double smallestPeriapsisDistance,
                                     const bool useEccentricityInsteadOfPericenter,
                                     const double speedTolerance,
                                     RootFinderPointer rootFinder )
{
    using basic_mathematics::UnivariateProxyPointer;
    using basic_mathematics::UnivariateProxy;
    
    // Compute incoming and outgoing hyperbolic excess velocity.
    const Eigen::Vector3d incomingHyperbolicExcessVelocity
            = incomingVelocity - centralBodyVelocity;
    const Eigen::Vector3d outgoingHyperbolicExcessVelocity
            = outgoingVelocity - centralBodyVelocity;
    
    // Compute absolute values of the hyperbolic excess velocities.
    const double absoluteIncomingExcessVelocity = incomingHyperbolicExcessVelocity.norm( );
    const double absoluteOutgoingExcessVelocity = outgoingHyperbolicExcessVelocity.norm( );
    
    // Compute bending angle.
    double bendingAngle = linear_algebra::computeAngleBetweenVectors(
                incomingHyperbolicExcessVelocity, outgoingHyperbolicExcessVelocity );
    
    // Compute maximum achievable bending angle.
    const double maximumBendingAngle =
            calculateUnpoweredGravityAssistMaximumBendingAngle(
                smallestPeriapsisDistance, absoluteIncomingExcessVelocity,
                absoluteOutgoingExcessVelocity, centralBodyGravitationalParameter );
    
    // Initialize bending effect deltaV, which is zero, unless extra bending angle is required.
    double bendingEffectDeltaV = 0.0;
    
    // Initialize velocity effect delta V parameter.
    double velocityEffectDeltaV = 0.0;
    
    // Check if an additional bending angle is required. If so, the additional bending angle
    // maneuver has to be performed. Also the pericenter radius will be the minimum pericenter
    // radius to obtain the largest possible bending angle 'for free'. Hence no root finding is
    // required for this case. As noted above, this may not be ideal for all cases. (for cases
    // in which the excess velocities are relatively small)
    if ( bendingAngle > maximumBendingAngle )
    {
        // Compute required extra bending angle that cannot be delivered by an unpowered swing-by.
        const double extraBendingAngle = bendingAngle - maximumBendingAngle;
        
        // Compute necessary delta-V due to bending-effect.
        bendingEffectDeltaV = 2.0 * std::min( absoluteIncomingExcessVelocity,
                                              absoluteOutgoingExcessVelocity ) *
                std::sin( extraBendingAngle / 2.0 );
        
        // This means the pericenter radius is now equal to the smallest pericenter radius, to
        // ensure the largest possible bending angle.
        const double pericenterRadius = smallestPeriapsisDistance;
        
        // Compute semi-major axis of hyperbolic legs.
        const double incomingSemiMajorAxis = -centralBodyGravitationalParameter /
                ( absoluteIncomingExcessVelocity * absoluteIncomingExcessVelocity );
        const double outgoingSemiMajorAxis = -centralBodyGravitationalParameter /
                ( absoluteOutgoingExcessVelocity * absoluteOutgoingExcessVelocity );
        
        // Compute incoming hyperbolic leg eccentricity.
        const double incomingEccentricity = 1.0 - pericenterRadius / incomingSemiMajorAxis;
        
        // Compute outgoing hyperbolic leg eccentricity.
        const double outgoingEccentricity = 1.0 - pericenterRadius / outgoingSemiMajorAxis;
        
        // Compute incoming and outgoing velocities at periapsis.
        const double incomingVelocityAtPeriapsis = absoluteIncomingExcessVelocity *
                std::sqrt( ( incomingEccentricity + 1.0 ) / ( incomingEccentricity - 1.0 ) );
        const double outgoingVelocityAtPeriapsis = absoluteOutgoingExcessVelocity *
                std::sqrt( ( outgoingEccentricity + 1.0 ) / ( outgoingEccentricity - 1.0 ) );
        
        // Compute necessary delta-V due to velocity-effect.
        velocityEffectDeltaV = std::fabs( incomingVelocityAtPeriapsis -
                                          outgoingVelocityAtPeriapsis );
    }
    
    else if ( ( std::fabs( absoluteIncomingExcessVelocity - absoluteOutgoingExcessVelocity )
                <= speedTolerance ) )
    {
        // In this case no maneuver has to be performed. Hence no iteration is performed, and the
        // delta V is simply kept at 0.0.
    }
    
    // Here the required maneuver to patch the incoming and outgoing excess velocities is
    // calculated. In this implementation, the eccentricity will be used as iteration parameter.
    else if ( useEccentricityInsteadOfPericenter )
    {
//<<<<<<< HEAD
        // Compute semi-major axis of hyperbolic legs.
        const double incomingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter /
                                             absoluteIncomingExcessVelocity /
                                             absoluteIncomingExcessVelocity;
        const double outgoingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter /
                                             absoluteOutgoingExcessVelocity /
                                             absoluteOutgoingExcessVelocity;

        // Set the gravity assist function with the variables to perform root finder calculations.
        EccentricityFindingFunctions eccentricityFindingFunctions( incomingSemiMajorAxis,
                                                                   outgoingSemiMajorAxis,
                                                                   bendingAngle );

        // Create an object containing the function of which we whish to obtain the root from.
        UnivariateProxyPointer rootFunction = std::make_shared< UnivariateProxy >(
                    std::bind( &EccentricityFindingFunctions::
                                 computeIncomingEccentricityFunction,
                                 eccentricityFindingFunctions, std::placeholders::_1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding( -1, std::bind(
                                      &EccentricityFindingFunctions::
                                      computeFirstDerivativeIncomingEccentricityFunction,
                                      eccentricityFindingFunctions, std::placeholders::_1 ) );

        // Initialize incoming eccentricity.
        double incomingEccentricity = TUDAT_NAN;

        // Set initial guess of the variable computed in Newton-Rapshon method.
        if ( ( absoluteOutgoingExcessVelocity / absoluteIncomingExcessVelocity ) < 100.0 )
        {
            // In these cases the very low estimate (which is given under else) may in some cases
            // result in no convergence. Hence a higher value of 1.01 is necessary. This will not
            // result in 'going through' 1.0 as mentioned below, because the eccentricity in these
            // cases is always high!
            try
            {
                incomingEccentricity = rootFinder->execute( rootFunction, 1.0 + 1.0e-2 );
            }
            catch(std::runtime_error&)
            {
                root_finders::RootFinderPointer rootFinder_temp = root_finders::createRootFinder(
                                        root_finders::bisectionRootFinderSettings(
                                            TUDAT_NAN, 10E-12,
                                            TUDAT_NAN, 1000 ) );

                incomingEccentricity = rootFinder_temp->execute( rootFunction, 1.0 + 1.0e-2 );

            }
        }

        else
        {
            // This is set to a value that is close to 1.0. This is more robust than higher values,
            // because for those higher values Newton Raphson sometimes 'goes through' 1.0. This
            // results in NaN values for the derivative of the eccentricity finding function.
            try
            {
                incomingEccentricity = rootFinder->execute( rootFunction, 1.0 + 1.0e-10 );
            }
            catch(std::runtime_error&)
            {
                root_finders::RootFinderPointer rootFinder_temp
                  = std::make_shared< root_finders::Bisection< > >( 1.0e-12, 1000 ) ;
                incomingEccentricity = rootFinder_temp->execute( rootFunction, 1.0 + 1.0e-10 );

            }
        }

        // Compute outgoing hyperbolic leg eccentricity.
        const double outgoingEccentricity = 1.0 - ( incomingSemiMajorAxis /
                                                    outgoingSemiMajorAxis ) *
                                            ( 1.0 - incomingEccentricity );

        // Compute incoming and outgoing velocities at periapsis.
        const double incomingVelocityAtPeriapsis = absoluteIncomingExcessVelocity *
                std::sqrt( ( incomingEccentricity + 1.0 ) / ( incomingEccentricity - 1.0 ) );
        const double outgoingVelocityAtPeriapsis = absoluteOutgoingExcessVelocity *
                std::sqrt( ( outgoingEccentricity + 1.0 ) / ( outgoingEccentricity - 1.0 ) );

//=======
//>>>>>>> feature/mga_dsm_refactor
        // Compute necessary delta-V due to velocity-effect.
        velocityEffectDeltaV = calculateGravityAssistDeltaVThroughEccentricity(
                    centralBodyGravitationalParameter, absoluteIncomingExcessVelocity, absoluteOutgoingExcessVelocity,
                    bendingAngle, rootFinder );
    }
    
    // Here the required maneuver to patch the incoming and outgoing excess velocities is
    // calculated. In this implementation, the pericenter radius will be used as iteration
    // parameter.
    else
    {
        // Compute necessary delta-V due to velocity-effect.
        velocityEffectDeltaV = calculateGravityAssistDeltaVThroughPericenter(
                    centralBodyGravitationalParameter, absoluteIncomingExcessVelocity, absoluteOutgoingExcessVelocity,
                    bendingAngle, smallestPeriapsisDistance, rootFinder );
    }
    
    // Compute and return the total delta-V.
    return bendingEffectDeltaV + velocityEffectDeltaV;
}

//! Propagate an unpowered gravity assist.
Eigen::Vector3d calculateUnpoweredGravityAssistOutgoingVelocity(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& centralBodyVelocity,
        const Eigen::Vector3d& incomingVelocity,
        const double outgoingVelocityRotationAngle,
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
                                                       std::cos(outgoingVelocityRotationAngle ) * unitVector2 + std::sin(bendingAngle ) *
                                                                                                                std::sin(outgoingVelocityRotationAngle ) * unitVector3 );
    
    // Add the relative outgoing velocity to the swing-by body velocity and return it.
    return centralBodyVelocity + relativeOutgoingVelocity;
}

//! Propagate a powered gravity assist.
Eigen::Vector3d calculatePoweredGravityAssistOutgoingVelocity(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& centralBodyVelocity,
        const Eigen::Vector3d& incomingVelocity,
        const double outgoingVelocityRotationAngle,
        const double pericenterRadius,
        const double deltaV )
{
    // Calculate the incoming velocity.
    const Eigen::Vector3d relativeIncomingVelocity = incomingVelocity - centralBodyVelocity;
    const double absoluteRelativeIncomingVelocity = relativeIncomingVelocity.norm( );

    if ( absoluteRelativeIncomingVelocity == 0 )
    {
        throw std::runtime_error( "Incoming excess velocity at swingby must be different from 0. " );
    }
    if ( pericenterRadius <= 0 )
    {
        throw std::runtime_error( "Error when computing powered swingby: pericenter radius (" + std::to_string(pericenterRadius) +
            ") must be larger than zero.");
    }
    
    // Calculate the incoming eccentricity and bending angle.
    const double incomingEccentricity = 1.0 + pericenterRadius /
            centralBodyGravitationalParameter *
            absoluteRelativeIncomingVelocity *
            absoluteRelativeIncomingVelocity;
    const double incomingBendingAngle = std::asin ( 1.0 / incomingEccentricity );
    
    // Calculate the pericenter velocities.
    // Deal with limit case where eccentricity is infinity by setting the velocity to absoluteRelativeIncomingVelocity
    double incomingPericenterVelocity;
    if ( incomingEccentricity < std::numeric_limits< double >::infinity() )
    {
        incomingPericenterVelocity = std::sqrt( absoluteRelativeIncomingVelocity *
                                                 absoluteRelativeIncomingVelocity *
                                                 ( incomingEccentricity + 1.0 ) /
                                                 ( incomingEccentricity - 1.0 ) );
    }
    else
    {
        incomingPericenterVelocity = absoluteRelativeIncomingVelocity;
    }

    const double outgoingPericenterVelocity = incomingPericenterVelocity + deltaV;
    
    // Calculate magnitude of the absolute relative outgoing velocity.
    const double absoluteRelativeOutgoingVelocity =
            std::sqrt( outgoingPericenterVelocity * outgoingPericenterVelocity -
                       2.0 * centralBodyGravitationalParameter / pericenterRadius );

    if ( !(absoluteRelativeOutgoingVelocity == absoluteRelativeOutgoingVelocity) )
    {
        throw std::runtime_error( "Invalid gravity assist: there is no feasible relative incoming/outgoing velocity." );
    }
    
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
            ( std::cos( bendingAngle ) * unitVector1 +
              std::sin( bendingAngle ) *
              std::cos(outgoingVelocityRotationAngle ) * unitVector2 +
              std::sin( bendingAngle ) *
              std::sin(outgoingVelocityRotationAngle ) * unitVector3 );
    
    // Add the relative outgoing velocity to the swing-by body velocity and return it.
    return centralBodyVelocity + relativeOutgoingVelocity;
}

//! Backward propagate a powered gravity assist.
Eigen::Vector3d calculatePoweredGravityAssistIncomingVelocity(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& centralBodyVelocity,
        const Eigen::Vector3d& outgoingVelocity,
        const double incomingVelocityRotationAngle,
        const double pericenterRadius,
        const double deltaV )
{
    // To compute the incoming velocity one uses the same equations as for computing the outgoing velocity, the only
    // difference being that the deltaV is negative.
    // Evidently, the rotationAngle has a different physical meaning with respect to the one used when computing the
    // outgoing velocity, since the two are defined with respect to two different frames.
    Eigen::Vector3d incomingVelocity;
    incomingVelocity = calculatePoweredGravityAssistOutgoingVelocity(
            centralBodyGravitationalParameter, centralBodyVelocity, outgoingVelocity, incomingVelocityRotationAngle,
            pericenterRadius, -deltaV);
    return incomingVelocity;
}

//! Compute pericenter radius function.
double PericenterFindingFunctions::computePericenterRadiusFunction( const double pericenterRadius )
{
    return std::asin( absoluteIncomingSemiMajorAxis_ / ( absoluteIncomingSemiMajorAxis_ +
                                                         pericenterRadius ) ) +
            std::asin( absoluteOutgoingSemiMajorAxis_ / ( absoluteOutgoingSemiMajorAxis_ +
                                                          pericenterRadius ) ) - bendingAngle_;
}

//! Compute first-derivative of the pericenter radius function.
double PericenterFindingFunctions::computeFirstDerivativePericenterRadiusFunction(
        const double pericenterRadius )
{
    return -absoluteIncomingSemiMajorAxis_ / ( absoluteIncomingSemiMajorAxis_ + pericenterRadius )
            / std::sqrt( ( pericenterRadius + 2.0 * absoluteIncomingSemiMajorAxis_ )
                         * pericenterRadius ) -
            absoluteOutgoingSemiMajorAxis_ / ( absoluteOutgoingSemiMajorAxis_ + pericenterRadius )
            / std::sqrt( ( pericenterRadius + 2.0 * absoluteOutgoingSemiMajorAxis_ ) *
                         pericenterRadius );
}

//! Compute incoming eccentricity function.
double EccentricityFindingFunctions::computeIncomingEccentricityFunction(
        const double incomingEccentricity )
{
    return std::asin( 1.0 / incomingEccentricity )
            + std::asin( 1.0 / ( 1.0 - incomingSemiMajorAxis_ / outgoingSemiMajorAxis_ *
                                 ( 1.0 - incomingEccentricity ) ) ) - bendingAngle_;
}

//! Compute first-derivative of the incoming eccentricity function.
double EccentricityFindingFunctions::computeFirstDerivativeIncomingEccentricityFunction(
        const double incomingEccentricity )
{
    const double eccentricitySquareMinusOne_ = incomingEccentricity * incomingEccentricity - 1.0;
    const double semiMajorAxisRatio_ = incomingSemiMajorAxis_ / outgoingSemiMajorAxis_ ;
    const double bParameter_ = 1.0 - semiMajorAxisRatio_ * ( 1.0 - incomingEccentricity );
    
    return -1.0 / ( incomingEccentricity * std::sqrt( eccentricitySquareMinusOne_ ) ) -
            semiMajorAxisRatio_ / ( bParameter_ * std::sqrt( bParameter_ * bParameter_ - 1.0 ) );
}

} // namespace mission_segments
} // namespace tudat
