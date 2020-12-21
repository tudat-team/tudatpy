#include "gravityAssist.h"
#include <tudat/basics/basicTypedefs.h>
#include <tudat/math/root_finders.h>
#include <tudat/math/basic.h>

using namespace tudat::root_finders;
namespace trf = tudat::root_finders;

namespace tudat {

namespace prototype {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////  MODIFIED POWERED GRAVITY ASSIST                                            ///////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Calculate deltaV of a gravity assist.
GravityAssistParameters calculateGravityAssistParameters(
    const double centralBodyGravitationalParameter,
    const Eigen::Vector3d &centralBodyVelocityOnEntrySOI,
    const Eigen::Vector3d &centralBodyVelocityOnExitSOI,
    const Eigen::Vector3d &incomingVelocity,
    const Eigen::Vector3d &outgoingVelocity,
    const double smallestPeriapsisDistance,
    const bool useEccentricityInsteadOfPericenter,
    const double speedTolerance,
    RootFinderPointer rootFinder)

{
  GravityAssistParameters saveParameters;
  using tudat::basic_mathematics::UnivariateProxy;
  using tudat::basic_mathematics::UnivariateProxyPointer;

  // Compute incoming and outgoing hyperbolic excess velocity.
  const Eigen::Vector3d incomingHyperbolicExcessVelocity = incomingVelocity - centralBodyVelocityOnEntrySOI;
  const Eigen::Vector3d outgoingHyperbolicExcessVelocity = outgoingVelocity - centralBodyVelocityOnExitSOI;

  // Compute absolute values of the hyperbolic excess velocities.
  const double absoluteIncomingExcessVelocity = incomingHyperbolicExcessVelocity.norm();
  const double absoluteOutgoingExcessVelocity = outgoingHyperbolicExcessVelocity.norm();

  // Compute bending angle.
  double bendingAngle = tudat::linear_algebra::computeAngleBetweenVectors(
      incomingHyperbolicExcessVelocity, outgoingHyperbolicExcessVelocity);

  // Compute maximum achievable bending angle.
  const double maximumBendingAngle =
      std::asin(1.0 / (1.0 + (smallestPeriapsisDistance * absoluteIncomingExcessVelocity * absoluteIncomingExcessVelocity / centralBodyGravitationalParameter))) + std::asin(1.0 / (1.0 + (smallestPeriapsisDistance * absoluteOutgoingExcessVelocity * absoluteOutgoingExcessVelocity / centralBodyGravitationalParameter)));

  // Initialize bending effect deltaV, which is zero, unless extra bending angle is required.
  double bendingEffectDeltaV = 0.0;

  // Initialize velocity effect delta V parameter.
  double velocityEffectDeltaV = 0.0;

  // Check if an additional bending angle is required. If so, the additional bending angle
  // maneuver has to be performed. Also the pericenter radius will be the minimum pericenter
  // radius to obtain the largest possible bending angle 'for free'. Hence no root finding is
  // required for this case. As noted above, this may not be ideal for all cases. (for cases
  // in which the excess velocities are relatively small)
  if (bendingAngle > maximumBendingAngle) {
    // Compute required extra bending angle that cannot be delivered by an unpowered swing-by.
    const double extraBendingAngle = bendingAngle - maximumBendingAngle;

    // Compute necessary delta-V due to bending-effect.
    bendingEffectDeltaV = 2.0 * std::min(absoluteIncomingExcessVelocity, absoluteOutgoingExcessVelocity) * std::sin(extraBendingAngle / 2.0);

    // This means the pericenter radius is now equal to the smallest pericenter radius, to
    // ensure the largest possible bending angle.
    const double pericenterRadius = smallestPeriapsisDistance;

    // Compute semi-major axis of hyperbolic legs.
    const double incomingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter / absoluteIncomingExcessVelocity / absoluteIncomingExcessVelocity;
    const double outgoingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter / absoluteOutgoingExcessVelocity / absoluteOutgoingExcessVelocity;

    // Compute incoming hyperbolic leg eccentricity.
    const double incomingEccentricity = 1 - pericenterRadius / incomingSemiMajorAxis;

    // Compute outgoing hyperbolic leg eccentricity.
    const double outgoingEccentricity = 1 - pericenterRadius / outgoingSemiMajorAxis;

    // Compute incoming and outgoing velocities at periapsis.
    const double incomingVelocityAtPeriapsis = absoluteIncomingExcessVelocity * std::sqrt((incomingEccentricity + 1.0) / (incomingEccentricity - 1.0));
    const double outgoingVelocityAtPeriapsis = absoluteOutgoingExcessVelocity * std::sqrt((outgoingEccentricity + 1.0) / (outgoingEccentricity - 1.0));

    // Compute necessary delta-V due to velocity-effect.
    velocityEffectDeltaV = std::fabs(incomingVelocityAtPeriapsis - outgoingVelocityAtPeriapsis);

    // SAVE TO POWERED GRAVITY ASSIST PARAMETERS
    saveParameters.incomingSemiMajorAxis = incomingSemiMajorAxis;
    saveParameters.outgoingSemiMajorAxis = outgoingSemiMajorAxis;

    saveParameters.incomingVelocityAtPeriapsis = incomingVelocityAtPeriapsis;
    saveParameters.outgoingVelocityAtPeriapsis = outgoingVelocityAtPeriapsis;

    saveParameters.incomingEccentricity = incomingEccentricity;
    saveParameters.outgoingEccentricity = outgoingEccentricity;

  } else if ((std::fabs(absoluteIncomingExcessVelocity - absoluteOutgoingExcessVelocity)
              <= speedTolerance)) {
    // In this case no maneuver has to be performed. Hence no iteration is performed, and the
    // delta V is simply kept at 0.0.

    // SAVE TO POWERED GRAVITY ASSIST PARAMETERS
    saveParameters.incomingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter / absoluteIncomingExcessVelocity / absoluteIncomingExcessVelocity;

    // SAVE TO POWERED GRAVITY ASSIST PARAMETERS
    saveParameters.outgoingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter / absoluteIncomingExcessVelocity / absoluteIncomingExcessVelocity;

    // SAVE TO POWERED GRAVITY ASSIST PARAMETERS
    saveParameters.incomingEccentricity =
        1 - smallestPeriapsisDistance / saveParameters.incomingSemiMajorAxis;

    // SAVE TO POWERED GRAVITY ASSIST PARAMETERS
    saveParameters.outgoingEccentricity =
        1 - smallestPeriapsisDistance / saveParameters.outgoingSemiMajorAxis;

    // Compute incoming and outgoing velocities at periapsis.
    saveParameters.incomingVelocityAtPeriapsis = absoluteIncomingExcessVelocity * std::sqrt((saveParameters.incomingEccentricity + 1.0) / (saveParameters.incomingEccentricity - 1.0));
    saveParameters.outgoingVelocityAtPeriapsis = absoluteOutgoingExcessVelocity * std::sqrt((saveParameters.outgoingEccentricity + 1.0) / (saveParameters.outgoingEccentricity - 1.0));

    //              // TODO: THIS WILL BE AN ISSUE, What is the purpose of having a speed tolerance?
  }
  //
  //                // Here the required maneuver to patch the incoming and outgoing excess velocities is
  //                // calculated. In this implementation, the eccentricity will be used as iteration parameter.
  else if (useEccentricityInsteadOfPericenter) {

    // Compute semi-major axis of hyperbolic legs.
    const double incomingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter / absoluteIncomingExcessVelocity / absoluteIncomingExcessVelocity;
    const double outgoingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter / absoluteOutgoingExcessVelocity / absoluteOutgoingExcessVelocity;

    // Set the gravity assist function with the variables to perform root finder calculations.
    tudat::mission_segments::EccentricityFindingFunctions eccentricityFindingFunctions(incomingSemiMajorAxis,
                                                                                       outgoingSemiMajorAxis,
                                                                                       bendingAngle);

    // Create an object containing the function of which we whish to obtain the root from.
    UnivariateProxyPointer rootFunction = std::make_shared<UnivariateProxy>(
        std::bind(&tudat::mission_segments::EccentricityFindingFunctions::computeIncomingEccentricityFunction,
                  eccentricityFindingFunctions, std::placeholders::_1));

    // Add the first derivative of the root function.
    rootFunction->addBinding(-1, std::bind(&tudat::mission_segments::EccentricityFindingFunctions::computeFirstDerivativeIncomingEccentricityFunction, eccentricityFindingFunctions, std::placeholders::_1));

    // Initialize incoming eccentricity.
    double incomingEccentricity = TUDAT_NAN;

    // Set initial guess of the variable computed in Newton-Rapshon method.
    if ((absoluteOutgoingExcessVelocity / absoluteIncomingExcessVelocity) < 100.0) {
      // In these cases the very low estimate (which is given under else) may in some cases
      // result in no convergence. Hence a higher value of 1.01 is necessary. This will not
      // result in 'going through' 1.0 as mentioned below, because the eccentricity in these
      // cases is always high!
      try {
        incomingEccentricity = rootFinder->execute(rootFunction, 1.0 + 1.0e-2);
      } catch (std::runtime_error) {
        tudat::root_finders::RootFinderPointer rootFinder_temp = std::make_shared<tudat::root_finders::Bisection>(1.0e-12, 1000);
        incomingEccentricity = rootFinder_temp->execute(rootFunction, 1.0 + 1.0e-2);
      }
    } else {
      // This is set to a value that is close to 1.0. This is more robust than higher values,
      // because for those higher values Newton Raphson sometimes 'goes through' 1.0. This
      // results in NaN values for the derivative of the eccentricity finding function.
      try {
        incomingEccentricity = rootFinder->execute(rootFunction, 1.0 + 1.0e-10);
      } catch (std::runtime_error) {
        tudat::root_finders::RootFinderPointer rootFinder_temp = std::make_shared<tudat::root_finders::Bisection>(1.0e-12, 1000);
        incomingEccentricity = rootFinder_temp->execute(rootFunction, 1.0 + 1.0e-10);
      }
    }

    // Compute outgoing hyperbolic leg eccentricity.
    const double outgoingEccentricity = 1.0 - (incomingSemiMajorAxis / outgoingSemiMajorAxis) * (1.0 - incomingEccentricity);

    // Compute incoming and outgoing velocities at periapsis.
    const double incomingVelocityAtPeriapsis = absoluteIncomingExcessVelocity * std::sqrt((incomingEccentricity + 1.0) / (incomingEccentricity - 1.0));
    const double outgoingVelocityAtPeriapsis = absoluteOutgoingExcessVelocity * std::sqrt((outgoingEccentricity + 1.0) / (outgoingEccentricity - 1.0));

    // Compute necessary delta-V due to velocity-effect.
    velocityEffectDeltaV = std::fabs(outgoingVelocityAtPeriapsis - incomingVelocityAtPeriapsis);

    // SAVE TO POWERED GRAVITY ASSIST PARAMETERS
    saveParameters.incomingSemiMajorAxis = incomingSemiMajorAxis;
    saveParameters.outgoingSemiMajorAxis = outgoingSemiMajorAxis;
    saveParameters.incomingVelocityAtPeriapsis = incomingVelocityAtPeriapsis;
    saveParameters.outgoingVelocityAtPeriapsis = outgoingVelocityAtPeriapsis;
    saveParameters.incomingEccentricity = incomingEccentricity;
    saveParameters.outgoingEccentricity = outgoingEccentricity;
  }
  //
  //                // Here the required maneuver to patch the incoming and outgoing excess velocities is
  //                // calculated. In this implementation, the pericenter radius will be used as iteration
  //                // parameter.
  else {
    // Compute semi-major axis of hyperbolic legs. This is the absolute semi-major axis, because
    // it will otherwisely result in the root of a negative function for various cases during
    // the rootfinding process.
    const double absoluteIncomingSemiMajorAxis = 1.0 * centralBodyGravitationalParameter / absoluteIncomingExcessVelocity / absoluteIncomingExcessVelocity;
    const double absoluteOutgoingSemiMajorAxis = 1.0 * centralBodyGravitationalParameter / absoluteOutgoingExcessVelocity / absoluteOutgoingExcessVelocity;
    // Set the gravity assist function with the variables to perform root finder calculations.
    tudat::mission_segments::PericenterFindingFunctions pericenterFindingFunctions(absoluteIncomingSemiMajorAxis,
                                                                                   absoluteOutgoingSemiMajorAxis,
                                                                                   bendingAngle);

    // Create an object containing the function of which we whish to obtain the root from.
    tudat::basic_mathematics::UnivariateProxyPointer rootFunction = std::make_shared<tudat::basic_mathematics::UnivariateProxy>(
        std::bind(&tudat::mission_segments::PericenterFindingFunctions::computePericenterRadiusFunction,
                  pericenterFindingFunctions, std::placeholders::_1));

    // Add the first derivative of the root function.
    rootFunction->addBinding(-1, std::bind(&tudat::mission_segments::PericenterFindingFunctions::computeFirstDerivativePericenterRadiusFunction, pericenterFindingFunctions, std::placeholders::_1));

    // Set pericenter radius based on result of Newton-Raphson root-finding algorithm.
    const double pericenterRadius = rootFinder->execute(rootFunction,
                                                        smallestPeriapsisDistance);

    // Compute incoming hyperbolic leg eccentricity.
    const double incomingEccentricity = 1.0 + pericenterRadius / absoluteIncomingSemiMajorAxis;

    // Compute outgoing hyperbolic leg eccentricity.
    const double outgoingEccentricity = 1.0 + pericenterRadius / absoluteOutgoingSemiMajorAxis;

    // Compute incoming and outgoing velocities at periapsis.
    const double incomingVelocityAtPeriapsis = absoluteIncomingExcessVelocity * std::sqrt((incomingEccentricity + 1.0) / (incomingEccentricity - 1.0));
    const double outgoingVelocityAtPeriapsis = absoluteOutgoingExcessVelocity * std::sqrt((outgoingEccentricity + 1.0) / (outgoingEccentricity - 1.0));

    // Compute necessary delta-V due to velocity-effect.
    velocityEffectDeltaV = std::fabs(incomingVelocityAtPeriapsis - outgoingVelocityAtPeriapsis);

    // SAVE TO POWERED GRAVITY ASSIST PARAMETERS
    saveParameters.incomingSemiMajorAxis = -absoluteIncomingSemiMajorAxis;
    saveParameters.outgoingSemiMajorAxis = -absoluteOutgoingSemiMajorAxis;
    saveParameters.incomingVelocityAtPeriapsis = incomingVelocityAtPeriapsis;
    saveParameters.outgoingVelocityAtPeriapsis = outgoingVelocityAtPeriapsis;
    saveParameters.incomingEccentricity = incomingEccentricity;
    saveParameters.outgoingEccentricity = outgoingEccentricity;
  }

  // Compute and return the total delta-V.
  saveParameters.totalDeltaV = velocityEffectDeltaV;
  saveParameters.incomingHyperbolicExcessVelocity = incomingHyperbolicExcessVelocity;
  saveParameters.outgoingHyperbolicExcessVelocity = outgoingHyperbolicExcessVelocity;

  return saveParameters;
}

}// namespace prototype
}// namespace tudat