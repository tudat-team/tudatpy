#ifndef TUDATPY_PLANETARY_FLYBY_H
#define TUDATPY_PLANETARY_FLYBY_H

#include <tudat/astro/mission_segments.h>
#include <tudat/basics/basicTypedefs.h>
#include <tudat/math/root_finders.h>

using namespace tudat::root_finders;
namespace trf = tudat::root_finders;

namespace tudat {

namespace prototype {

struct GravityAssistParameters {
 public:
  double radiusOfPeriapsis;

  // Keplerian Parameters of Incoming Hyperbolic trajectory.
  double incomingSemiMajorAxis;
  double incomingEccentricity;
  double incomingVelocityAtPeriapsis;

  // Keplerian Parameters of Outgoing Hyperbolic trajectory.
  double outgoingSemiMajorAxis;
  double outgoingEccentricity;
  double outgoingVelocityAtPeriapsis;

  // Manoeuvre parameters.
  double totalDeltaV;
  Eigen::Vector3d incomingHyperbolicExcessVelocity;
  Eigen::Vector3d outgoingHyperbolicExcessVelocity;

  // optional
  std::string bodyName;
  std::string dateString;
};

GravityAssistParameters calculateGravityAssistParameters(
    const double centralBodyGravitationalParameter,
    const Eigen::Vector3d &centralBodyVelocityOnEntrySOI,
    const Eigen::Vector3d &centralBodyVelocityOnExitSOI,
    const Eigen::Vector3d &incomingVelocity,
    const Eigen::Vector3d &outgoingVelocity,
    const double smallestPeriapsisDistance,
    const bool useEccentricityInsteadOfPericenter = true,
    const double speedTolerance = 1e-6,
    RootFinderPointer rootFinder = std::make_shared<trf::NewtonRaphson>(1.0e-12, 1000));

}// namespace prototype
}// namespace tudat

#endif//TUDATPY_PLANETARY_FLYBY_H