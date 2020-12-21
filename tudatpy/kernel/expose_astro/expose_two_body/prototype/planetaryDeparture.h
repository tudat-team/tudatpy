#ifndef TUDATPY_PLANETARY_DEPARTURE_H
#define TUDATPY_PLANETARY_DEPARTURE_H

#include <tudat/astro/mission_segments.h>
#include <tudat/basics/basicTypedefs.h>

#include "auxilliary.h"
#include "state.h"
#include <Eigen/Geometry>
#include <tudat/astro/mission_segments.h>

namespace tms = tudat::mission_segments;

namespace tudat {

namespace prototype {

class PlanetaryDeparture {
 public:
  PlanetaryDeparture(
      Eigen::Vector3d outgoingVelocity,
      std::shared_ptr<BaseState> centralBodyState,
      double periapsisDistance,
      double gravitationalParameter,
      bool progradeOrbit = true);

  void catchSphereOfInfluenceCalculable();

  Eigen::Vector3d getHyperbolicExcessVelocity();

  double getHyperbolicExcessSpeed();

  double getEccentricity();

  double getAngularMomentum();

  double getPeriapsisSpeed();

  double getParkingOrbitSpeed();

  double getRequiredDeltaV();

  double getBetaAngle();

  Eigen::Vector3d getOrbitalPlaneVector();

  Eigen::Vector3d getPeriapsisVelocityUnitVector();

  Eigen::Vector3d getPeriapsisOutgoingVelocity();

  Eigen::Vector3d getPeriapsisIncomingVelocity();

  Eigen::Vector3d getPeriapsisPosition();

  double getSphereOfInfluenceDistanceAtExit();

  double getTrueAnomalyExitAsymptote();

  double getSemiMajorAxisOutgoing();

  double getHyperbolicAnomalyAtExit();

  double getTimeOfFlightToExit();

  std::shared_ptr<CartesianState> getOutgoingCartesianState();

  std::shared_ptr<CartesianState> getIncomingCartesianState();

  std::shared_ptr<CartesianState> getCentralBodyState();

  Eigen::Vector3d getOutgoingVelocity();


 private:
  // constructor values
  std::shared_ptr<CartesianState> centralBodyCartesianState_;
  Eigen::Vector3d outgoingVelocity_;
  double gravitationalParameter_;
  double periapsisDistance_;
  bool progradeOrbit_;

  // derived values
  Eigen::Vector3d hyperbolicExcessVelocity_;
  bool hyperbolicExcessVelocityCalculated_ = false;
  double hyperbolicExcessSpeed_;
  bool hyperbolicExcessSpeedCalculated_ = false;
  double eccentricity_;
  bool eccentricityCalculated_ = false;
  double angularMomentum_;
  bool angularMomentumCalculated_ = false;
  double periapsisSpeed_;
  bool periapsisSpeedCalculated_ = false;
  double parkingOrbitSpeed_;
  bool parkingOrbitSpeedCalculated_ = false;
  double requiredDeltaV_;
  bool requiredDeltaVCalculated_ = false;
  double betaAngle_;
  bool betaAngleCalculated_ = false;
  Eigen::Vector3d orbitalPlaneVector_;
  bool orbitalPlaneVectorCalculated_ = false;
  Eigen::Vector3d periapsisVelocityUnitVector_;
  bool periapsisVelocityUnitVectorCalculated_ = false;
  Eigen::Vector3d periapsisOutgoingVelocity_;
  bool periapsisOutgoingVelocityCalculated_ = false;
  Eigen::Vector3d periapsisIncomingVelocity_;
  bool periapsisIncomingVelocityCalculated_ = false;
  Eigen::Vector3d periapsisPosition_;
  bool periapsisPositionCalculated_ = false;
  double sphereOfInfluenceDistanceAtExit_;
  bool sphereOfInfluenceDistanceAtExitCalculated_ = false;
  double trueAnomalyExitAsymptote_;
  bool trueAnomalyExitAsymptoteCalculated_ = false;
  double semiMajorAxisOutgoing_;
  bool semiMajorAxisOutgoingCalculated_ = false;
  double hyperbolicAnomalyAtExit_;
  bool hyperbolicAnomalyAtExitCalculated_ = false;
  double timeOfFlightToExit_;
  bool timeOfFlightToExitCalculated_ = false;
};

struct DepartureParameters {
 public:
  double radiusOfPeriapsis;

  // Keplerian Parameters of Incoming Hyperbolic trajectory.
  double initialVelocityAtPeriapsis;

  // Keplerian Parameters of Outgoing Hyperbolic trajectory.
  double outgoingSemiMajorAxis;
  double outgoingEccentricity;
  double finalVelocityAtPeriapsis;

  // Manoeuvre parameters.
  double totalDeltaV;
  Eigen::Vector3d outgoingHyperbolicExcessVelocity;

  // optional
  std::string bodyName;
  std::string dateString;
};

DepartureParameters calculateDepartureParameters(const double centralBodyGravitationalParameter,
                                                 const Eigen::Vector3d &centralBodyVelocityOnExitSOI,
                                                 const Eigen::Vector3d &outgoingVelocity,
                                                 const double departurePeriapsisDistance);

}// namespace prototype
}// namespace tudat

#endif//TUDATPY_PLANETARY_DEPARTURE_H