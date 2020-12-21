#ifndef TUDATPY_PLANETARY_RENDEZVOUS_H
#define TUDATPY_PLANETARY_RENDEZVOUS_H

#include <tudat/astro/mission_segments.h>
#include <tudat/basics/basicTypedefs.h>

#include "auxilliary.h"
#include "state.h"
#include <Eigen/Geometry>
#include <tudat/astro/mission_segments.h>

namespace tudat {

namespace prototype {

class PlanetaryRendezvous {
 public:
  PlanetaryRendezvous(
      Eigen::Vector3d incomingVelocity,
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

  double getTurningAngle();

  Eigen::Vector3d getOrbitalPlaneVector();

  Eigen::Vector3d getPeriapsisVelocityUnitVector();

  Eigen::Vector3d getPeriapsisOutgoingVelocity();

  Eigen::Vector3d getPeriapsisIncomingVelocity();

  Eigen::Vector3d getPeriapsisPosition();

  double getSphereOfInfluenceDistanceAtEntrance();

  double getTrueAnomalyEntranceAsymptote();

  double getSemiMajorAxisIncoming();

  double getHyperbolicAnomalyAtEntrance();

  double getTimeOfFlightFromEntrance();

  std::shared_ptr<CartesianState> getOutgoingCartesianState();

  std::shared_ptr<CartesianState> getIncomingCartesianState();

  std::shared_ptr<CartesianState> getCentralBodyState();

  Eigen::Vector3d getIncomingVelocity();

 private:
  // constructor values
  std::shared_ptr<CartesianState> centralBodyCartesianState_;
  Eigen::Vector3d incomingVelocity_;
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
  double turningAngle_;
  double betaAngle_;
  bool turningAngleCalculated_ = false;
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

struct ArrivalParameters {
 public:
  double radiusOfPeriapsis;

  // Keplerian Parameters of Incoming Hyperbolic trajectory.
  double initialVelocityAtPeriapsis;
  double incomingSemiMajorAxis;
  double incomingEccentricity;

  // Keplerian Parameters of Target trajectory.
  double finalVelocityAtPeriapsis;

  // Manoeuvre parameters.
  double totalDeltaV;
  Eigen::Vector3d incomingHyperbolicExcessVelocity;

  // optional
  std::string bodyName;
  std::string dateString;
};

ArrivalParameters calculateArrivalParameters(const double centralBodyGravitationalParameter,
                                             const Eigen::Vector3d &centralBodyVelocityOnEntrySOI,
                                             const Eigen::Vector3d &incomingVelocity,
                                             const double arrivalPeriapsisDistance);

}// namespace prototype
}// namespace tudat

#endif//TUDATPY_PLANETARY_RENDEZVOUS_H
