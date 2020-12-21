#include "planetaryRendezvous.h"
#include "auxilliary.h"
#include "state.h"
#include <Eigen/Geometry>
#include <iostream>
#include <tudat/astro/mission_segments.h>

namespace tudat {

namespace prototype {

PlanetaryRendezvous::PlanetaryRendezvous(
    Eigen::Vector3d incomingVelocity,
    std::shared_ptr<BaseState> centralBodyState,
    double periapsisDistance,
    double gravitationalParameter,
    bool progradeOrbit)
    : incomingVelocity_(incomingVelocity),
      centralBodyCartesianState_(centralBodyState->toCartesian()),
      periapsisDistance_(periapsisDistance),
      gravitationalParameter_(gravitationalParameter),
      progradeOrbit_(progradeOrbit){};

void PlanetaryRendezvous::catchSphereOfInfluenceCalculable() {
  if (centralBodyCartesianState_->getGravitationalParameter() == TUDAT_NAN) {
    throw std::logic_error("Gravitational parameter of the parent body "
                           "(e.g. Sun) must be defined in the central body "
                           "(e.g. Earth) state for sphere of influence calculations.");
  };
}

std::shared_ptr<CartesianState> PlanetaryRendezvous::getCentralBodyState() {
  return centralBodyCartesianState_;
}

Eigen::Vector3d PlanetaryRendezvous::getIncomingVelocity() {
  return incomingVelocity_;
}

Eigen::Vector3d PlanetaryRendezvous::getHyperbolicExcessVelocity() {
  if (!hyperbolicExcessVelocityCalculated_) {
    hyperbolicExcessVelocityCalculated_ = true;
    hyperbolicExcessVelocity_ = incomingVelocity_ - centralBodyCartesianState_->getStateVector().tail<3>();
  }
  return hyperbolicExcessVelocity_;
}

double PlanetaryRendezvous::getHyperbolicExcessSpeed() {
  if (!hyperbolicExcessSpeedCalculated_) {
    hyperbolicExcessSpeedCalculated_ = true;
    hyperbolicExcessSpeed_ = getHyperbolicExcessVelocity().norm();
  }
  return hyperbolicExcessSpeed_;
}

double PlanetaryRendezvous::getEccentricity() {
  if (!eccentricityCalculated_) {
    // curtis equation (8.38)
    eccentricityCalculated_ = true;
    eccentricity_ = 1 + (periapsisDistance_ * pow(getHyperbolicExcessSpeed(), 2)) / (gravitationalParameter_);
  }
  return eccentricity_;
}

double PlanetaryRendezvous::getAngularMomentum() {
  if (!angularMomentumCalculated_) {
    // curtis equation (8.39)
    angularMomentumCalculated_ = true;
    angularMomentum_ = periapsisDistance_ * sqrt(pow(getHyperbolicExcessSpeed(), 2) + (2 * gravitationalParameter_) / (periapsisDistance_));
  }
  return angularMomentum_;
}

double PlanetaryRendezvous::getPeriapsisSpeed() {
  if (!periapsisSpeedCalculated_) {
    // curtis equation (8.40)
    periapsisSpeedCalculated_ = true;
    periapsisSpeed_ = getAngularMomentum() / periapsisDistance_;
  }
  return periapsisSpeed_;
}

double PlanetaryRendezvous::getParkingOrbitSpeed() {
  if (!parkingOrbitSpeedCalculated_) {
    // curtis equation (8.41)
    parkingOrbitSpeedCalculated_ = true;
    parkingOrbitSpeed_ = sqrt(gravitationalParameter_ / periapsisDistance_);
  }
  return parkingOrbitSpeed_;
}

double PlanetaryRendezvous::getRequiredDeltaV() {
  if (!requiredDeltaVCalculated_) {
    // curtis equation (8.41)
    requiredDeltaVCalculated_ = true;
    requiredDeltaV_ = abs(getPeriapsisSpeed() - getParkingOrbitSpeed());
  }
  return requiredDeltaV_;
}

double PlanetaryRendezvous::getTurningAngle() {
  if (!turningAngleCalculated_) {
    // curtis equation (8.54)
    turningAngleCalculated_ = true;
    turningAngle_ = 2 * asin(1 / getEccentricity());
  }
  return turningAngle_;
}

double PlanetaryRendezvous::getBetaAngle() {
  if (!betaAngleCalculated_) {
    // curtis equation (8.41)
    betaAngleCalculated_ = true;
    betaAngle_ = acos(1 / getEccentricity());
  }
  return betaAngle_;
}

Eigen::Vector3d PlanetaryRendezvous::getOrbitalPlaneVector() {
  if (!orbitalPlaneVectorCalculated_) {
    orbitalPlaneVectorCalculated_ = true;
    double _sign;
    if (progradeOrbit_) {
      _sign = -1.0;
    } else {
      _sign = 1.0;
    }
    orbitalPlaneVector_ = _sign * ((centralBodyCartesianState_->getStateVector().head<3>()).cross(getHyperbolicExcessVelocity())).normalized();
  }
  return orbitalPlaneVector_;
}

Eigen::Vector3d PlanetaryRendezvous::getPeriapsisVelocityUnitVector() {
  //
  if (!periapsisVelocityUnitVectorCalculated_) {
    periapsisVelocityUnitVectorCalculated_ = true;
    double _rotationAngle = M_PI / 2 - getBetaAngle();

    // rotate about the orbital plane vector
    periapsisVelocityUnitVector_ = (Eigen::AngleAxisd(
                                        _rotationAngle,
                                        getOrbitalPlaneVector())
                                    * getHyperbolicExcessVelocity())
                                       .normalized();
  }
  return periapsisVelocityUnitVector_;
}

Eigen::Vector3d PlanetaryRendezvous::getPeriapsisOutgoingVelocity() {
  if (!periapsisOutgoingVelocityCalculated_) {
    periapsisOutgoingVelocityCalculated_ = true;
    periapsisOutgoingVelocity_ = getPeriapsisVelocityUnitVector() * getParkingOrbitSpeed();
  }
  return periapsisOutgoingVelocity_;
}

Eigen::Vector3d PlanetaryRendezvous::getPeriapsisIncomingVelocity() {
  if (!periapsisIncomingVelocityCalculated_) {
    periapsisIncomingVelocityCalculated_ = true;
    periapsisIncomingVelocity_ = getPeriapsisVelocityUnitVector() * getPeriapsisSpeed();
  }
  return periapsisIncomingVelocity_;
}

Eigen::Vector3d PlanetaryRendezvous::getPeriapsisPosition() {

  if (!periapsisPositionCalculated_) {
    periapsisPositionCalculated_ = true;
    double rotationAngle = -M_PI / 2;
    periapsisPosition_ = Eigen::AngleAxisd(rotationAngle, getOrbitalPlaneVector()) * getPeriapsisVelocityUnitVector() * periapsisDistance_;
  }
  return periapsisPosition_;
}

double PlanetaryRendezvous::getSphereOfInfluenceDistanceAtEntrance() {
  // check that gravitational parameter is defined for parent body (e.g. Sun)
  catchSphereOfInfluenceCalculable();
  if (!sphereOfInfluenceDistanceAtExitCalculated_) {
    sphereOfInfluenceDistanceAtExitCalculated_ = true;
    sphereOfInfluenceDistanceAtExit_ = calculateSphereOfInfluence(
        centralBodyCartesianState_->getStateVector().segment(0, 3).norm(),
        gravitationalParameter_,
        centralBodyCartesianState_->getGravitationalParameter());
  }
  return sphereOfInfluenceDistanceAtExit_;
}

double PlanetaryRendezvous::getTrueAnomalyEntranceAsymptote() {
  if (!trueAnomalyExitAsymptoteCalculated_) {
    trueAnomalyExitAsymptoteCalculated_ = true;
    double p = pow(getAngularMomentum(), 2) / gravitationalParameter_;
    double r = getSphereOfInfluenceDistanceAtEntrance();
    double ecc = getEccentricity();
    //    trueAnomalyExitAsymptote_ = acos(1 - pow(getAngularMomentum(), 2) / gravitationalParameter_ / getSphereOfInfluenceDistanceAtExit() / getEccentricity());
    trueAnomalyExitAsymptote_ = acos((p / r - 1) / ecc);
    //    trueAnomalyExitAsymptote_ = acos(1 / (-getEccentricity()));
  }
  return trueAnomalyExitAsymptote_;
}

double PlanetaryRendezvous::getSemiMajorAxisIncoming() {
  if (!semiMajorAxisOutgoingCalculated_) {
    semiMajorAxisOutgoingCalculated_ = true;
    semiMajorAxisOutgoing_ = -gravitationalParameter_ / pow(getHyperbolicExcessSpeed(), 2);
  }
  return semiMajorAxisOutgoing_;
}

double PlanetaryRendezvous::getHyperbolicAnomalyAtEntrance() {
  if (!hyperbolicAnomalyAtExitCalculated_) {
    hyperbolicAnomalyAtExitCalculated_ = true;
    hyperbolicAnomalyAtExit_ = acosh(
        (getEccentricity() + cos(getTrueAnomalyEntranceAsymptote()))
        / (1 + getEccentricity() * cos(getTrueAnomalyEntranceAsymptote())));
  }
  return hyperbolicAnomalyAtExit_;
}

double PlanetaryRendezvous::getTimeOfFlightFromEntrance() {
  if (!timeOfFlightToExitCalculated_) {
    timeOfFlightToExitCalculated_ = true;
    double meanAnomaly = getEccentricity() * sinh(getHyperbolicAnomalyAtEntrance()) - getHyperbolicAnomalyAtEntrance();
    //    std::cout << "--" << meanAnomaly;
    timeOfFlightToExit_ = -sqrt(
                              pow(-getSemiMajorAxisIncoming(), 3) / gravitationalParameter_)
        * meanAnomaly;
    //        abs(
    //            sqrt(pow((getSemiMajorAxisOutgoing()), 3) / gravitationalParameter_)
    //            * (getEccentricity() * sinh(getHyperbolicAnomalyAtExit()) - getHyperbolicAnomalyAtExit()));
  }
  return timeOfFlightToExit_;
}

std::shared_ptr<CartesianState> PlanetaryRendezvous::getOutgoingCartesianState() {
  return std::make_shared<CartesianState>(
      getPeriapsisPosition(),
      getPeriapsisOutgoingVelocity(),
      gravitationalParameter_
      );
};

std::shared_ptr<CartesianState> PlanetaryRendezvous::getIncomingCartesianState() {
  return std::make_shared<CartesianState>(
      getPeriapsisPosition(),
      getPeriapsisIncomingVelocity(),
      gravitationalParameter_);
};

ArrivalParameters calculateArrivalParameters(const double centralBodyGravitationalParameter,
                                             const Eigen::Vector3d &centralBodyVelocityOnEntrySOI,
                                             const Eigen::Vector3d &incomingVelocity,
                                             const double arrivalPeriapsisDistance) {

  ArrivalParameters arrivalParameters;
  const Eigen::Vector3d incomingHyperbolicExcessVelocity = incomingVelocity - centralBodyVelocityOnEntrySOI;

  // compute incoming hyperbolic leg
  const double incomingExcessVelocityNorm = incomingHyperbolicExcessVelocity.norm();
  const double incomingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter / incomingExcessVelocityNorm / incomingExcessVelocityNorm;
  const double incomingEccentricity = 1 - arrivalPeriapsisDistance / incomingSemiMajorAxis;
  const double incomingVelocityAtPeriapsis = incomingExcessVelocityNorm * std::sqrt((incomingEccentricity + 1.0) / (incomingEccentricity - 1.0));

  // compute outgoing circular orbit
  const double outgoingVelocityAtPeriapsis = sqrt(centralBodyGravitationalParameter / arrivalPeriapsisDistance);

  // Compute necessary delta-V due to velocity-effect.
  const double poweredDeltaV = std::fabs(outgoingVelocityAtPeriapsis - incomingVelocityAtPeriapsis);

  // SAVE TO POWERED GRAVITY ASSIST PARAMETERS
  arrivalParameters.incomingSemiMajorAxis = incomingSemiMajorAxis;
  arrivalParameters.radiusOfPeriapsis = arrivalPeriapsisDistance;
  arrivalParameters.finalVelocityAtPeriapsis = outgoingVelocityAtPeriapsis;
  arrivalParameters.initialVelocityAtPeriapsis = incomingVelocityAtPeriapsis;
  arrivalParameters.incomingEccentricity = incomingEccentricity;

  // Compute and return the total delta-V.
  arrivalParameters.totalDeltaV = poweredDeltaV;
  arrivalParameters.incomingHyperbolicExcessVelocity = incomingHyperbolicExcessVelocity;

  return arrivalParameters;
}

}// namespace prototype
}// namespace tudat