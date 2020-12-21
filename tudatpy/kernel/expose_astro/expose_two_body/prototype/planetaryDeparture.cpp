
#include "planetaryDeparture.h"
#include "auxilliary.h"
#include "state.h"
#include <Eigen/Geometry>
#include <iostream>
#include <tudat/astro/mission_segments.h>

namespace tms = tudat::mission_segments;

namespace tudat {

namespace prototype {

PlanetaryDeparture::PlanetaryDeparture(
    Eigen::Vector3d outgoingVelocity,
    std::shared_ptr<BaseState> centralBodyState,
    double periapsisDistance,
    double gravitationalParameter,
    bool progradeOrbit)
    : outgoingVelocity_(outgoingVelocity),
      centralBodyCartesianState_(centralBodyState->toCartesian()),
      periapsisDistance_(periapsisDistance),
      gravitationalParameter_(gravitationalParameter),
      progradeOrbit_(progradeOrbit){};

void PlanetaryDeparture::catchSphereOfInfluenceCalculable() {
  if (centralBodyCartesianState_->getGravitationalParameter() == TUDAT_NAN) {
    throw std::logic_error("Gravitational parameter of the parent body "
                           "(e.g. Sun) must be defined in the central body "
                           "(e.g. Earth) state for sphere of influence calculations.");
  };
}

std::shared_ptr<CartesianState> PlanetaryDeparture::getCentralBodyState() {
  return centralBodyCartesianState_;
}

Eigen::Vector3d PlanetaryDeparture::getOutgoingVelocity() {
  return outgoingVelocity_;
}

Eigen::Vector3d PlanetaryDeparture::getHyperbolicExcessVelocity() {
  if (!hyperbolicExcessVelocityCalculated_) {
    hyperbolicExcessVelocityCalculated_ = true;
    hyperbolicExcessVelocity_ = outgoingVelocity_ - centralBodyCartesianState_->getStateVector().tail<3>();
  }
  // [bug departure] verified
  return hyperbolicExcessVelocity_;
}

double PlanetaryDeparture::getHyperbolicExcessSpeed() {
  if (!hyperbolicExcessSpeedCalculated_) {
    hyperbolicExcessSpeedCalculated_ = true;
    hyperbolicExcessSpeed_ = getHyperbolicExcessVelocity().norm();
  }
  // [bug departure] verified
  return hyperbolicExcessSpeed_;
}

double PlanetaryDeparture::getEccentricity() {
  if (!eccentricityCalculated_) {
    // curtis equation (8.38)
    eccentricityCalculated_ = true;
    eccentricity_ = 1 + (periapsisDistance_ * pow(getHyperbolicExcessSpeed(), 2)) / (gravitationalParameter_);
  }
  // [bug departure] verified
  return eccentricity_;
}

double PlanetaryDeparture::getAngularMomentum() {
  if (!angularMomentumCalculated_) {
    // curtis equation (8.39)
    angularMomentumCalculated_ = true;
    angularMomentum_ = periapsisDistance_ * sqrt(pow(getHyperbolicExcessSpeed(), 2) + (2 * gravitationalParameter_) / (periapsisDistance_));
  }
  // [bug departure] verified
  return angularMomentum_;
}

double PlanetaryDeparture::getPeriapsisSpeed() {
  if (!periapsisSpeedCalculated_) {
    // curtis equation (8.40)
    periapsisSpeedCalculated_ = true;
    periapsisSpeed_ = getAngularMomentum() / periapsisDistance_;
  }
  return periapsisSpeed_;
}

double PlanetaryDeparture::getParkingOrbitSpeed() {
  if (!parkingOrbitSpeedCalculated_) {
    // curtis equation (8.41)
    parkingOrbitSpeedCalculated_ = true;
    parkingOrbitSpeed_ = sqrt(gravitationalParameter_ / periapsisDistance_);
  }
  return parkingOrbitSpeed_;
}

double PlanetaryDeparture::getRequiredDeltaV() {
  if (!requiredDeltaVCalculated_) {
    // curtis equation (8.41)
    requiredDeltaVCalculated_ = true;
    requiredDeltaV_ = abs(getPeriapsisSpeed() - getParkingOrbitSpeed());
  }
  return requiredDeltaV_;
}

double PlanetaryDeparture::getBetaAngle() {
  if (!betaAngleCalculated_) {
    // curtis equation (8.41)
    betaAngleCalculated_ = true;
    betaAngle_ = acos(1 / getEccentricity());
  }
  return betaAngle_;
}

Eigen::Vector3d PlanetaryDeparture::getOrbitalPlaneVector() {
  if (!orbitalPlaneVectorCalculated_) {
    orbitalPlaneVectorCalculated_ = true;
    double _sign;
    if (progradeOrbit_) {
      _sign = 1.0;
    } else {
      _sign = -1.0;
    }
    orbitalPlaneVector_ = _sign * ((centralBodyCartesianState_->getStateVector().head<3>()).cross(getHyperbolicExcessVelocity())).normalized();
  }
  return orbitalPlaneVector_;
}

Eigen::Vector3d PlanetaryDeparture::getPeriapsisPosition() {

  if (!periapsisPositionCalculated_) {
    periapsisPositionCalculated_ = true;
    double rotationAngle = getBetaAngle();
    periapsisPosition_ = -(Eigen::AngleAxisd(rotationAngle, getOrbitalPlaneVector()) * getHyperbolicExcessVelocity()).normalized() * periapsisDistance_;
  }
  return periapsisPosition_;
}

Eigen::Vector3d PlanetaryDeparture::getPeriapsisVelocityUnitVector() {
  //
  if (!periapsisVelocityUnitVectorCalculated_) {
    periapsisVelocityUnitVectorCalculated_ = true;
    double _rotationAngle = M_PI / 2;

    // rotate about the orbital plane vector
    periapsisVelocityUnitVector_ = (Eigen::AngleAxisd(
                                         _rotationAngle,
                                         getOrbitalPlaneVector())
                                     * getPeriapsisPosition())
                                        .normalized();
  }
  return periapsisVelocityUnitVector_;
}

Eigen::Vector3d PlanetaryDeparture::getPeriapsisOutgoingVelocity() {
  if (!periapsisOutgoingVelocityCalculated_) {
    periapsisOutgoingVelocityCalculated_ = true;
    periapsisOutgoingVelocity_ = getPeriapsisVelocityUnitVector() * getPeriapsisSpeed();
  }
  return periapsisOutgoingVelocity_;
}

Eigen::Vector3d PlanetaryDeparture::getPeriapsisIncomingVelocity() {
  if (!periapsisIncomingVelocityCalculated_) {
    periapsisIncomingVelocityCalculated_ = true;
    periapsisIncomingVelocity_ = getPeriapsisVelocityUnitVector() * getParkingOrbitSpeed();
  }
  return periapsisIncomingVelocity_;
}

double PlanetaryDeparture::getSphereOfInfluenceDistanceAtExit() {
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

double PlanetaryDeparture::getTrueAnomalyExitAsymptote() {
  if (!trueAnomalyExitAsymptoteCalculated_) {
    trueAnomalyExitAsymptoteCalculated_ = true;
    double p = pow(getAngularMomentum(), 2) / gravitationalParameter_;
    double r = getSphereOfInfluenceDistanceAtExit();
    double ecc = getEccentricity();
    //    trueAnomalyExitAsymptote_ = acos(1 - pow(getAngularMomentum(), 2) / gravitationalParameter_ / getSphereOfInfluenceDistanceAtExit() / getEccentricity());
    trueAnomalyExitAsymptote_ = acos((p / r - 1) / ecc);
    //    trueAnomalyExitAsymptote_ = acos(1 / (-getEccentricity()));
  }
  return trueAnomalyExitAsymptote_;
}

double PlanetaryDeparture::getSemiMajorAxisOutgoing() {
  if (!semiMajorAxisOutgoingCalculated_) {
    semiMajorAxisOutgoingCalculated_ = true;
    semiMajorAxisOutgoing_ = -gravitationalParameter_ / pow(getHyperbolicExcessSpeed(), 2);
  }
  return semiMajorAxisOutgoing_;
}

double PlanetaryDeparture::getHyperbolicAnomalyAtExit() {
  if (!hyperbolicAnomalyAtExitCalculated_) {
    hyperbolicAnomalyAtExitCalculated_ = true;
    hyperbolicAnomalyAtExit_ = acosh(
        (getEccentricity() + cos(getTrueAnomalyExitAsymptote()))
        / (1 + getEccentricity() * cos(getTrueAnomalyExitAsymptote())));
  }
  return hyperbolicAnomalyAtExit_;
}

double PlanetaryDeparture::getTimeOfFlightToExit() {
  if (!timeOfFlightToExitCalculated_) {
    timeOfFlightToExitCalculated_ = true;
    double meanAnomaly = getEccentricity() * sinh(getHyperbolicAnomalyAtExit()) - getHyperbolicAnomalyAtExit();
    //    std::cout << "--" << meanAnomaly;
    timeOfFlightToExit_ = sqrt(
                              pow(-getSemiMajorAxisOutgoing(), 3) / gravitationalParameter_)
        * meanAnomaly;
    //        abs(
    //            sqrt(pow((getSemiMajorAxisOutgoing()), 3) / gravitationalParameter_)
    //            * (getEccentricity() * sinh(getHyperbolicAnomalyAtExit()) - getHyperbolicAnomalyAtExit()));
  }
  return timeOfFlightToExit_;
}

std::shared_ptr<CartesianState> PlanetaryDeparture::getOutgoingCartesianState() {
  return std::make_shared<CartesianState>(
      getPeriapsisPosition(),
      getPeriapsisOutgoingVelocity(),
      gravitationalParameter_);
};

std::shared_ptr<CartesianState> PlanetaryDeparture::getIncomingCartesianState() {
  return std::make_shared<CartesianState>(
      getPeriapsisPosition(),
      getPeriapsisIncomingVelocity(),
      gravitationalParameter_);
};

DepartureParameters calculateDepartureParameters(const double centralBodyGravitationalParameter,
                                                 const Eigen::Vector3d &centralBodyVelocityOnExitSOI,
                                                 const Eigen::Vector3d &outgoingVelocity,
                                                 const double departurePeriapsisDistance) {

  DepartureParameters saveParameters;
  const Eigen::Vector3d outgoingHyperbolicExcessVelocity = outgoingVelocity - centralBodyVelocityOnExitSOI;
  const double absoluteOutgoingExcessVelocity = outgoingHyperbolicExcessVelocity.norm();
  const double outgoingSemiMajorAxis = -1.0 * centralBodyGravitationalParameter / absoluteOutgoingExcessVelocity / absoluteOutgoingExcessVelocity;

  // Compute outgoing hyperbolic leg eccentricity.
  const double outgoingEccentricity = 1 - departurePeriapsisDistance / outgoingSemiMajorAxis;
  const double outgoingVelocityAtPeriapsis = absoluteOutgoingExcessVelocity * std::sqrt((outgoingEccentricity + 1.0) / (outgoingEccentricity - 1.0));
  const double initialVelocityAtPeriapsis = sqrt(
      centralBodyGravitationalParameter / departurePeriapsisDistance);

  // Compute necessary delta-V due to velocity-effect.
  const double velocityEffectDeltaV = std::fabs(initialVelocityAtPeriapsis - outgoingVelocityAtPeriapsis);

  // SAVE TO POWERED GRAVITY ASSIST PARAMETERS
  saveParameters.outgoingSemiMajorAxis = outgoingSemiMajorAxis;
  saveParameters.finalVelocityAtPeriapsis = outgoingVelocityAtPeriapsis;
  saveParameters.initialVelocityAtPeriapsis = initialVelocityAtPeriapsis;
  saveParameters.outgoingEccentricity = outgoingEccentricity;
  saveParameters.radiusOfPeriapsis = departurePeriapsisDistance;

  // Compute and return the total delta-V.
  saveParameters.totalDeltaV = velocityEffectDeltaV;
  saveParameters.outgoingHyperbolicExcessVelocity = outgoingHyperbolicExcessVelocity;

  return saveParameters;
}

}// namespace prototype
}// namespace tudat