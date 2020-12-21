/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "orbit.h"
#include <memory>
using namespace tudat::bodies;
using namespace tudat::prototype;

Orbit::Orbit() : isSet_(false){};// default constructor used for default arguments using Orbit.

Orbit::Orbit(
    std::shared_ptr<SimpleBody> body,
    std::shared_ptr<BaseState> state,
    double epoch,
    std::shared_ptr<ReferenceFrame> referenceFrame)
    : state_(state), epoch_(epoch), referenceFrame_(referenceFrame), isSet_(true){};

//
bool Orbit::isSet() { return isSet_; }

// static class method
std::shared_ptr<Orbit> Orbit::fromSpice(double epoch, std::string body_name){
    std::shared_ptr<BaseState> state = BaseState::fromSpice(std::string bodyName,
      epoch, // ephemerisTime
      double gravitationalParameter,
      ReferenceFrame referenceFrame,
      std::string aberrationCorrections = "none");

    return std::make_shared<Orbit>()};

// static class method
std::shared_ptr<Orbit> Orbit::fromBody(){

};

//// static class method
//std::shared_ptr<Orbit> Orbit::fromEphemeris(std::shared_ptr<Ephemeris> ephemeris){
//
//};

// static class method
std::shared_ptr<Orbit> Orbit::fromState(std::shared_ptr<BaseState> state, double epoch, std::shared_ptr<SimpleBody> body);

// from state object
//void Orbit::propagate(double deltaTransient, std::string transient = "time", std::string method = "keplerian");

// from state object / cartesian
Eigen::Vector3d Orbit::getPosition() { return state_->toCartesian()->getPositionVector(); };
Eigen::Vector3d Orbit::getVelocity() { return state_->toCartesian()->getVelocityVector(); };

// from state object / keplerian
double Orbit::getSemiMajorAxis() { return state_->toKeplerian()->getSemiMajorAxis(); };
double Orbit::getEccentricity() { return state_->toKeplerian()->getEccentricity(); };
double Orbit::getRightAscension() { return state_->toKeplerian()->getLongitudeAscendingNode(); };
double Orbit::getArgOfPeriapsis() { return state_->toKeplerian()->getArgumentOfPeriapsis(); };
double Orbit::getInclination() { return state_->toKeplerian()->getInclination(); };
double Orbit::getTrueAnomaly() { return state_->toKeplerian()->getTrueAnomaly(); };

// from state object
Eigen::Vector6d Orbit::getStateVector(std::string which) {
  Eigen::Vector6d stateVector;
  if (which == "cartesian") {
    stateVector = state_->toCartesian()->getStateVector();
  } else if (which == "keplerian") {
    stateVector = state_->toKeplerian()->getStateVector();
  } else if (which == "equinoctial") {
    stateVector = state_->toEquinoctial()->getStateVector();
  }
  return stateVector;
};

bool Orbit::isOpen() {
  return (getEccentricity() >= 1);
};

Eigen::Vector6d Orbit::sample(
    double initial_true_anomaly,
    double final_true_anomaly,
    std::shared_ptr<ReferenceFrame> referenceFrame) {
  Eigen::Vector6d sampled;
  if (referenceFrame->isSet()) {

  } else {
  }
  if (isOpen()) {

  } else {
  }
  return sampled;
};
