/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "simpleBody.h"

namespace tudat {

namespace bodies {

SimpleBody::SimpleBody() : isSet_(false) {}

SimpleBody::SimpleBody(
    std::shared_ptr<SimpleBody> parent,
    std::shared_ptr<Ephemeris> ephemeris,
    double gravitationalParameter,
    std::string name,
    std::string symbol,
    double radius,
    double polarRadius,
    double meanRadius,
    double rotationalPeriod,
    double J2,
    double J3,
    double mass,
    double atmosphereDistance,
    double flybyLimit)
    : isSet_(true),
      parent_(parent),
      gravitationalParameter_(gravitationalParameter),
      name_(name),
      symbol_(symbol),
      radius_(radius),
      polarRadius_(polarRadius),
      meanRadius_(meanRadius),
      rotationalPeriod_(rotationalPeriod),
      J2_(J2),
      J3_(J3),
      mass_(mass),
      atmosphereDistance_(atmosphereDistance),
      flybyLimit_(flybyLimit){};

SimpleBody::SimpleBody(
    double gravitationalParameter,
    std::string name,
    std::string symbol,
    double radius,
    double polarRadius,
    double meanRadius,
    double rotationalPeriod,
    double J2,
    double J3,
    double mass,
    double atmosphereDistance,
    double flybyLimit)
    : isSet_(true),
      gravitationalParameter_(gravitationalParameter),
      name_(name),
      symbol_(symbol),
      radius_(radius),
      polarRadius_(polarRadius),
      meanRadius_(meanRadius),
      rotationalPeriod_(rotationalPeriod),
      J2_(J2),
      J3_(J3),
      mass_(mass),
      atmosphereDistance_(atmosphereDistance),
      flybyLimit_(flybyLimit){};

bool SimpleBody::isSet() { return isSet_; };
void SimpleBody::addChild(std::shared_ptr<SimpleBody> body) { children_.push_back(body); };
std::vector<std::shared_ptr<SimpleBody>> SimpleBody::getChildren() { return children_; };
std::shared_ptr<SimpleBody> SimpleBody::getParent() { return parent_; };
double SimpleBody::getMeanRadius() { return meanRadius_; };
double SimpleBody::getPolarRadius() { return polarRadius_; };
double SimpleBody::getRadius() { return radius_; };
double SimpleBody::getRotationalPeriod() { return rotationalPeriod_; };
double SimpleBody::getJ2() { return J2_; };
double SimpleBody::getJ3() { return J3_; };
double SimpleBody::getMass() { return mass_; };
std::string SimpleBody::getName() { return name_; };
std::string SimpleBody::getSymbol() { return symbol_; };
double SimpleBody::getGravitationalParameter() { return gravitationalParameter_; };
double SimpleBody::getAtmosphereDistance() { return atmosphereDistance_; };
double SimpleBody::getFlybyLimit() { return flybyLimit_; };

SimpleSystemOfBodies::SimpleSystemOfBodies(
    std::shared_ptr<SimpleBody> parent,
    std::shared_ptr<ReferenceFrame> referenceFrame,
    std::vector<std::shared_ptr<SimpleBody>> children)
    : parent_(parent), referenceFrame_(referenceFrame), children_(children) {
  bodyMap_.insert({parent->getName(), parent});
  for (auto child : children) {
    parent->addChild(child);
    bodyMap_.insert({child->getName(), child});
  }
}

std::map<std::string, std::shared_ptr<SimpleBody>> SimpleSystemOfBodies::getBodyMap() { return bodyMap_; };
std::shared_ptr<ReferenceFrame> SimpleSystemOfBodies::getReferenceFrame() { return referenceFrame_; };
std::shared_ptr<SimpleBody> SimpleSystemOfBodies::getParent() { return parent_; };

double SimpleSystemOfBodies::getGravitationalParameter() {
  double gravitationalParameterSum = parent_->getGravitationalParameter();
  for (auto child : children_) {
    gravitationalParameterSum += child->getGravitationalParameter();
  }
  return gravitationalParameterSum;
}

};// namespace bodies

}// namespace tudat
