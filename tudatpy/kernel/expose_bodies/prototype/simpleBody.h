/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_SIMPLEBODY_H
#define TUDATPY_SIMPLEBODY_H

#include <tudat/astro/ephemerides.h>




#include "../../expose_astro/prototype/frames.h"
#include <memory>
#include <string>
#include <tudat/astro/ephemerides.h>
#include <vector>

namespace tudat {

using namespace ephemerides;

namespace bodies {

class SimpleBody : public std::enable_shared_from_this<SimpleBody> {
 public:
  SimpleBody();

//  SimpleBody(
//      std::shared_ptr<SimpleBody> parent,
//      std::shared_ptr<Ephemeris> ephemeris,
//      double gravitationalParameter,
//      std::string name,
//      std::string symbol = "",
//      double radius = 0,
//      double polarRadius = 0,
//      double meanRadius = 0,
//      double rotationalPeriod = 0,
//      double J2 = 0,
//      double J3 = 0,
//      double mass = 0,
//      double atmosphereDistance = 0,
//      double flybyLimit = 0);

  SimpleBody(
      std::shared_ptr<SimpleBody> parent,
      double gravitationalParameter,
      std::string name,
      std::string symbol = "",
      double radius = 0,
      double polarRadius = 0,
      double meanRadius = 0,
      double rotationalPeriod = 0,
      double J2 = 0,
      double J3 = 0,
      double mass = 0,
      double atmosphereDistance = 0,
      double flybyLimit = 0,
      double sphereOfInfluenceDistance = 0);

  SimpleBody(
      double gravitationalParameter,
      std::string name,
      std::string symbol = "",
      double radius = 0,
      double polarRadius = 0,
      double meanRadius = 0,
      double rotationalPeriod = 0,
      double J2 = 0,
      double J3 = 0,
      double mass = 0,
      double atmosphereDistance = 0,
      double flybyLimit = 0);

  bool isSet();
  void addChild(std::shared_ptr<SimpleBody> body);
  std::vector<std::shared_ptr<SimpleBody>> getChildren();
  std::shared_ptr<SimpleBody> getParent();

  double getMeanRadius();
  double getPolarRadius();
  double getRadius();
  double getRotationalPeriod();
  double getJ2();
  double getJ3();
  double getMass();
  std::string getName();
  std::string getSymbol();
  double getGravitationalParameter();
  double getAtmosphereDistance();
  double getFlybyLimit();

  // methods for complying with existing tudat
  void setEphemeris(std::shared_ptr<Ephemeris> ephemeris);

 private:
  bool isSet_;
  std::shared_ptr<SimpleBody> parent_;
  double gravitationalParameter_;
  std::string name_;
  std::string symbol_;
  double radius_;
  double polarRadius_;
  double meanRadius_;
  double rotationalPeriod_;
  double J2_;
  double J3_;
  double mass_;
  double flybyLimit_;
  double atmosphereDistance_;
  std::vector<std::shared_ptr<SimpleBody>> children_;
};

class SimpleSystemOfBodies {
 public:
  SimpleSystemOfBodies(
      std::shared_ptr<SimpleBody> parent,
      std::shared_ptr<ReferenceFrame> referenceFrame,
      std::vector<std::shared_ptr<SimpleBody>> children);

  double getGravitationalParameter();
//  std::shared_ptr<Ephemeris> getEphemeris();

  std::map<std::string, std::shared_ptr<SimpleBody>> getBodyMap();
  std::shared_ptr<ReferenceFrame> getReferenceFrame();
  std::shared_ptr<SimpleBody> getParent();


 private:
  std::shared_ptr<SimpleBody> parent_;
  std::shared_ptr<ReferenceFrame> referenceFrame_;
  std::vector<std::shared_ptr<SimpleBody>> children_;
  std::map<std::string, std::shared_ptr<SimpleBody>> bodyMap_;
};

}// namespace bodies
}// namespace tudat
#endif//TUDATPY_SIMPLEBODY_H
