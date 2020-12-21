/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_ORBIT_H
#define TUDATPY_ORBIT_H

#include "../../../expose_bodies/prototype/simpleBody.h"
#include "../../prototype/frames.h"
#include <tudat/astro/ephemerides.h>
#include "state.h"

using namespace tudat::bodies;
using namespace tudat::prototype;

class Orbit {
  // An orbit is a state, around a body, at a given epoch.

  Orbit();// default constructor used for default arguments using Orbit.

  Orbit(std::shared_ptr<SimpleBody> body,
        std::shared_ptr<BaseState> state,
        double epoch,
        std::shared_ptr<ReferenceFrame> referenceFrame);

  // static classmethod
  static std::shared_ptr<Orbit> fromSpice(
      std::shared_ptr<SimpleBody> body
      double epoch);

  // static classmethod
  static std::shared_ptr<Orbit> fromEphemeris(
      std::shared_ptr<SimpleBody> body,
      std::shared_ptr<Ephemeris> ephemeris
      double epoch);

  // static classmethod
  static std::shared_ptr<Orbit> fromState(
      std::shared_ptr<SimpleBody> body,
      std::shared_ptr<BaseState> state,
      double epoch);

  std::shared_ptr<Orbit> applyDeltaV
      orbit.sample()
      orbit_plotter = OrbitPlotter()
      for deltav in [100, 200, 300]:
        new_orbit = orbit.propagate(100).apply_deltaV(deltaV)
        orbit_plotter.plot(new_orbit, label=f"deltav = {}")

  orbit_1 = depature_orbit.propagate(t_ougoing).update_attractor(Sun).propagate(t_lambert).change_attractor(arrival_body).propagate(t_incoming)
  PlanetaryDeparture()
  PlanetaryRendezvous()
  PlanetaryFlyby
  .orbit[0].apply_maneuvre()
  PlanetaryFlyby.orbit[1]

  // status of orbit
  bool isSet();

  // from state object
  //  void propagate(double delta_transient, std::string transient = "time", std::string method = "keplerian");

  // from state object / cartesian
  Eigen::Vector3d getPosition();
  Eigen::Vector3d getVelocity();

  // from state object / keplerian
  double getSemiMajorAxis();
  double getEccentricity();
  double getRightAscension();
  double getArgOfPeriapsis();
  double getInclination();
  double getTrueAnomaly();

  // from state object / equinoctial
  Eigen::Vector6d getStateVector(std::string which = "cartesian");

  bool isOpen();

  Eigen::Vector6d sample(
      double initial_true_anomaly,
      double final_true_anomaly,
      std::shared_ptr<ReferenceFrame> referenceFrame = std::make_shared<ReferenceFrame>());

 private:
  // variables defining orbit
  bool isSet_;
  double epoch_;
  std::shared_ptr<BaseState> state_;
  std::shared_ptr<ReferenceFrame> referenceFrame_;

  // internal for caching return variables
  //  std::unordered_map<std::string, double> doubleVariableCache_;
  //  std::unordered_map<std::string, Eigen::Vector3d> vector3dVariableCache_;
  //  std::unordered_map<std::string, Eigen::Vector6d> vector6dVariableCache_;
};

#endif//TUDATPY_ORBIT_H
