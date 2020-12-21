/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_PLANETARYNODE_H
#define TUDATPY_PLANETARYNODE_H

namespace tudat {

namespace prototype{

class PlanetaryNode {
 public:
  PlanetaryNode(
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

  Eigen::Vector3d getEntranceVelocity();


 private:
  // constructor values
  std::shared_ptr<CartesianState> centralBodyCartesianState_;
  Eigen::Vector3d incomingVelocity;
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


}
}
#endif//TUDATPY_PLANETARYNODE_H
