#ifndef TUDATPY_STATE_H
#define TUDATPY_STATE_H

#include "../../../expose_bodies/prototype/simpleBody.h"
#include "../../prototype/frames.h"
// NOTE(Geoffrey): This will become <tudat/astro/frames.h>
#include <memory>
#include <tudat/astro/mission_segments.h>
#include <tudat/basics/basicTypedefs.h>

using namespace tudat::bodies;

namespace tudat {

namespace prototype {

class KeplerianState;
class EquinoctialState;
class CartesianState;
class SphericalState;

class BaseState {
 public:
  BaseState(double gravitationalParameter = TUDAT_NAN,
            std::shared_ptr<ReferenceFrame> referenceFrame = std::make_shared<ReferenceFrame>());

  BaseState(std::shared_ptr<SimpleBody> body = std::make_shared<SimpleBody>(),
            std::shared_ptr<ReferenceFrame> referenceFrame = std::make_shared<ReferenceFrame>());

  double getGravitationalParameter();
  std::shared_ptr<SimpleBody> getBody();
  std::shared_ptr<ReferenceFrame> getReferenceFrame();
  virtual const Eigen::Vector6d getStateVector();
  virtual const std::shared_ptr<CartesianState> toCartesian();
  virtual const std::shared_ptr<KeplerianState> toKeplerian();
  virtual const std::shared_ptr<EquinoctialState> toEquinoctial();
  const double getMeanMotion();
  static const std::shared_ptr<BaseState> fromSpice(std::string bodyName,
                                                    double ephemerisTime,
                                                    double gravitationalParameter,
                                                    std::string frameOrigin = "SSB",
                                                    std::string frameOrientation = "J2000",
                                                    std::string aberrationCorrections = "none");

  static const std::shared_ptr<BaseState> fromSpice(std::string bodyName,
                                                    double ephemerisTime,
                                                    double gravitationalParameter,
                                                    std::shared_ptr<ReferenceFrame> referenceFrame,
                                                    std::string aberrationCorrections = "none");

 private:
  double gravitationalParameter_;
  std::shared_ptr<ReferenceFrame> referenceFrame_;
  std::shared_ptr<SimpleBody> body_;
};

class CartesianState : public BaseState, public std::enable_shared_from_this<CartesianState> {
 public:
  CartesianState(
      Eigen::Vector3d positionVector,
      Eigen::Vector3d velocityVector,
      double gravitationalParameter = TUDAT_NAN,
      std::shared_ptr<ReferenceFrame> referenceFrame = std::make_shared<ReferenceFrame>());

  CartesianState(
      Eigen::Vector3d positionVector,
      Eigen::Vector3d velocityVector,
      std::shared_ptr<SimpleBody> body = std::make_shared<SimpleBody>(),
      std::shared_ptr<ReferenceFrame> referenceFrame = std::make_shared<ReferenceFrame>());

  const Eigen::Vector3d getPositionVector();
  const Eigen::Vector3d getVelocityVector();
  const Eigen::Vector6d getStateVector();
  const std::string getString(int precision);
  const std::shared_ptr<KeplerianState> toKeplerian() override;
  const std::shared_ptr<EquinoctialState> toEquinoctial() override;
  const std::shared_ptr<CartesianState> toCartesian() override;

 private:
  Eigen::Vector3d positionVector_;
  Eigen::Vector3d velocityVector_;
};

class EquinoctialState : public BaseState, public std::enable_shared_from_this<EquinoctialState> {

 public:
  EquinoctialState(double fElement,
                   double gElement,
                   double hElement,
                   double kElement,
                   double semiParameter,
                   double trueLongitude,
                   double gravitationalParameter = TUDAT_NAN,
                   std::shared_ptr<ReferenceFrame> = std::make_shared<ReferenceFrame>());

  EquinoctialState(double fElement,
                   double gElement,
                   double hElement,
                   double kElement,
                   double semiParameter,
                   double trueLongitude,
                   std::shared_ptr<SimpleBody> body = std::make_shared<SimpleBody>(),
                   std::shared_ptr<ReferenceFrame> = std::make_shared<ReferenceFrame>());

  const double getfElement();
  const double getgElement();
  const double gethElement();
  const double getkElement();
  const double getSemiParameter();
  const double getTrueLongitude();
  const Eigen::Vector6d getStateVector();
  const std::string getString(int precision = 3);
  const std::shared_ptr<KeplerianState> toKeplerian() override;
  const std::shared_ptr<EquinoctialState> toEquinoctial() override;
  const std::shared_ptr<CartesianState> toCartesian() override;

 private:
  const double fElement_;
  const double gElement_;
  const double hElement_;
  const double kElement_;
  const double semiParameter_;
  const double trueLongitude_;
};

class KeplerianState : public BaseState, public std::enable_shared_from_this<KeplerianState> {
 public:
  KeplerianState(
      double semiMajorAxis,
      double eccentricity,
      double inclination,
      double longitudeAscendingNode,
      double argumentOfPeriapsis,
      double trueAnomaly,
      double gravitationalParameter = TUDAT_NAN,
      std::shared_ptr<ReferenceFrame> = std::make_shared<ReferenceFrame>());

  KeplerianState(
      double semiMajorAxis,
      double eccentricity,
      double inclination,
      double longitudeAscendingNode,
      double argumentOfPeriapsis,
      double trueAnomaly,
      std::shared_ptr<SimpleBody> body = std::make_shared<SimpleBody>(),
      std::shared_ptr<ReferenceFrame> = std::make_shared<ReferenceFrame>());

  const double getSemiMajorAxis();
  const double getEccentricity();
  const double getInclination();
  const double getLongitudeAscendingNode();
  const double getArgumentOfPeriapsis();
  const double getTrueAnomaly();
  const Eigen::Vector6d getStateVector();
  const std::string getString(int precision = 3);
  const std::shared_ptr<KeplerianState> toKeplerian() override;
  const std::shared_ptr<EquinoctialState> toEquinoctial() override;
  const std::shared_ptr<CartesianState> toCartesian() override;
  const std::shared_ptr<KeplerianState> propagate(
      double time,
      std::shared_ptr<tudat::root_finders::RootFinderCore<double>> rootFinder = std::shared_ptr<tudat::root_finders::RootFinderCore<double>>());

 private:
  const double semiMajorAxis_;
  const double eccentricity_;
  const double inclination_;
  const double longitudeAscendingNode_;
  const double argumentOfPeriapsis_;
  const double trueAnomaly_;
};

}// namespace prototype

}// namespace tudat

#endif// TUDATPY_STATE_H
