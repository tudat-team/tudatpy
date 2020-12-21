#include "state.h"
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/conversions.h>
#include <tudat/interface/spice.h>

namespace toec = tudat::orbital_element_conversions;

namespace tudat {

namespace prototype {

BaseState::BaseState(double gravitationalParameter,
                     std::shared_ptr<ReferenceFrame> referenceFrame)
    : referenceFrame_(referenceFrame),
      gravitationalParameter_(gravitationalParameter) {
}

BaseState::BaseState(std::shared_ptr<SimpleBody> body,
                     std::shared_ptr<ReferenceFrame> referenceFrame)
    : referenceFrame_(referenceFrame),
      body_(body) {
}

double BaseState::getGravitationalParameter() { return gravitationalParameter_; }

std::shared_ptr<SimpleBody> BaseState::getBody() { return body_; }

std::shared_ptr<ReferenceFrame> BaseState::getReferenceFrame() { return referenceFrame_; }

const Eigen::Vector6d BaseState::getStateVector() {
  std::logic_error("getStateVector() not Implemented in BaseState class!");
};

const std::shared_ptr<CartesianState> BaseState::toCartesian() {
  std::logic_error("toCartesian() not Implemented in BaseState class!");
};

const std::shared_ptr<KeplerianState> BaseState::toKeplerian() {
  std::logic_error("toKeplerian() not Implemented in BaseState class!");
};

const std::shared_ptr<EquinoctialState> BaseState::toEquinoctial() {
  std::logic_error("toEquinoctial() not Implemented in BaseState class!");
};

const double BaseState::getMeanMotion() {
  return sqrt(getGravitationalParameter() / abs(pow(toKeplerian()->getSemiMajorAxis(), 3)));
};

const std::shared_ptr<BaseState> BaseState::fromSpice(std::string bodyName,
                                                      double ephemerisTime,
                                                      double gravitationalParameter,
                                                      std::string frameOrigin,
                                                      std::string frameOrientation,
                                                      std::string aberrationCorrections) {
  Eigen::Vector6d stateVector = tudat::spice_interface::getBodyCartesianStateAtEpoch(
      bodyName,
      frameOrigin,
      frameOrientation,
      aberrationCorrections,
      ephemerisTime);
  return std::make_shared<CartesianState>(
      stateVector.segment(0, 3),
      stateVector.segment(3, 3),
      gravitationalParameter,
      std::make_shared<ReferenceFrame>(frameOrigin, frameOrientation));
};

const std::shared_ptr<BaseState> BaseState::fromSpice(std::string bodyName,
                                                      double ephemerisTime,
                                                      double gravitationalParameter,
                                                      std::shared_ptr<ReferenceFrame> referenceFrame,
                                                      std::string aberrationCorrections) {
  Eigen::Vector6d stateVector = tudat::spice_interface::getBodyCartesianStateAtEpoch(
      bodyName,
      referenceFrame->getOrigin(),
      referenceFrame->getOrientation(),
      aberrationCorrections,
      ephemerisTime);
  return std::make_shared<CartesianState>(
      stateVector.segment(0, 3),
      stateVector.segment(3, 3),
      gravitationalParameter,
      referenceFrame);
};

CartesianState::CartesianState(
    Eigen::Vector3d positionVector,
    Eigen::Vector3d velocityVector,
    double gravitationalParameter,
    std::shared_ptr<ReferenceFrame> referenceFrame)
    : BaseState(gravitationalParameter, referenceFrame),
      positionVector_(positionVector),
      velocityVector_(velocityVector){};

CartesianState::CartesianState(
    Eigen::Vector3d positionVector,
    Eigen::Vector3d velocityVector,
    std::shared_ptr<SimpleBody> body,
    std::shared_ptr<ReferenceFrame> referenceFrame)
    : BaseState(body, referenceFrame),
      positionVector_(positionVector),
      velocityVector_(velocityVector){};

const Eigen::Vector3d CartesianState::getPositionVector() { return positionVector_; }

const Eigen::Vector3d CartesianState::getVelocityVector() { return velocityVector_; }

const Eigen::Vector6d CartesianState::getStateVector() {
  Eigen::Vector6d cartesianStateVector;
  cartesianStateVector << positionVector_, velocityVector_;
  return cartesianStateVector;
}

const std::string CartesianState::getString(int precision) {
  // __str__ is used to find the “informal” (readable) string representation of an object
  std::ostringstream out;
  out << "---------------"
      << "\n";
  out << "Cartesian State"
      << "\n";
  out << "---------------"
      << "\n";
  out << "position_x = " << positionVector_(0) << " [m]"
      << "\n";
  out << "position_y = " << positionVector_(1) << " [m]"
      << "\n";
  out << "position_z = " << positionVector_(2) << " [m]"
      << "\n";
  out << "velocity_x = " << velocityVector_(0) << " [m/s]"
      << "\n";
  out << "velocity_y = " << velocityVector_(1) << " [m/s]"
      << "\n";
  out << "velocity_z = " << velocityVector_(2) << " [m/s]"
      << "\n";
  return out.str();
}

const std::shared_ptr<EquinoctialState> CartesianState::toEquinoctial() {
  Eigen::Vector6d equinoctialStateVector = toec::convertCartesianToModifiedEquinoctialElements<double>(
      this->getStateVector(),
      this->getGravitationalParameter());
  return std::make_shared<EquinoctialState>(
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::fElementIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::gElementIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::hElementIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::kElementIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::semiParameterIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::trueLongitudeIndex),
      this->getGravitationalParameter(),
      this->getReferenceFrame());
}

const std::shared_ptr<KeplerianState> CartesianState::toKeplerian() {
  Eigen::Vector6d keplerianStateVector = toec::convertCartesianToKeplerianElements<double>(
      this->getStateVector(),
      this->getGravitationalParameter());
  return std::make_shared<KeplerianState>(
      keplerianStateVector(toec::KeplerianElementIndices::semiMajorAxisIndex),
      keplerianStateVector(toec::KeplerianElementIndices::eccentricityIndex),
      keplerianStateVector(toec::KeplerianElementIndices::inclinationIndex),
      keplerianStateVector(toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex),
      keplerianStateVector(toec::KeplerianElementIndices::argumentOfPeriapsisIndex),
      keplerianStateVector(toec::KeplerianElementIndices::trueAnomalyIndex),
      this->getGravitationalParameter(),
      this->getReferenceFrame());
}

const std::shared_ptr<KeplerianState> KeplerianState::propagate(double time,
                                                                std::shared_ptr<tudat::root_finders::RootFinderCore<double>> rootFinder) {
  Eigen::Vector6d keplerianStateVector = toec::propagateKeplerOrbit(
      this->getStateVector(),
      time,
      this->getGravitationalParameter(),
      rootFinder);
  return std::make_shared<KeplerianState>(
      keplerianStateVector(toec::KeplerianElementIndices::semiMajorAxisIndex),
      keplerianStateVector(toec::KeplerianElementIndices::eccentricityIndex),
      keplerianStateVector(toec::KeplerianElementIndices::inclinationIndex),
      keplerianStateVector(toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex),
      keplerianStateVector(toec::KeplerianElementIndices::argumentOfPeriapsisIndex),
      keplerianStateVector(toec::KeplerianElementIndices::trueAnomalyIndex),
      this->getGravitationalParameter(),
      this->getReferenceFrame());
}

const std::shared_ptr<CartesianState> CartesianState::toCartesian() {
  return shared_from_this();
}

KeplerianState::KeplerianState(
    double semiMajorAxis,
    double eccentricity,
    double inclination,
    double longitudeAscendingNode,
    double argumentOfPeriapsis,
    double trueAnomaly,
    double gravitationalParameter,
    std::shared_ptr<ReferenceFrame> referenceFrame)
    : BaseState(gravitationalParameter, referenceFrame),
      semiMajorAxis_(semiMajorAxis),
      eccentricity_(eccentricity),
      inclination_(inclination),
      longitudeAscendingNode_(longitudeAscendingNode),
      argumentOfPeriapsis_(argumentOfPeriapsis),
      trueAnomaly_(trueAnomaly) {}

KeplerianState::KeplerianState(
    double semiMajorAxis,
    double eccentricity,
    double inclination,
    double longitudeAscendingNode,
    double argumentOfPeriapsis,
    double trueAnomaly,
    std::shared_ptr<SimpleBody> body,
    std::shared_ptr<ReferenceFrame> referenceFrame)
    : BaseState(body, referenceFrame),
      semiMajorAxis_(semiMajorAxis),
      eccentricity_(eccentricity),
      inclination_(inclination),
      longitudeAscendingNode_(longitudeAscendingNode),
      argumentOfPeriapsis_(argumentOfPeriapsis),
      trueAnomaly_(trueAnomaly) {}

const double KeplerianState::getSemiMajorAxis() { return semiMajorAxis_; }
const double KeplerianState::getEccentricity() { return eccentricity_; }
const double KeplerianState::getInclination() { return inclination_; }
const double KeplerianState::getLongitudeAscendingNode() { return longitudeAscendingNode_; }
const double KeplerianState::getArgumentOfPeriapsis() { return argumentOfPeriapsis_; }
const double KeplerianState::getTrueAnomaly() { return trueAnomaly_; }

const Eigen::Vector6d KeplerianState::getStateVector() {
  Eigen::Vector6d keplerianStateVector;
  keplerianStateVector(toec::KeplerianElementIndices::semiMajorAxisIndex) = getSemiMajorAxis();
  keplerianStateVector(toec::KeplerianElementIndices::eccentricityIndex) = getEccentricity();
  keplerianStateVector(toec::KeplerianElementIndices::inclinationIndex) = getInclination();
  keplerianStateVector(toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex) = getLongitudeAscendingNode();
  keplerianStateVector(toec::KeplerianElementIndices::argumentOfPeriapsisIndex) = getArgumentOfPeriapsis();
  keplerianStateVector(toec::KeplerianElementIndices::trueAnomalyIndex) = getTrueAnomaly();
  return keplerianStateVector;
}

const std::string KeplerianState::getString(int precision) {
  // __str__ is used to find the “informal” (readable) string representation of an object
  int scale = pow(10, precision);
  std::ostringstream out;
  out << "---------------"
      << "\n";
  out << "Keplerian State"
      << "\n";
  out << "---------------"
      << "\n";
  out << "semi_major_axis = " << round(semiMajorAxis_ * scale) / scale << " [m]"
      << "\n";
  out << "eccentricity = " << round(eccentricity_ * scale) / scale << " [rad]"
      << "\n";
  out << "inclination = " << round(inclination_ * scale) / scale << " [rad]"
      << "\n";
  out << "longitude_ascending_node = " << round(longitudeAscendingNode_ * scale) / scale << " [rad]"
      << "\n";
  out << "argument_of_periapsis = " << round(argumentOfPeriapsis_ * scale) / scale << " [rad]"
      << "\n";
  out << "true_anomaly = " << round(trueAnomaly_ * scale) / scale << " [rad]"
      << "\n";
  return out.str();
}
//// __repr__ is used to find the “official” string representation of an object
//std::ostringstream out;
//std::string repr = "hi";
//return repr;

const std::shared_ptr<CartesianState> KeplerianState::toCartesian() {
  Eigen::Vector6d cartesianStateVector = toec::convertKeplerianToCartesianElements<double>(
      this->getStateVector(),
      this->getGravitationalParameter());
  return std::make_shared<CartesianState>(
      cartesianStateVector.head(3),
      cartesianStateVector.tail<3>(),
      getBody(),
      this->getReferenceFrame());
}

const std::shared_ptr<EquinoctialState> KeplerianState::toEquinoctial() {
  Eigen::Vector6d equinoctialStateVector = toec::convertKeplerianToModifiedEquinoctialElements<double>(this->getStateVector());
  return std::make_shared<EquinoctialState>(
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::fElementIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::gElementIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::hElementIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::kElementIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::semiParameterIndex),
      equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::trueLongitudeIndex),
      this->getGravitationalParameter(),
      this->getReferenceFrame());
}
const std::shared_ptr<KeplerianState> KeplerianState::toKeplerian() {
  return shared_from_this();
}

EquinoctialState::EquinoctialState(double fElement,
                                   double gElement,
                                   double hElement,
                                   double kElement,
                                   double semiParameter,
                                   double trueLongitude,
                                   double gravitationalParameter,
                                   std::shared_ptr<ReferenceFrame> referenceFrame)
    : BaseState(gravitationalParameter, referenceFrame),
      fElement_(fElement),
      gElement_(gElement),
      hElement_(hElement),
      kElement_(kElement),
      semiParameter_(semiParameter),
      trueLongitude_(trueLongitude) {
}

EquinoctialState::EquinoctialState(double fElement,
                                   double gElement,
                                   double hElement,
                                   double kElement,
                                   double semiParameter,
                                   double trueLongitude,
                                   std::shared_ptr<SimpleBody> body,
                                   std::shared_ptr<ReferenceFrame> referenceFrame)
    : BaseState(body, referenceFrame),
      fElement_(fElement),
      gElement_(gElement),
      hElement_(hElement),
      kElement_(kElement),
      semiParameter_(semiParameter),
      trueLongitude_(trueLongitude) {
}

const double EquinoctialState::getfElement() { return fElement_; }
const double EquinoctialState::getgElement() { return gElement_; }
const double EquinoctialState::gethElement() { return hElement_; }
const double EquinoctialState::getkElement() { return kElement_; }
const double EquinoctialState::getSemiParameter() { return semiParameter_; }
const double EquinoctialState::getTrueLongitude() { return trueLongitude_; }

const Eigen::Vector6d EquinoctialState::getStateVector() {
  Eigen::Vector6d equinoctialStateVector;
  equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::fElementIndex) = fElement_;
  equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::gElementIndex) = gElement_;
  equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::hElementIndex) = hElement_;
  equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::kElementIndex) = kElement_;
  equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::semiParameterIndex) = semiParameter_;
  equinoctialStateVector(toec::ModifiedEquinoctialElementVectorIndices::trueLongitudeIndex) = trueLongitude_;
  return equinoctialStateVector;
}

const std::shared_ptr<CartesianState> EquinoctialState::toCartesian() {
  Eigen::Vector6d cartesianStateVector = toec::convertModifiedEquinoctialToCartesianElements<double>(
      this->getStateVector(),
      this->getGravitationalParameter(),
      false);
  return std::make_shared<CartesianState>(
      cartesianStateVector.head(3),
      cartesianStateVector.tail<3>(),
      this->getGravitationalParameter(),
      this->getReferenceFrame());
}

const std::shared_ptr<KeplerianState> EquinoctialState::toKeplerian() {
  Eigen::Vector6d keplerianStateVector = toec::convertModifiedEquinoctialToKeplerianElements<double>(
      this->getStateVector(),
      false);
  return std::make_shared<KeplerianState>(
      keplerianStateVector(toec::KeplerianElementIndices::semiMajorAxisIndex),
      keplerianStateVector(toec::KeplerianElementIndices::eccentricityIndex),
      keplerianStateVector(toec::KeplerianElementIndices::inclinationIndex),
      keplerianStateVector(toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex),
      keplerianStateVector(toec::KeplerianElementIndices::argumentOfPeriapsisIndex),
      keplerianStateVector(toec::KeplerianElementIndices::trueAnomalyIndex),
      this->getGravitationalParameter(),
      this->getReferenceFrame());
}

const std::shared_ptr<EquinoctialState> EquinoctialState::toEquinoctial() {
  return shared_from_this();
}

const std::string EquinoctialState::getString(int precision) {
  // __str__ is used to find the “informal” (readable) string representation of an object
  std::ostringstream out;
  out << "-----------------"
      << "\n";
  out << "Equinoctial State"
      << "\n";
  out << "-----------------"
      << "\n";
  out << "f_element = " << getfElement() << "\n";
  out << "g_element = " << getgElement() << "\n";
  out << "h_element = " << gethElement() << "\n";
  out << "k_element = " << getkElement() << "\n";
  out << "semi_parameter = " << getSemiParameter() << "\n";
  out << "true_longitude = " << getTrueLongitude() << "\n";
  return out.str();
}

}// namespace prototype

}// namespace tudat