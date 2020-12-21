/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_state.h"
#include "prototype/state.h"
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace tudat::prototype;

namespace tudatpy {

namespace trampoline {
class PyBaseState : public BaseState {
 public:
  /* Inherit the constructors */
  using BaseState::BaseState;

  const Eigen::Vector6d getStateVector() override {
    PYBIND11_OVERLOAD_PURE(
        const Eigen::Vector6d, /* Return type */
        BaseState,             /* Parent class */
        getStateVector,        /* Name of function in C++ (must match Python name) */
                               /* Argument(s) */
    );
  }

  const std::shared_ptr<CartesianState> toCartesian() override {
    PYBIND11_OVERLOAD_PURE(
        const std::shared_ptr<CartesianState>, /* Return type */
        BaseState,                             /* Parent class */
        toCartesian,                           /* Name of function in C++ (must match Python name) */
                                               /* Argument(s) */
    );
  }

  const std::shared_ptr<KeplerianState> toKeplerian() override {
    PYBIND11_OVERLOAD_PURE(
        const std::shared_ptr<KeplerianState>, /* Return type */
        BaseState,                             /* Parent class */
        toKeplerian,                           /* Name of function in C++ (must match Python name) */
                                               /* Argument(s) */
    );
  }

  const std::shared_ptr<EquinoctialState> toEquinoctial() override {
    PYBIND11_OVERLOAD_PURE(
        const std::shared_ptr<EquinoctialState>, /* Return type */
        BaseState,                               /* Parent class */
        toEquinoctial,                           /* Name of function in C++ (must match Python name) */
                                                 /* Argument(s) */
    );
  };

  const std::shared_ptr<BaseState> fromSpice(std::string bodyName,
                                             double ephemerisTime,
                                             double gravitationalParameter,
                                             std::string frameOrigin = "SSB",
                                             std::string frameOrientation = "J2000",
                                             std::string aberrationCorrections = "none") {
    PYBIND11_OVERLOAD_PURE(
        const std::shared_ptr<BaseState>,                                                                     /* Return type */
        BaseState,                                                                                            /* Parent class */
        fromSpice,                                                                                            /* Name of function in C++ (must match Python name) */
        bodyName, ephemerisTime, gravitationalParameter, frameOrigin, frameOrientation, aberrationCorrections /* Argument(s) */
    );
  }

  const std::shared_ptr<BaseState> fromSpice(std::string bodyName,
                                             double ephemerisTime,
                                             double gravitationalParameter,
                                             ReferenceFrame referenceFrame,
                                             std::string aberrationCorrections = "none") {
    PYBIND11_OVERLOAD_PURE(
        const std::shared_ptr<BaseState>,                                                      /* Return type */
        BaseState,                                                                             /* Parent class */
        fromSpice,                                                                             /* Name of function in C++ (must match Python name) */
        bodyName, ephemerisTime, gravitationalParameter, referenceFrame, aberrationCorrections /* Argument(s) */
    );
  }
};
}// namespace trampoline

void expose_state(py::module &m) {

  py::class_<BaseState, std::shared_ptr<BaseState>, trampoline::PyBaseState>(m, "_BaseState")
      .def_property_readonly("state_vector", &BaseState::getStateVector)
      .def_property_readonly("mean_motion", &BaseState::getMeanMotion)
      .def_property_readonly("gravitational_parameter", &BaseState::getGravitationalParameter)
      .def_static("from_spice", py::overload_cast<std::string, double, double, std::string, std::string, std::string>(&BaseState::fromSpice),
                  py::arg("body_name"),
                  py::arg("ephemeris_time"),
                  py::arg("gravitational_parameter"),
                  py::arg("frame_origin") = "SSB",
                  py::arg("frame_orientation") = "J2000",
                  py::arg("aberration_corrections") = "none")
      .def_static("from_spice", py::overload_cast<std::string, double, double, std::shared_ptr<ReferenceFrame>, std::string>(&BaseState::fromSpice),
                  py::arg("body_name"),
                  py::arg("ephemeris_time"),
                  py::arg("gravitational_parameter"),
                  py::arg("reference_frame"),
                  py::arg("aberration_corrections") = "none");

  py::class_<KeplerianState, std::shared_ptr<KeplerianState>, BaseState>(m, "KeplerianState")
      .def(py::init<double, double, double, double, double, double, double, std::shared_ptr<ReferenceFrame>>(),
           py::arg("semi_major_axis"),
           py::arg("eccentricity"),
           py::arg("inclination"),
           py::arg("longitude_ascending_node"),
           py::arg("argument_of_periapsis"),
           py::arg("true_anomaly"),
           py::arg("gravitational_parameter") = TUDAT_NAN,
           py::arg("reference_frame") = std::make_shared<ReferenceFrame>())
      .def(py::init<double, double, double, double, double, double, double, std::shared_ptr<ReferenceFrame>>(),
           py::arg("semi_major_axis"),
           py::arg("eccentricity"),
           py::arg("inclination"),
           py::arg("longitude_ascending_node"),
           py::arg("argument_of_periapsis"),
           py::arg("true_anomaly"),
           py::arg("central_body") = std::make_shared<SimpleBody>(),
           py::arg("reference_frame") = std::make_shared<ReferenceFrame>())
      .def_property_readonly("semi_major_axis", &KeplerianState::getSemiMajorAxis)
      .def_property_readonly("eccentricity", &KeplerianState::getEccentricity)
      .def_property_readonly("inclination", &KeplerianState::getInclination)
      .def_property_readonly("longitude_ascending_node", &KeplerianState::getLongitudeAscendingNode)
      .def_property_readonly("argument_of_periapsis", &KeplerianState::getArgumentOfPeriapsis)
      .def_property_readonly("true_anomaly", &KeplerianState::getTrueAnomaly)
      .def("__str__", &KeplerianState::getString, py::arg("precision") = 1)
      .def("to_cartesian", &KeplerianState::toCartesian)
      .def("to_keplerian", &KeplerianState::toKeplerian)
      .def("to_equinoctial", &KeplerianState::toEquinoctial)
      .def("propagate", &KeplerianState::propagate,
           py::arg("time"),
           py::arg("root_finder") = std::shared_ptr<tudat::root_finders::RootFinderCore<double>>());

  py::class_<EquinoctialState, std::shared_ptr<EquinoctialState>, BaseState>(m, "EquinoctialState")
      .def(py::init<double, double, double, double, double, double, double, std::shared_ptr<ReferenceFrame>>(),
           py::arg("f_element"),
           py::arg("g_element"),
           py::arg("h_element"),
           py::arg("k_element"),
           py::arg("semi_parameter"),
           py::arg("true_longitude"),
           py::arg("gravitational_parameter") = TUDAT_NAN,
           py::arg("reference_frame") = std::make_shared<ReferenceFrame>())
      .def(py::init<double, double, double, double, double, double, double, std::shared_ptr<ReferenceFrame>>(),
           py::arg("f_element"),
           py::arg("g_element"),
           py::arg("h_element"),
           py::arg("k_element"),
           py::arg("semi_parameter"),
           py::arg("true_longitude"),
           py::arg("central_body") = std::make_shared<SimpleBody>(),
           py::arg("reference_frame") = std::make_shared<ReferenceFrame>())
      .def_property_readonly("f_element", &EquinoctialState::getGravitationalParameter)
      .def_property_readonly("g_element", &EquinoctialState::getgElement)
      .def_property_readonly("h_element", &EquinoctialState::gethElement)
      .def_property_readonly("k_element", &EquinoctialState::getkElement)
      .def_property_readonly("semi_parameter", &EquinoctialState::getSemiParameter)
      .def_property_readonly("true_longitude", &EquinoctialState::getTrueLongitude)
      .def("__str__", &EquinoctialState::getString, py::arg("precision") = 1)
      .def("to_cartesian", &EquinoctialState::toCartesian)
      .def("to_equinoctial", &EquinoctialState::toEquinoctial)
      .def("to_keplerian", &EquinoctialState::toKeplerian);

  py::class_<CartesianState, std::shared_ptr<CartesianState>, BaseState>(m, "CartesianState")
      .def(py::init<Eigen::Vector3d, Eigen::Vector3d, double, std::shared_ptr<ReferenceFrame>>(),
           py::arg("position_vector"),
           py::arg("velocity_vector"),
           py::arg("gravitational_parameter") = TUDAT_NAN,
           py::arg("reference_frame") = std::make_shared<ReferenceFrame>())
      .def(py::init<Eigen::Vector3d, Eigen::Vector3d, double, std::shared_ptr<ReferenceFrame>>(),
           py::arg("position_vector"),
           py::arg("velocity_vector"),
           py::arg("central_body") = std::make_shared<SimpleBody>(),
           py::arg("reference_frame") = std::make_shared<ReferenceFrame>())
      .def_property_readonly("velocity_vector", &CartesianState::getVelocityVector)
      .def_property_readonly("position_vector", &CartesianState::getPositionVector)
      .def("__str__", &CartesianState::getString, py::arg("precision") = 1)
      .def("to_keplerian", &CartesianState::toKeplerian)
      .def("to_cartesian", &CartesianState::toCartesian)
      .def("to_equinoctial", &CartesianState::toEquinoctial);
}
}// namespace tudatpy