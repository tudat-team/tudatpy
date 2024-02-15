/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_ephemerides.h"

#include <tudat/astro/ephemerides.h>
#include <tudat/simulation/simulation.h>// TODO: EphemerisType should be in <tudat/astro/ephemerides.h>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace te = tudat::ephemerides;
namespace tss = tudat::simulation_setup;

namespace tudatpy {

namespace trampoline {

class PyEphemeris : public te::Ephemeris {
 public:
  /* Inherit the constructors */
  using te::Ephemeris::Ephemeris;

  Eigen::Vector6d getCartesianState(const double seconds_since_epoch) override {
    PYBIND11_OVERLOAD_PURE(
        Eigen::Vector6d,    /* Return type */
        te::Ephemeris,      /* Parent class */
        getCartesianState,  /* Name of function in C++ (must match Python name) */
        seconds_since_epoch /* Argument(s) */
    );
  }

  //  Eigen::Matrix<long double, 6, 1> getCartesianLongState(const double seconds_since_epoch) override {
  //    PYBIND11_OVERLOAD_PURE(
  //        Eigen::Matrix<long double, 6, 1>, /* Return type */
  //        te::Ephemeris,                    /* Parent class */
  //        getCartesianLongState,            /* Name of function in C++ (must match Python name) */
  //        seconds_since_epoch               /* Argument(s) */
  //    );
  //  }

  //  Eigen::Matrix<long double, 6, 1> getCartesianStateFromExtendedTime() override {
  //    PYBIND11_OVERLOAD_PURE(
  //        Eigen::Matrix<long double, 6, 1>,  /* Return type */
  //        te::Ephemeris,                     /* Parent class */
  //        getCartesianStateFromExtendedTime, /* Name of function in C++ (must match Python name) */
  //        const Time &                       /* Argument(s) */
  //    );
  //  }

  //  Eigen::Matrix<long double, 6, 1> getCartesianStateFromExtendedTime() override {
  //    PYBIND11_OVERLOAD_PURE(
  //        Eigen::Matrix<long double, 6, 1>,  /* Return type */
  //        te::Ephemeris,                     /* Parent class */
  //        getCartesianStateFromExtendedTime, /* Name of function in C++ (must match Python name) */
  //        const double                       /* Argument(s) */
  //    );
};

}// namespace trampoline

void expose_ephemerides(py::module &m) {
  /*
   * ephemerides/
   *  ├── approximatePlanetPositionsBase.h
   *  ├── approximatePlanetPositionsCircularCoplanar.h
   *  ├── approximatePlanetPositionsDataContainer.h
   *  ├── approximatePlanetPositions.h
   *  ├── cartesianStateExtractor.h
   *  ├── compositeEphemeris.h
   *  ├── constantEphemeris.h
   *  ├── constantRotationalEphemeris.h
   *  ├── customEphemeris.h
   *  ├── ephemeris.h
   *  ├── frameManager.h
   *  ├── fullPlanetaryRotationModel.h
   *  ├── itrsToGcrsRotationModel.h
   *  ├── keplerEphemeris.h
   *  ├── keplerStateExtractor.h
   *  ├── multiArcEphemeris.h
   *  ├── rotationalEphemeris.h
   *  ├── simpleRotationalEphemeris.h
   *  ├── synchronousRotationalEphemeris.h
   *  ├── tabulatedEphemeris.h
   *  ├── tabulatedRotationalEphemeris.h
   *  └── tleEphemeris.h
   *
   * ephemerides/
   *  ├── approximatePlanetPositionsBase.cpp
   *  ├── approximatePlanetPositionsCircularCoplanar.cpp
   *  ├── approximatePlanetPositions.cpp
   *  ├── cartesianStateExtractor.cpp
   *  ├── CMakeLists.txt
   *  ├── compositeEphemeris.cpp
   *  ├── ephemeris.cpp
   *  ├── frameManager.cpp
   *  ├── fullPlanetaryRotationModel.cpp
   *  ├── keplerEphemeris.cpp
   *  ├── keplerStateExtractor.cpp
   *  ├── rotationalEphemeris.cpp
   *  ├── simpleRotationalEphemeris.cpp
   *  ├── synchronousRotationalEphemeris.cpp
   *  ├── tabulatedEphemeris.cpp
   *  ├── tabulatedRotationalEphemeris.cpp
   *  └── tleEphemeris.cpp
   *
   */
  //////////////////////////////////////////////////////////////////////////////
  //
  //////////////////////////////////////////////////////////////////////////////
  py::class_<te::Ephemeris, std::shared_ptr<te::Ephemeris>>(m, "Ephemeris")
      .def("get_cartesian_state", &te::Ephemeris::getCartesianState, py::arg("seconds_since_epoch") = 0.0)
      .def("get_cartesian_position", &te::Ephemeris::getCartesianPosition, py::arg("seconds_since_epoch") = 0.0)
      .def("get_cartesian_velocity", &te::Ephemeris::getCartesianVelocity, py::arg("seconds_since_epoch") = 0.0);

  py::enum_<tss::EphemerisType>(m.attr("Ephemeris"), "EphemerisType")
      .value("approximate_planet_positions", tss::approximate_planet_positions)
      .value("direct_spice_ephemeris", tss::direct_spice_ephemeris)
      .value("interpolated_spice", tss::interpolated_spice)
      .value("constant_ephemeris", tss::constant_ephemeris)
      .value("kepler_ephemeris", tss::kepler_ephemeris)
      .value("direct_tle_ephemeris", tss::direct_tle_ephemeris)
      .value("interpolated_tle_ephemeris", tss::interpolated_tle_ephemeris)
      .value("custom_ephemeris", tss::custom_ephemeris);

  //////////////////////////////////////////////////////////////////////////////
  // constantEphemeris.h
  //////////////////////////////////////////////////////////////////////////////
  py::class_<te::ConstantEphemeris,
             std::shared_ptr<te::ConstantEphemeris>,
             te::Ephemeris>(
      m, "ConstantEphemeris")
      .def(py::init<
               const std::function<Eigen::Vector6d()>,//<pybind11/functional.h>,<pybind11/eigen.h>
               const std::string &,
               const std::string &>(),
           py::arg("constant_state_function"),
           py::arg("reference_frame_origin") = "SSB",
           py::arg("reference_frame_orientation") = "ECLIPJ2000")
      .def(py::init<
               const Eigen::Vector6d,//<pybind11/eigen.h>
               const std::string &,
               const std::string &>(),
           py::arg("constant_state"),
           py::arg("reference_frame_origin") = "SSB",
           py::arg("reference_frame_orientation") = "ECLIPJ2000")
      .def("update_constant_state", &te::ConstantEphemeris::updateConstantState,
           py::arg("new_state"));

  //////////////////////////////////////////////////////////////////////////////
  // keplerEphemeris.h
  //////////////////////////////////////////////////////////////////////////////

  py::class_<te::KeplerEphemeris,
             std::shared_ptr<te::KeplerEphemeris>,
             te::Ephemeris>(
      m, "KeplerEphemeris");

  //////////////////////////////////////////////////////////////////////////////
  // tleEphemeris.h / tleEphemeris.cpp
  //////////////////////////////////////////////////////////////////////////////
  py::class_<te::Tle, std::shared_ptr<te::Tle>>(m, "Tle")
      .def(py::init<//ctor 1
               const std::string &>(),
           py::arg("lines"))
      .def(py::init<//ctor 2
               const std::string &,
               const std::string &>(),
           py::arg("line_1"),
           py::arg("line_2"))
      .def("get_epoch", &te::Tle::getEpoch)
      .def("get_b_star", &te::Tle::getBStar)
      .def("get_epoch", &te::Tle::getEpoch)
      .def("get_inclination", &te::Tle::getInclination)
      .def("get_right_ascension", &te::Tle::getRightAscension)
      .def("get_eccentricity", &te::Tle::getEccentricity)
      .def("get_arg_of_perigee", &te::Tle::getArgOfPerigee)
      .def("get_mean_anomaly", &te::Tle::getMeanAnomaly)
      .def("get_mean_motion", &te::Tle::getMeanMotion);

  py::class_<te::TleEphemeris,
             std::shared_ptr<te::TleEphemeris>,
             te::Ephemeris>(m, "TleEphemeris")
      .def(py::init<
               const std::string &,
               const std::string &,
               const std::shared_ptr<te::Tle>,
               const bool>(),
           py::arg("frame_origin") = "Earth",
           py::arg("frame_orientation") = "J2000",
           py::arg("tle") = nullptr,
           py::arg("use_sdp") = false);

  //////////////////////////////////////////////////////////////////////////////
  // tabulatedEphemeris.h
  //////////////////////////////////////////////////////////////////////////////
  py::class_<te::TabulatedCartesianEphemeris< double, double >,
             std::shared_ptr<te::TabulatedCartesianEphemeris< double, double > >,
             te::Ephemeris>(m, "TabulatedEphemeris")
             .def("reset_interpolator", &te::TabulatedCartesianEphemeris< double, double >::resetInterpolator,
                  py::arg("interpolator") );

};
}// namespace tudatpy
