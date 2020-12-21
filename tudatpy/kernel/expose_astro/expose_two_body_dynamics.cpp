/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_two_body_dynamics.h"

#include "expose_two_body/expose_iod.h"
#include "expose_two_body/expose_state.h"

#include "expose_two_body/prototype/gravityAssist.h"
//#include "expose_two_body/prototype/orbit.h"
#include "expose_two_body/prototype/planetaryDeparture.h"
#include "expose_two_body/prototype/planetaryRendezvous.h"
#include "expose_two_body/prototype/state.h"
#include "prototype/frames.h"

// TODO: Serious thought needs to be placed into whether a simple body can be consistently converted into a "simulation body".
#include "../expose_bodies/prototype/simpleBody.h"

using namespace tudat::bodies;

#include <tudat/astro/conversions.h>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <docstrings/astro/two_body.hpp>

#include <tudat/astro/mission_segments.h>

// needed for std::optional

// INCLUDES STOLEN FROM PGA
#include "tudat/astro/mission_segments/gravityAssist.h"
#include "tudat/math/basic/functionProxy.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/root_finders/bisection.h"
#include <Eigen/Dense>
#include <boost/bind.hpp>
#include <cmath>

namespace py = pybind11;
namespace tms = tudat::mission_segments;
namespace trf = tudat::root_finders;
namespace toec = tudat::orbital_element_conversions;

using namespace tudat::prototype;

namespace tudatpy {

namespace trampoline {

class PyLambertTargeter : public tms::LambertTargeter {
 public:
  /* Inherit the constructors */
  using tms::LambertTargeter::LambertTargeter;

  void execute() override {
    PYBIND11_OVERLOAD_PURE(
        void,                 /* Return type */
        tms::LambertTargeter, /* Parent class */
        execute,              /* Name of function in C++ (must match Python name) */
                              /* Argument(s) */
    );
  }

  /* Trampoline (need one for each virtual function) */
  Eigen::Vector3d get_departure_velocity() {
    PYBIND11_OVERLOAD(
        Eigen::Vector3d,                /* Return type */
        tms::LambertTargeter,           /* Parent class */
        getInertialVelocityAtDeparture, /* Name of function in C++ (must match Python name) */
                                        /* Argument(s) */
    );
  }
};
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

void expose_two_body_dynamics(py::module &m) {
  /*
    * mission_segments/
    * ├── escapeAndCapture.cpp
    * ├── gravityAssist.cpp
    * ├── improvedInversePolynomialWall.cpp
    * ├── lambertRoutines.cpp
    * ├── lambertTargeterGooding.cpp
    * ├── lambertTargeterIzzo.cpp
    * ├── multiRevolutionLambertTargeterIzzo.cpp
    * ├── oscillatingFunctionNovak.cpp
    * └── zeroRevolutionLambertTargeterIzzo.cpp
    *
    * mission_segments
    * ├── escapeAndCapture.h
    * ├── gravityAssist.h
    * ├── improvedInversePolynomialWall.h
    * ├── lambertRoutines.h
    * ├── lambertTargeterGooding.h
    * ├── lambertTargeter.h
    * ├── lambertTargeterIzzo.h
    * ├── multiRevolutionLambertTargeterIzzo.h
    * ├── oscillatingFunctionNovak.h
    * └── zeroRevolutionLambertTargeterIzzo.h
    */

//  py::class_<Orbit, std::shared_ptr<Orbit>>(m, "Orbit")
//      .def(py::init<
//               std::shared_ptr<SimpleBody>,
//               std::shared_ptr<BaseState>,
//               double,
//               std::shared_ptr<ReferenceFrame>>(),
//           py::arg("central_body"),
//           py::arg("state"),
//           py::arg("epoch"),
//           py::arg("reference_frame"));

  //////////////////////////////////////////////////////////////////////
  // prototype [submodule] / state
  //////////////////////////////////////////////////////////////////////
//  auto state = m.def_submodule("state");
//  expose_state(state);

  //////////////////////////////////////////////////////////////////////
  // prototype [submodule: iod] / planetaryDeparture
  //////////////////////////////////////////////////////////////////////
//  auto iod = m.def_submodule("iod");
//  expose_iod(iod);

  py::class_<ArrivalParameters>(m, "ArrivalParameters")
      .def_property_readonly("radius_of_periapsis", [](ArrivalParameters &a) { return a.radiusOfPeriapsis; })
      .def_property_readonly("initial_velocity_at_periapsis", [](ArrivalParameters &a) { return a.initialVelocityAtPeriapsis; })
      .def_property_readonly("incoming_semi_major_axis", [](ArrivalParameters &a) { return a.incomingSemiMajorAxis; })
      .def_property_readonly("incoming_eccentricity", [](ArrivalParameters &a) { return a.incomingEccentricity; })
      .def_property_readonly("final_velocity_at_periapsis", [](ArrivalParameters &a) { return a.finalVelocityAtPeriapsis; })
      .def_property_readonly("total_delta_v", [](ArrivalParameters &a) { return a.totalDeltaV; })
      .def_property_readonly("incoming_hyperbolic_excess_velocity", [](ArrivalParameters &a) { return a.incomingHyperbolicExcessVelocity; });

  py::class_<DepartureParameters>(m, "DepartureParameters")
      .def("__str__", [](const DepartureParameters &a) {
        // __str__ is used to find the “informal” (readable) string representation of an object
        std::ostringstream out;
        out << "Departure(rp=" << round(a.radiusOfPeriapsis * 10) / 10 << ", sma=";
        return out.str();
      })
      .def("__repr__", [](const DepartureParameters &a) {
        // __repr__ is used to find the “official” string representation of an object
        std::ostringstream out;
        std::string repr = "hi";
        return repr;
      })
      .def_property_readonly("radius_of_periapsis", [](DepartureParameters &a) { return a.radiusOfPeriapsis; })
      .def_property_readonly("initial_velocity_at_periapsis", [](DepartureParameters &a) { return a.initialVelocityAtPeriapsis; })
      .def_property_readonly("outgoing_semi_major_axis", [](DepartureParameters &a) { return a.outgoingSemiMajorAxis; })
      .def_property_readonly("outgoing_eccentricity", [](DepartureParameters &a) { return a.outgoingEccentricity; })
      .def_property_readonly("final_velocity_at_periapsis", [](DepartureParameters &a) { return a.finalVelocityAtPeriapsis; })
      .def_property_readonly("total_delta_v", [](DepartureParameters &a) { return a.totalDeltaV; })
      .def_property_readonly("outgoing_hyperbolic_excess_velocity", [](DepartureParameters &a) { return a.outgoingHyperbolicExcessVelocity; });

  py::class_<GravityAssistParameters>(m, "GravityAssistParameters")
      .def_property_readonly("radius_of_periapsis", [](GravityAssistParameters &a) { return a.radiusOfPeriapsis; })
      .def_property_readonly("incoming_semi_major_axis", [](GravityAssistParameters &a) { return a.incomingSemiMajorAxis; })
      .def_property_readonly("incoming_eccentricity", [](GravityAssistParameters &a) { return a.incomingEccentricity; })
      .def_property_readonly("incoming_velocity_at_periapsis", [](GravityAssistParameters &a) { return a.incomingVelocityAtPeriapsis; })
      .def_property_readonly("outgoing_semi_major_axis", [](GravityAssistParameters &a) { return a.outgoingSemiMajorAxis; })
      .def_property_readonly("outgoing_eccentricity", [](GravityAssistParameters &a) { return a.outgoingEccentricity; })
      .def_property_readonly("outgoing_velocity_at_periapsis", [](GravityAssistParameters &a) { return a.outgoingVelocityAtPeriapsis; })
      .def_property_readonly("total_delta_v", [](GravityAssistParameters &a) { return a.totalDeltaV; })
      .def_property_readonly("incoming_hyperbolic_excess_velocity", [](GravityAssistParameters &a) { return a.incomingHyperbolicExcessVelocity; })
      .def_property_readonly("outgoing_hyperbolic_excess_velocity", [](GravityAssistParameters &a) { return a.outgoingHyperbolicExcessVelocity; });

  m.def("calculate_departure_parameters", &calculateDepartureParameters,
        py::arg("central_body_gravitational_parameter"),
        py::arg("central_body_velocity_on_exit_soi"),
        py::arg("outgoing_velocity"),
        py::arg("departure_periapsis_distance"),
        R"mydelimiter(
Planetary departure is concerned when a spacecraft is situated in a near-circular parking orbit around a body and subsequently provides an impulsive manoeuvre to establish a hyperbolic trajectory ($e>1$). The Keplerian parameters are desired in order to analyse the geometry of the complete interplanetary trajectory. The desired parameters for a \textit{basic analysis} are the eccentricity ($e$) and the velocity required at the periapsis ($V_p$).

The orbital parameters are determined using the equations derived by D. Curtis \cite{curtis_2014}. These scalar parameters allow for complete description of the perifocal coordinate system (PQW), otherwise known as the \textit{natural frame} of an orbit \cite{curtis_2014}.

.. math::
   :nowrap:

	\begin{equation}
	\label{eq:pd_e}
	e=1+\frac{r_p\cdot{V_\infty^2}}{\mu}
	\end{equation}

	\begin{equation}
	\label{eq:v_c}
	v_c = \sqrt{\frac{\mu}{r_p}}
	\end{equation}

Given the $v_\infty$ required for the interplanetary leg following departure, \autoref{eq:pd_e} is used to calculated the eccentricity which gives full geometrical description of the orbital conic shape with $r_p$ having been set by the parking orbit altitude.

.. math::
   :nowrap:

	\begin{equation}
	\label{eq:pd_dv}
	\Delta{V} = v_p-v_c
	\end{equation}

	\begin{equation}
	\label{eq:pd_beta}
	\beta = \arccos{\frac{1}{e}}
	\end{equation}

	\begin{equation}
	\label{eq:pd_delta}
	\Delta = a\sqrt{e^2 - 1}
	\end{equation}

	\begin{equation}
    \label{eq:pd_dv}
    \Delta{V} = V_{p,hyp}-V_c
    \end{equation}

Where $V_{p,hyp}$ is defined in \autoref{eq:visviva_hyp}. $V_c$, the orbital velocity of a circular orbit, the $\Delta{V}$ required at periapsis to provide $v_\infty$ at infinity can be calculated as shown in \autoref{eq:pd_dv}.
)mydelimiter");

  m.def("calculate_arrival_parameters", &calculateArrivalParameters,
        py::arg("central_body_gravitational_parameter"),
        py::arg("central_body_velocity_on_entry_soi"),
        py::arg("incoming_velocity"),
        py::arg("arrival_periapsis_distance"));

  m.def("calculate_gravity_assist_parameters", &calculateGravityAssistParameters,
        py::arg("central_body_gravitational_parameter"),
        py::arg("central_body_velocity_on_entry_soi"),
        py::arg("central_body_velocity_on_exit_soi"),
        py::arg("incoming_velocity"),
        py::arg("outgoing_velocity"),
        py::arg("smallest_periapsis_distance"),
        py::arg("use_eccentricity_instead_of_pericenter") = true,
        py::arg("speed_tolerance") = 1e-6,
        py::arg("root_finder") = std::make_shared<trf::NewtonRaphson>(1.0e-12, 1000));

  //////////////////////////////////////////////////////////////////////
  // escapeAndCapture.cpp
  //////////////////////////////////////////////////////////////////////
  m.def("compute_escape_or_capture_delta_v",
        &tms::computeEscapeOrCaptureDeltaV,
        py::arg("gravitational_param"),
        py::arg("semi_major_axis"),
        py::arg("eccentricity"),
        py::arg("excess_velocity"));

  //////////////////////////////////////////////////////////////////////
  //  gravityAssist.cpp
  //////////////////////////////////////////////////////////////////////
  m.def("gravity_assist",// overload 1 (returns double)
        py::overload_cast<const double,
                          const Eigen::Vector3d &,
                          const Eigen::Vector3d &,
                          const Eigen::Vector3d &,
                          const double,
                          const bool,
                          const double,
                          trf::RootFinderPointer>(&tms::gravityAssist),
        py::arg("gravitational_param"),
        py::arg("central_body_velocity"),
        py::arg("incoming_velocity"),
        py::arg("outgoing_velocity"),
        py::arg("smallest_periapsis_distance"),
        py::arg("use_eccentricity_over_pericenter") = true,
        py::arg("speed_tolerance") = 1.0e-6,
        py::arg("root_finder") = std::make_shared<trf::NewtonRaphson>(1.0e-12, 1000));

  m.def("gravity_assist",// overload 2: unassisted (returns 3 dim vector)
        py::overload_cast<const double,
                          const Eigen::Vector3d &,
                          const Eigen::Vector3d &,
                          const double,
                          const double>(&tms::gravityAssist),
        py::arg("gravitational_param"),
        py::arg("central_body_velocity"),
        py::arg("incoming_velocity"),
        py::arg("rotation_angle"),
        py::arg("pericenter_radius"));

  m.def("gravity_assist",// overload 3: assisted (returns 3 dim vector)
        py::overload_cast<const double,
                          const Eigen::Vector3d &,
                          const Eigen::Vector3d &,
                          const double,
                          const double,
                          const double>(&tms::gravityAssist),
        py::arg("gravitational_param"),
        py::arg("central_body_velocity"),
        py::arg("incoming_velocity"),
        py::arg("rotation_angle"),
        py::arg("pericenter_radius"),
        py::arg("delta_v"));

  py::class_<tms::PericenterFindingFunctions,
             std::shared_ptr<tms::PericenterFindingFunctions>>(m, "PericenterFindingFunctions")
      .def(py::init<const double,
                    const double,
                    const double>(),
           py::arg("absolute_incoming_semi_major_axis"),
           py::arg("absolute_outgoing_semi_major_axis"),
           py::arg("bending_angle"))
      .def("compute_pericenter_radius_fn",
           &tms::PericenterFindingFunctions::computePericenterRadiusFunction)
      .def("compute_derivative_pericenter_radius_fn",
           &tms::PericenterFindingFunctions::computeFirstDerivativePericenterRadiusFunction);

  py::class_<tms::EccentricityFindingFunctions,
             std::shared_ptr<tms::EccentricityFindingFunctions>>(m, "EccentricityFindingFunctions")
      .def(py::init<const double,
                    const double,
                    const double>(),
           py::arg("absolute_incoming_semi_major_axis"),
           py::arg("absolute_outgoing_semi_major_axis"),
           py::arg("bending_angle"))
      .def("compute_incoming_eccentricity_fn",
           &tms::EccentricityFindingFunctions::computeIncomingEccentricityFunction)
      .def("compute_derivative_incoming_eccentricity_fn",
           &tms::EccentricityFindingFunctions::computeFirstDerivativeIncomingEccentricityFunction);

  //////////////////////////////////////////////////////////////////////
  //  improvedInversePolynomialWall.cpp
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  //  lambertRoutines.cpp
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  //  lambertTargeter.cpp (complete)
  //////////////////////////////////////////////////////////////////////
  py::class_<tms::LambertTargeter,
             std::shared_ptr<tms::LambertTargeter>,
             trampoline::PyLambertTargeter>(m, "LambertTargeter")
      // This class required a trampoline class to inherit from due to
      // the virtual ~execute member.
      .def(py::init<const Eigen::Vector3d &,
                    const Eigen::Vector3d &,
                    const double,
                    const double>(),
           py::arg("departure_position"),
           py::arg("arrival_position"),
           py::arg("time_of_flight"),
           py::arg("gravitational_parameter"), lambert_targeter_ctor_docstring().c_str())
      .def("get_departure_velocity", &tms::LambertTargeter::getInertialVelocityAtDeparture)
      .def("get_arrival_velocity", &tms::LambertTargeter::getInertialVelocityAtArrival)
      .def("get_velocity_vectors", &tms::LambertTargeter::getInertialVelocityVectors);

  //////////////////////////////////////////////////////////////////////
  //  lambertTargeterGooding.cpp (complete)
  //////////////////////////////////////////////////////////////////////
  py::class_<tms::LambertTargeterGooding,
             std::shared_ptr<tms::LambertTargeterGooding>,
             tms::LambertTargeter>(m, "LambertTargeterGooding")
      .def(py::init<const Eigen::Vector3d &,
                    const Eigen::Vector3d &,
                    const double,
                    const double,
                    trf::RootFinderPointer>(),
           py::arg("departure_position"),
           py::arg("arrival_position"),
           py::arg("time_of_flight"),
           py::arg("gravitational_parameter"),
           py::arg("root_finder") = trf::RootFinderPointer())
      .def("get_radial_departure_velocity", &tms::LambertTargeterGooding::getRadialVelocityAtDeparture)
      .def("get_radial_arrival_velocity", &tms::LambertTargeterGooding::getRadialVelocityAtArrival)
      .def("get_transverse_departure_velocity", &tms::LambertTargeterGooding::getTransverseVelocityAtDeparture)
      .def("get_transverse_arrival_velocity", &tms::LambertTargeterGooding::getTransverseVelocityAtArrival)
      .def("get_semi_major_axis", &tms::LambertTargeterGooding::getSemiMajorAxis);

  //////////////////////////////////////////////////////////////////////
  //  lambertTargeterIzzo.cpp (complete)
  //////////////////////////////////////////////////////////////////////
  py::class_<tms::LambertTargeterIzzo,
             std::shared_ptr<tms::LambertTargeterIzzo>,
             tms::LambertTargeter>(m, "LambertTargeterIzzo")
      .def(py::init<const Eigen::Vector3d &,
                    const Eigen::Vector3d &,
                    const double,
                    const double,
                    const bool,
                    const double,
                    const int>(),
           py::arg("departure_position"),
           py::arg("arrival_position"),
           py::arg("time_of_flight"),
           py::arg("gravitational_parameter"),
           py::arg("is_retrograde") = false,
           py::arg("tolerance") = 1e-9,
           py::arg("max_iter") = 50)
      .def("get_radial_departure_velocity", &tms::LambertTargeterIzzo::getRadialVelocityAtDeparture)
      .def("get_radial_arrival_velocity", &tms::LambertTargeterIzzo::getRadialVelocityAtArrival)
      .def("get_transverse_departure_velocity", &tms::LambertTargeterIzzo::getTransverseVelocityAtDeparture)
      .def("get_transverse_arrival_velocity", &tms::LambertTargeterIzzo::getTransverseVelocityAtArrival)
      .def("get_semi_major_axis", &tms::LambertTargeterIzzo::getSemiMajorAxis);

  //////////////////////////////////////////////////////////////////////
  //  oscillatingFunctionNovak.cpp
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  //  zeroRevolutionLambertTargeterIzzo.cpp (complete)
  //////////////////////////////////////////////////////////////////////
  py::class_<tms::ZeroRevolutionLambertTargeterIzzo,
             std::shared_ptr<tms::ZeroRevolutionLambertTargeterIzzo>,
             tms::LambertTargeter>(m, "ZeroRevolutionLambertTargeterIzzo")
      .def(py::init<const Eigen::Vector3d &,
                    const Eigen::Vector3d &,
                    const double,
                    const double,
                    const bool,
                    const double,
                    const int>(),
           py::arg("departure_position"),
           py::arg("arrival_position"),
           py::arg("time_of_flight"),
           py::arg("gravitational_parameter"),
           py::arg("is_retrograde") = false,
           py::arg("tolerance") = 1e-9,
           py::arg("max_iter") = 50)
      .def("get_radial_departure_velocity", &tms::ZeroRevolutionLambertTargeterIzzo::getRadialVelocityAtDeparture)
      .def("get_radial_arrival_velocity", &tms::ZeroRevolutionLambertTargeterIzzo::getRadialVelocityAtArrival)
      .def("get_transverse_departure_velocity", &tms::ZeroRevolutionLambertTargeterIzzo::getTransverseVelocityAtDeparture)
      .def("get_transverse_arrival_velocity", &tms::ZeroRevolutionLambertTargeterIzzo::getTransverseVelocityAtArrival)
      .def("get_semi_major_axis", &tms::ZeroRevolutionLambertTargeterIzzo::getSemiMajorAxis);

  //////////////////////////////////////////////////////////////////////
  //  multiRevolutionLambertTargeterIzzo.cpp (complete)
  //////////////////////////////////////////////////////////////////////
  py::class_<tms::MultiRevolutionLambertTargeterIzzo,
             std::shared_ptr<tms::MultiRevolutionLambertTargeterIzzo>,
             tms::ZeroRevolutionLambertTargeterIzzo>(m, "MultiRevolutionLambertTargeterIzzo")
      .def(py::init<const Eigen::Vector3d &,
                    const Eigen::Vector3d &,
                    const double,
                    const double,
                    const int,
                    const bool,
                    const bool,
                    const double,
                    const int>(),
           py::arg("departure_position"),
           py::arg("arrival_position"),
           py::arg("time_of_flight"),
           py::arg("gravitational_parameter"),
           py::arg("n_revolutions") = 0,
           py::arg("is_right_branch") = false,
           py::arg("is_retrograde") = false,
           py::arg("tolerance") = 1e-9,
           py::arg("max_iter") = 50)
      //        .def("NO_MAXIMUM_REVOLUTIONS", &tms::MultiRevolutionLambertTargeterIzzo::NO_MAXIMUM_REVOLUTIONS)
      .def("compute_for_revolutions_and_branch", &tms::MultiRevolutionLambertTargeterIzzo::getRadialVelocityAtArrival)
      .def("get_max_n_revolutions", &tms::MultiRevolutionLambertTargeterIzzo::getRadialVelocityAtArrival);
};

}// namespace tudatpy
