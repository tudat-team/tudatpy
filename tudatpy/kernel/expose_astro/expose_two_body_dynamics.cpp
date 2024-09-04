/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "docstrings.h"

#include "expose_two_body_dynamics.h"

#include <tudat/astro/mission_segments.h>
#include <tudat/astro/basic_astro.h>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace tms = tudat::mission_segments;
namespace trf = tudat::root_finders;
namespace toec = tudat::orbital_element_conversions;

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

    /* Trampoline (need one for each virtual function) */
    //  Eigen::Vector3d getInertialVelocityAtArrival() {
    //    PYBIND11_OVERLOAD(
    //        Eigen::Vector3d,              /* Return type */
    //        tms::LambertTargeter,         /* Parent class */
    //        getInertialVelocityAtArrival, /* Name of function in C++ (must match Python name) */
    //                                      /* Argument(s) */
    //    );
    //  }

    //  /* Trampoline (need one for each virtual function) */
    //  std::pair<Eigen::Vector3d, Eigen::Vector3d> getInertialVelocityVectors() {
    //    PYBIND11_OVERLOAD(
    //        std::pair<Eigen::Vector3d, Eigen::Vector3d>, /* Return type */
    //        tms::LambertTargeter,                        /* Parent class */
    //        getInertialVelocityVectors,                  /* Name of function in C++ (must match Python name) */
    //                                                     /* Argument(s) */
    //    );
    //  }
};

}// namespace trampoline

namespace astro {
namespace two_body_dynamics {
void expose_two_body_dynamics(py::module &m) {


    m.def("compute_escape_or_capture_delta_v",
          &tms::computeEscapeOrCaptureDeltaV,
          py::arg("gravitational_param"),
          py::arg("semi_major_axis"),
          py::arg("eccentricity"),
          py::arg("excess_velocity"));

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
                 py::arg("gravitational_parameter"))
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

    //////////////////////////////////////////////////////////////////////
    //  keplerPropagator.h (complete)
    /////////////////////////////////////////////////////////////////

    m.def("propagate_kepler_orbit",
          &toec::propagateKeplerOrbit< double >,
          py::arg("initial_kepler_elements"),
          py::arg("propagation_time"),
          py::arg("gravitational_parameter"),
          py::arg("root_finder") = trf::RootFinderPointer( ),
          get_docstring("propagate_kepler_orbit").c_str( ) );


}
} // namespace two_body_dynamics
} // namespace astro
}// namespace tudatpy
