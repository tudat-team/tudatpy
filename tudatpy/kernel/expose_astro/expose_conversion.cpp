/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_conversion.h"

#include <tudat/astro/conversions.h>
#include <tudat/astro/ephemerides/rotationalEphemeris.h>

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace toec = tudat::orbital_element_conversions;
namespace te = tudat::ephemerides;

namespace tudatpy {

void expose_conversion(py::module &m) {

    m.def("transform_to_inertial_orientation",
          &te::transformStateToInertialOrientation<double, double>,
          py::arg("state_in_body_fixed_frame"),
          py::arg("current_time"),
          py::arg("rotational_ephemeris"));

    py::enum_<toec::KeplerianElementIndices>(m, "KeplerianElementIndices")
            .value("semi_major_axis_index", toec::KeplerianElementIndices::semiMajorAxisIndex)
            .value("eccentricity_index", toec::KeplerianElementIndices::eccentricityIndex)
            .value("inclination_index", toec::KeplerianElementIndices::inclinationIndex)
            .value("argument_of_periapsis_index", toec::KeplerianElementIndices::argumentOfPeriapsisIndex)
            .value("longitude_of_ascending_node_index", toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex)
            .value("true_anomaly_index", toec::KeplerianElementIndices::trueAnomalyIndex)
            .value("semi_latus_rectum_index", toec::KeplerianElementIndices::semiLatusRectumIndex)
            .export_values();

    py::enum_<toec::SphericalOrbitalStateElementIndices>(m, "SphericalOrbitalStateElementIndices")
            .value("radius_index", toec::SphericalOrbitalStateElementIndices::radiusIndex)
            .value("latitude_index", toec::SphericalOrbitalStateElementIndices::latitudeIndex)
            .value("longitude_index", toec::SphericalOrbitalStateElementIndices::longitudeIndex)
            .value("speed_index", toec::SphericalOrbitalStateElementIndices::speedIndex)
            .value("flight_path_index", toec::SphericalOrbitalStateElementIndices::flightPathIndex)
            .value("heading_angle_index", toec::SphericalOrbitalStateElementIndices::headingAngleIndex)
            .export_values();

    m.def("cartesian_to_keplerian",
          &toec::convertCartesianToKeplerianElements< double >,
          py::arg("cartesian_elements"),
          py::arg("gravitational_parameter"));

    m.def("keplerian_to_cartesian",
          py::overload_cast< const Eigen::Vector6d&, double >(
              &toec::convertKeplerianToCartesianElements< double > ),
          py::arg("keplerian_elements"),
          py::arg("gravitational_parameter"));

    m.def("keplerian_to_cartesian",
          py::overload_cast<
          double, double, double, double, double, double, double >(
              &toec::convertKeplerianToCartesianElements< double > ),
          py::arg("semi_major_axis"),
          py::arg("eccentricity"),
          py::arg("inclination"),
          py::arg("argument_of_periapsis"),
          py::arg("longitude_of_ascending_node"),
          py::arg("true_anomaly"),
          py::arg("gravitational_parameter"));

    m.def("spherical_to_cartesian",
          py::overload_cast< const Eigen::Vector6d& >(
              &toec::convertSphericalOrbitalToCartesianState< double > ),
          py::arg("spherical_orbital_state"));

    m.def("spherical_to_cartesian",
          py::overload_cast<
          double, double, double, double, double, double >(
              &toec::convertSphericalOrbitalToCartesianState< double > ),
          py::arg("radial_distance"),
          py::arg("latitude"),
          py::arg("longitude"),
          py::arg("speed"),
          py::arg("flight_path_angle"),
          py::arg("heading_angle"));

};

}// namespace tudatpy
