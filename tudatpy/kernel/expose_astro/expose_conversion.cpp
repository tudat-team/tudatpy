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

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace toec = tudat::orbital_element_conversions;

namespace tudatpy {

void expose_conversion(py::module &m) {

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

  m.def("convert_keplerian_to_cartesian_elements",
        &toec::convertKeplerianToCartesianElements<>,
        py::arg("keplerian_elements"),
        py::arg("central_body_gravitational_parameter"));

  m.def("convert_spherical_orbital_to_cartesian_state",
        &toec::convertSphericalOrbitalToCartesianState<>,
        py::arg("spherical_orbital_state"));

};

}// namespace tudatpy
