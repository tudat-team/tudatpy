/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "docstrings.h"

#include "expose_polyhedron_utilities.h"

#include <tudat/astro/basic_astro/polyhedronFuntions.h>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

namespace tba = tudat::basic_astrodynamics;

namespace py = pybind11;

namespace tudatpy {

namespace astro {
namespace polyhedron_utilities {

void expose_polyhedron_utilities(py::module &m) {

    m.def("surface_area",
          &tba::computePolyhedronSurfaceArea,
          py::arg("vertices_coordinates"),
          py::arg("vertices_defining_each_facet"),
          get_docstring("surface_area").c_str());

    m.def("volume",
          &tba::computePolyhedronVolume,
          py::arg("vertices_coordinates"),
          py::arg("vertices_defining_each_facet"),
          get_docstring("volume").c_str());

    m.def("centroid",
          &tba::computePolyhedronCentroidPosition,
          py::arg("vertices_coordinates"),
          py::arg("vertices_defining_each_facet"),
          get_docstring("centroid").c_str());

    m.def("modify_centroid",
          &tba::modifyPolyhedronCentroidPosition,
          py::arg("vertices_coordinates"),
          py::arg("vertices_defining_each_facet"),
          py::arg("desired_centroid"),
          get_docstring("modify_centroid").c_str());

    m.def("inertia_tensor_from_density",
          py::overload_cast<
          const Eigen::MatrixXd&,
          const Eigen::MatrixXi&,
          const double > (&tba::computePolyhedronInertiaTensor),
          py::arg("vertices_coordinates"),
          py::arg("vertices_defining_each_facet"),
          py::arg("density"),
          get_docstring("inertia_tensor_from_density").c_str());

    m.def("inertia_tensor_from_gravitational_parameter",
          py::overload_cast<
          const Eigen::MatrixXd&,
          const Eigen::MatrixXi&,
          const double,
          const double > (&tba::computePolyhedronInertiaTensor),
          py::arg("vertices_coordinates"),
          py::arg("vertices_defining_each_facet"),
          py::arg("gravitational_parameter"),
          py::arg("gravitational_constant"),
          get_docstring("inertia_tensor_from_gravitational_parameter").c_str());

}

} // namespace polyhedron_utilities
} // namespace astro
}// namespace tudatpy
