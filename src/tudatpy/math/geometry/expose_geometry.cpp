/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

// #include "expose_geometry.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/basics/basicTypedefs.h>
#include <tudat/math/geometric.h>

namespace py = pybind11;

namespace tgs = tudat::geometric_shapes;

namespace tudatpy {

    PYBIND11_MODULE(expose_geometry, m) {
        py::class_<tudat::SurfaceGeometry,
                   std::shared_ptr<tudat::SurfaceGeometry>>(m,
                                                            "SurfaceGeometry");

        py::class_<tgs::CompositeSurfaceGeometry,
                   std::shared_ptr<tgs::CompositeSurfaceGeometry>,
                   tudat::SurfaceGeometry>(m, "CompositeSurfaceGeometry");

        py::class_<tgs::Capsule, std::shared_ptr<tgs::Capsule>,
                   tgs::CompositeSurfaceGeometry>(m, "Capsule")
            .def(py::init<const double, const double, const double,
                          const double, const double>(),
                 py::arg("nose_radius"), py::arg("middle_radius"),
                 py::arg("rear_length"), py::arg("rear_angle"),
                 py::arg("side_radius"))
            .def_property_readonly("middle_radius",
                                   &tgs::Capsule::getMiddleRadius)
            .def_property_readonly("volume", &tgs::Capsule::getVolume)
            .def_property_readonly("length", &tgs::Capsule::getLength);
    };

}  // namespace tudatpy
