/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/basic_astro.h>

namespace py = pybind11;
namespace tmg = tudat::mission_geometry;

namespace tudatpy {
    namespace astro {
        namespace fundamentals {

            PYBIND11_MODULE(expose_fundamentals, m) {
                m.def("compute_shadow_function", &tmg::computeShadowFunction,
                      py::arg("occulted_body_position"),
                      py::arg("occulted_body_radius"),
                      py::arg("occulting_body_position"),
                      py::arg("occulting_body_radius"),
                      py::arg("satellite_position"));
            }

        }  // namespace fundamentals
    }  // namespace astro
}  // namespace tudatpy
