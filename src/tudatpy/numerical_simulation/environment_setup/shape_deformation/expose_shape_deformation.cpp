/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/simulation/environment_setup/createBodyDeformationModel.h>

#include "tudatpy/docstrings.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {
            namespace shape_deformation {

                PYBIND11_MODULE(expose_shape_deformation, m) {
                    py::class_<tss::BodyDeformationSettings,
                               std::shared_ptr<tss::BodyDeformationSettings>>(
                        m, "BodyDeformationSettings",
                        get_docstring("BodyDeformationSettings").c_str());

                    py::class_<
                        tss::BasicSolidBodyDeformationSettings,
                        std::shared_ptr<tss::BasicSolidBodyDeformationSettings>,
                        tss::BodyDeformationSettings>(
                        m, "BasicSolidBodyDeformationSettings",
                        get_docstring("BasicSolidBodyDeformationSettings")
                            .c_str());


                    m.def("basic_solid_body_tidal",
                          &tss::basicTidalBodyShapeDeformation,
                          py::arg("tide_raising_bodies"),
                          py::arg("displacement_love_numbers"),
                          py::arg("reference_radius") = TUDAT_NAN,
                          get_docstring("basic_solid_body_tidal").c_str());

                    m.def("degree_two_basic_solid_body_tidal",
                          &tss::degreeTwoBasicTidalBodyShapeDeformation,
                          py::arg("tide_raising_bodies"),
                          py::arg("love_number"), py::arg("shida_number"),
                          py::arg("reference_radius") = TUDAT_NAN,
                          get_docstring("degree_two_basic_solid_body_tidal")
                              .c_str());

                    m.def("iers_2010_solid_body_tidal",
                          &tss::iers2010TidalBodyShapeDeformation,
                          get_docstring("iers_2010_solid_body_tidal").c_str());
                }

            }  // namespace shape_deformation
        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
