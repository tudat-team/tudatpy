/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_gravity_field_variation_setup.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/environment_setup.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

//#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
//#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tg = tudat::gravitation;

namespace tudatpy {
namespace numerical_simulation {
namespace environment_setup {
namespace gravity_field_variation {

    void expose_gravity_field_variation_setup(py::module &m) {
        py::enum_<tg::BodyDeformationTypes>(
                m, "BodyDeformationTypes", "<no_doc>")
                .value("basic_solid_body", tg::basic_solid_body)
                .value("tabulated_deformation", tg::tabulated_variation)
                .export_values();

        py::class_<tss::GravityFieldVariationSettings,
                std::shared_ptr<tss::GravityFieldVariationSettings>>(m, "GravityFieldVariationSettings");


// fixedSingleDegreeLoveNumberGravityFieldVariationSettingsSimplified (originally exposed function) not declared
// if followed next best declaration, it results in duplicate exposition (compare underneath)
//        m.def("solid_body_tide_simplified",
//              py::overload_cast<const std::string, const double, const int>(
//                    &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
//              py::arg("tide_raising_body"),
//              py::arg("love_number"),
//              py::arg("degree"));


        m.def("solid_body_tide",
              py::overload_cast<const std::string, const double, const int>(
                      &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
              py::arg("tide_raising_body"),
              py::arg("love_number"),
              py::arg("degree"));

        m.def("solid_body_tide",
              py::overload_cast<const std::string, const std::complex<double>, const int>(
                      &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
              py::arg("tide_raising_body"),
              py::arg("love_number"),
              py::arg("degree"));

        m.def("solid_body_tide",
              py::overload_cast<const std::string, std::map<int, double> >(
                      &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
              py::arg("tide_raising_body"),
              py::arg("love_number_per_degree"));

        m.def("solid_body_tide",
              py::overload_cast<const std::string, std::map<int, std::complex<double> > >(
                      &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
              py::arg("tide_raising_body"),
              py::arg("love_number_per_degree"));

        m.def("solid_body_tide",
              py::overload_cast<const std::string, const std::vector<double>, const int,
                      const std::shared_ptr<tss::ModelInterpolationSettings> >(
                      &tss::orderVariableSingleDegreeLoveNumberGravityFieldVariationSettings),
              py::arg("tide_raising_body"),
              py::arg("love_number_per_order"),
              py::arg("degree"),
              py::arg("interpolation_settings") = nullptr);

        m.def("solid_body_tide",
              py::overload_cast<const std::string, const std::vector<std::complex<double> >, const int,
                      const std::shared_ptr<tss::ModelInterpolationSettings> >(
                      &tss::orderVariableSingleDegreeLoveNumberGravityFieldVariationSettings),
              py::arg("tide_raising_body"),
              py::arg("love_number_per_order"),
              py::arg("degree"),
              py::arg("interpolation_settings") = nullptr);

        m.def("solid_body_tide",
              py::overload_cast<const std::string, const std::map<int, std::vector<double> >,
                      const std::shared_ptr<tss::ModelInterpolationSettings> >(
                      &tss::degreeOrderVariableLoveNumberGravityFieldVariationSettings),
              py::arg("tide_raising_body"),
              py::arg("love_number_per_degree_and_order"),
              py::arg("interpolation_settings") = nullptr);

        m.def("solid_body_tide",
              py::overload_cast<const std::string, const std::map<int, std::vector<std::complex<double> > >,
                      const std::shared_ptr<tss::ModelInterpolationSettings> >(
                      &tss::degreeOrderVariableLoveNumberGravityFieldVariationSettings),
              py::arg("tide_raising_body"),
              py::arg("love_number_per_degree_and_order"),
              py::arg("interpolation_settings") = nullptr);

        m.def("tabulated",
              &tss::tabulatedGravityFieldVariationSettings,
              py::arg("cosine_variations_table"),
              py::arg("sine_variations_table"),
              py::arg("minimum_degree"),
              py::arg("minimum_order"),
              py::arg("interpolation_settings"));
    }

}// namespace gravity_field_variation
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy