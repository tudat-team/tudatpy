/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_gravitation.h"

#include <tudat/astro/gravitation.h>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
namespace tg = tudat::gravitation;

namespace tudatpy {

void expose_gravitation(py::module &m) {

    py::class_<tg::GravityFieldModel,
            std::shared_ptr<tg::GravityFieldModel>>(m, "GravityFieldModel")
            .def(py::init<
                 const double,
                 const std::function<void()>>(),
                 py::arg("gravitational_parameter"),
                 py::arg("update_inertia_tensor") = std::function<void()>()// <pybind11/functional.h>
            )
            .def("get_gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter)
            .def_property("gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter,
                          &tg::GravityFieldModel::resetGravitationalParameter);

    py::class_<tg::SphericalHarmonicsGravityField,
            std::shared_ptr<tg::SphericalHarmonicsGravityField >,
            tg::GravityFieldModel>(m, "SphericalHarmonicsGravityField")
            .def_property_readonly("reference_radius", &tg::SphericalHarmonicsGravityField::getReferenceRadius )
            .def_property_readonly("maximum_degree", &tg::SphericalHarmonicsGravityField::getDegreeOfExpansion )
            .def_property_readonly("maximum_order", &tg::SphericalHarmonicsGravityField::getOrderOfExpansion )
            .def_property("cosine_coefficients", &tg::SphericalHarmonicsGravityField::getCosineCoefficients,
                          &tg::SphericalHarmonicsGravityField::setCosineCoefficients)
            .def_property("sine_coefficients", &tg::SphericalHarmonicsGravityField::getSineCoefficients,
                          &tg::SphericalHarmonicsGravityField::setSineCoefficients);

    m.def("degree_spherical_harmonic_coefficients_from_inertia",
          py::overload_cast<
          const Eigen::Matrix3d,
          const double,
          const double,
          const int,
          const bool >( &tg::getDegreeTwoSphericalHarmonicCoefficients ),
          py::arg("inertia_tensor"),
          py::arg("gravitational_parameter"),
          py::arg("reference_radius"),
          py::arg("maximum_output_degree") = 2,
          py::arg("output_normalized_coefficients") = true );
};

}// namespace tudatpy
