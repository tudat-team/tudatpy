/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/gravitation.h>
#include <tudat/math/basic.h>

#include "tudatpy/docstrings.h"

namespace py = pybind11;
namespace tg = tudat::gravitation;
namespace tbm = tudat::basic_mathematics;

namespace tudat {
    namespace gravitation {
        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, double>
        getDegreeTwoSphericalHarmonicCoefficientsPy(
            const Eigen::Matrix3d inertiaTensor,
            const double bodyGravitationalParameter,
            const double referenceRadius,
            const bool useNormalizedCoefficients) {
            return tg::getDegreeTwoSphericalHarmonicCoefficients(
                inertiaTensor, bodyGravitationalParameter, referenceRadius, 2,
                useNormalizedCoefficients);
        }
    }  // namespace gravitation
}  // namespace tudat

namespace tudatpy {
    namespace astro {
        namespace gravitation {


            PYBIND11_MODULE(expose_gravitation, m) {
                m.def("legendre_normalization_factor",
                      &tbm::calculateLegendreGeodesyNormalizationFactor,
                      py::arg("degree"), py::arg("order"),
                      get_docstring("legendre_normalization_factor").c_str());

                m.def(
                    "normalize_spherical_harmonic_coefficients",
                    py::overload_cast<const Eigen::MatrixXd&,
                                      const Eigen::MatrixXd&>(
                        &tbm::
                            convertUnnormalizedToGeodesyNormalizedCoefficients),
                    py::arg("unnormalized_cosine_coefficients"),
                    py::arg("unnormalized_sine_coefficients"),
                    get_docstring("normalize_spherical_harmonic_coefficients")
                        .c_str());

                m.def(
                    "unnormalize_spherical_harmonic_coefficients",
                    py::overload_cast<const Eigen::MatrixXd&,
                                      const Eigen::MatrixXd&>(
                        &tbm::
                            convertGeodesyNormalizedToUnnormalizedCoefficients),
                    py::arg("normalized_cosine_coefficients"),
                    py::arg("normalized_sine_coefficients"),
                    get_docstring("unnormalize_spherical_harmonic_coefficients")
                        .c_str());

                m.def("spherical_harmonic_coefficients_from_inertia",
                      tg::getDegreeTwoSphericalHarmonicCoefficientsPy,
                      py::arg("inertia_tensor"),
                      py::arg("gravitational_parameter"),
                      py::arg("reference_radius"),
                      py::arg("output_normalized_coefficients") = true,
                      get_docstring(
                          "spherical_harmonic_coefficients_from_inertia")
                          .c_str());
            }


        }  // namespace gravitation
    }  // namespace astro
}  // namespace tudatpy
