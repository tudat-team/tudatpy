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
#include <tudat/math/basic.h>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace tg = tudat::gravitation;
namespace tbm = tudat::basic_mathematics;

namespace tudatpy {
<<<<<<< HEAD
namespace astro {
namespace gravitation {
=======
    namespace astro {
        namespace gravitation {
>>>>>>> 430786d34e5d83dfba61cf73c209ba1eca6b921e

void expose_gravitation(py::module &m) {

    m.def("legendre_normalization_factor",
          &tbm::calculateLegendreGeodesyNormalizationFactor,
          py::arg("degree"),
          py::arg("order") );

    m.def("normalize_spherical_harmonic_coefficients",
          py::overload_cast< const Eigen::MatrixXd&, const Eigen::MatrixXd& >(
              &tbm::convertUnnormalizedToGeodesyNormalizedCoefficients ),
          py::arg("unnormalized_cosine_coefficients"),
          py::arg("unnormalized_sine_coefficients") );

    m.def("unnormalize_spherical_harmonic_coefficients",
          py::overload_cast< const Eigen::MatrixXd&, const Eigen::MatrixXd& >(
              &tbm::convertGeodesyNormalizedToUnnormalizedCoefficients ),
          py::arg("normalized_cosine_coefficients"),
          py::arg("normalized_sine_coefficients") );

    m.def("spherical_harmonic_coefficients_from_inertia",
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
}

<<<<<<< HEAD
}
}
=======
        } // namespace gravitation
    } // namespace astro
>>>>>>> 430786d34e5d83dfba61cf73c209ba1eca6b921e
}// namespace tudatpy
