/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_shape_based_thrust.h"

#include <tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h>
#include <tudat/astro/low_thrust/shape_based/getRecommendedBaseFunctionsHodographicShaping.h>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;
namespace tsbm = tudat::shape_based_methods;

namespace tudatpy {
namespace astro {
namespace shape_based_thrust {


void expose_shape_based_thrust(py::module &m)
{
    py::class_<
            tsbm::BaseFunctionHodographicShaping,
            std::shared_ptr<tsbm::BaseFunctionHodographicShaping>
            >(m, "BaseFunctionHodographicShaping");

    m.def("recommended_radial_hodograph_functions",
          py::overload_cast< const double >(
              &tsbm::getRecommendedRadialVelocityBaseFunctions ),
          py::arg("time_of_flight") );

    m.def("recommended_normal_hodograph_functions",
          py::overload_cast< const double >(
              &tsbm::getRecommendedNormalBaseFunctions ),
          py::arg("time_of_flight") );


    m.def("recommended_axial_hodograph_functions",
          py::overload_cast< const double, const int >(
              &tsbm::getRecommendedAxialVelocityBaseFunctions ),
          py::arg("time_of_flight"),
          py::arg("number_of_revolutions") );


    m.def("hodograph_constant",
          &tsbm::hodographConstant );

    m.def("hodograph_sine",
          &tsbm::hodographSine,
          py::arg("frequency") );

    m.def("hodograph_cosine",
          &tsbm::hodographCosine,
          py::arg("frequency") );

    m.def("hodograph_exponential",
          &tsbm::hodographExponential,
          py::arg("exponent") );

    m.def("hodograph_scaled_exponential",
          &tsbm::hodographScaledExponential,
          py::arg("exponent"),
          py::arg("scale_factor"));

    m.def("hodograph_exponential_sine",
          &tsbm::hodographExponentialSine,
          py::arg("exponent"),
          py::arg("frequency") );

    m.def("hodograph_scaled_exponential_sine",
          &tsbm::hodographScaledExponentialSine,
          py::arg("exponent"),
          py::arg("frequency"),
          py::arg("scale_factor") );

    m.def("hodograph_exponential_cosine",
          &tsbm::hodographExponentialCosine,
          py::arg("exponent"),
          py::arg("frequency"));

    m.def("hodograph_scaled_exponential_cosine",
          &tsbm::hodographScaledExponentialCosine,
          py::arg("exponent"),
          py::arg("frequency"),
          py::arg("scale_factor") );

    m.def("hodograph_power",
          &tsbm::hodographPower,
          py::arg("exponent") );

    m.def("hodograph_scaled_power",
          &tsbm::hodographScaledPower,
          py::arg("exponent"),
          py::arg("scale_factor") );


    m.def("hodograph_power_sine",
          &tsbm::hodographPowerSine,
          py::arg("exponent"),
          py::arg("frequency") );

    m.def("hodograph_scaled_power_sine",
          &tsbm::hodographScaledPowerSine,
          py::arg("exponent"),
          py::arg("frequency"),
          py::arg("scale_factor") );

    m.def("hodograph_power_cosine",
          &tsbm::hodographPowerCosine,
          py::arg("exponent"),
          py::arg("frequency") );

    m.def("hodograph_scaled_power_cosine",
          &tsbm::hodographScaledPowerCosine,
          py::arg("exponent"),
          py::arg("frequency"),
          py::arg("scale_factor") );

}


}// namespace shape_based_thrust
}// namespace astro
}// namespace tudatpy
