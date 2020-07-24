/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_propagators.h"

#include <tudat/astro/propagators.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy {

void expose_propagators(py::module &m) {

  // Astrodynamics/Propagators/singeStateTypeDerivative.h
  m.def("get_single_integration_size",
        &tp::getSingleIntegrationSize,
        py::arg("state_type"));
  m.def("get_single_integration_differential_equation_order",
        &tp::getSingleIntegrationDifferentialEquationOrder,
        py::arg("state_type"));
  m.def("get_generalized_acceleration_size",
        &tp::getGeneralizedAccelerationSize,
        py::arg("state_type"));

  py::class_<
      tp::SingleStateTypeDerivative<double, double>,
      std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>>
      SingleStateTypeDerivative_(m, "SingleStateTypeDerivative");
  //                    .def(py::init<const tp::IntegratedStateType>(),
  //                         py::arg("integrated_state_type"));

  py::class_<
      tp::NBodyStateDerivative<double, double>,
      std::shared_ptr<tp::NBodyStateDerivative<double, double>>>
      NBodyStateDerivative_(m, "NBodyStateDerivative");

  py::class_<
      tp::NBodyCowellStateDerivative<double, double>,
      std::shared_ptr<tp::NBodyCowellStateDerivative<double, double>>>(m, "NBodyCowellStateDerivative")
      .def(py::init<
               const tudat::basic_astrodynamics::AccelerationMap &,
               const std::shared_ptr<tp::CentralBodyData<double, double>>,
               const std::vector<std::string> &>(),
           py::arg("acceleration_models_per_body"),
           py::arg("central_body_data"),
           py::arg("bodies_to_integrate"));

}

}// namespace tudatpy
