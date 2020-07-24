/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_numerical_integrators.h"

#include <tudat/math/integrators.h>

#include <pybind11/pybind11.h>

namespace tni = tudat::numerical_integrators;
namespace py = pybind11;

namespace tudatpy {

void expose_numerical_integrators(py::module &m) {

  py::enum_<tni::AvailableIntegrators>(m, "AvailableIntegrators")
      .value("euler", tni::AvailableIntegrators::euler)
      .value("rungeKutta4", tni::AvailableIntegrators::rungeKutta4)
      .value("rungeKuttaVariableStepSize", tni::AvailableIntegrators::rungeKuttaVariableStepSize)
      .value("bulirschStoer", tni::AvailableIntegrators::bulirschStoer)
      .value("adamsBashforthMoulton", tni::AvailableIntegrators::adamsBashforthMoulton)
      .export_values();

  py::class_<
      tni::IntegratorSettings<double>,
      std::shared_ptr<tni::IntegratorSettings<double>>>(m, "IntegratorSettings")
      .def(py::init<
               const tni::AvailableIntegrators,
               const double,
               const double,
               const int,
               const bool>(),
           py::arg("integrator_type"),
           py::arg("initial_time"),
           py::arg("initial_time_step"),
           py::arg("save_frequency") = 1,
            // TODO: Discuss length of this argument: assess_propagation_termination_condition_during_integration_substeps.
           py::arg("assess_propagation_termination_condition_during_integration_substeps") = false);
}

}// namespace tudatpy