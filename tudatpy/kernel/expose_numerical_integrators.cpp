//
// Created by ggarrett on 26-04-20.
//

#include "expose_numerical_integrators.h"

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

        py::class_ <
        tni::IntegratorSettings < double > ,
                std::shared_ptr < tni::IntegratorSettings < double >> > (m, "IntegratorSettings")
                        .def(py::init<
                                     const tni::AvailableIntegrators,
                                     const double,
                                     const double,
                                     const int,
                                     const bool
                             >(),
                             py::arg("integrator_type"),
                             py::arg("initial_time"),
                             py::arg("initial_time_step"),
                             py::arg("save_frequency") = 1,
                             py::arg("assess_propagation_termination_condition_during_integration_substeps") = false);

    }
}