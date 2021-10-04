/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_integrator_setup.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/propagation_setup.h>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;


namespace tudatpy {
namespace numerical_simulation {
namespace propagation_setup {
namespace integrator {

    void expose_integrator_setup(py::module &m) {

// ENUMS
            py::enum_<tni::AvailableIntegrators>(m, "AvailableIntegrators", get_docstring("AvailableIntegrators").c_str())
                    .value("euler_type", tni::AvailableIntegrators::euler,
                           get_docstring("AvailableIntegrators.euler").c_str())
                    .value("runge_kutta_4_type", tni::AvailableIntegrators::rungeKutta4,
                           get_docstring("AvailableIntegrators.runge_kutta_4").c_str())
                    .value("runge_kutta_variable_step_size_type", tni::AvailableIntegrators::rungeKuttaVariableStepSize,
                           get_docstring("AvailableIntegrators.runge_kutta_variable_step_size").c_str())
                    .value("bulirsch_stoer_type", tni::AvailableIntegrators::bulirschStoer,
                           get_docstring("AvailableIntegrators.bulirsch_stoer").c_str())
                    .value("adams_bashforth_moulton_type", tni::AvailableIntegrators::adamsBashforthMoulton,
                           get_docstring("AvailableIntegrators.adams_bashforth_moulton").c_str());

            py::enum_<tni::RungeKuttaCoefficients::CoefficientSets>(m, "RKCoefficientSets",
                                                                    get_docstring("RKCoefficientSets").c_str())
                    .value("rkf_45", tni::RungeKuttaCoefficients::rungeKuttaFehlberg45)
                    .value("rkf_56", tni::RungeKuttaCoefficients::rungeKuttaFehlberg56)
                    .value("rkf_78", tni::RungeKuttaCoefficients::rungeKuttaFehlberg78)
                    .value("rkdp_87", tni::RungeKuttaCoefficients::rungeKutta87DormandPrince)
                    .export_values();

            py::enum_<tni::ExtrapolationMethodStepSequences>(m, "ExtrapolationMethodStepSequences",
                                                             get_docstring("ExtrapolationMethodStepSequences").c_str())
                    .value("bulirsch_stoer_sequence", tni::ExtrapolationMethodStepSequences::bulirsch_stoer_sequence)
                    .value("deufelhard_sequence", tni::ExtrapolationMethodStepSequences::deufelhard_sequence)
                    .export_values();

// CLASSES
            py::class_<
                    tni::IntegratorSettings<double>,
                    std::shared_ptr<tni::IntegratorSettings<double>>>(m, "IntegratorSettings",
                                                                      get_docstring("IntegratorSettings").c_str())
//                .def(py::init<
//                             const tni::AvailableIntegrators,
//                             const double,
//                             const double,
//                             const int,
//                             const bool>(),
//                     py::arg("integrator_type"),
//                     py::arg("initial_time"),
//                     py::arg("initial_time_step"),
//                     py::arg("save_frequency") = 1,
//                        // TODO: Discuss length of this argument: assess_propagation_termination_condition_during_integration_substeps.
//                     py::arg("assess_propagation_termination_condition_during_integration_substeps") = false)
                    .def_readwrite("initial_time", &tni::IntegratorSettings<double>::initialTime_);

            py::class_<tni::RungeKuttaVariableStepSizeBaseSettings<double>,
                    std::shared_ptr<tni::RungeKuttaVariableStepSizeBaseSettings<double>>,
                    tni::IntegratorSettings<double>>(m, "RungeKuttaVariableStepSizeBaseSettings",
                                                     get_docstring("RungeKuttaVariableStepSizeBaseSettings").c_str());

            py::class_<tni::RungeKuttaVariableStepSizeSettingsVectorTolerances<double>,
                    std::shared_ptr<tni::RungeKuttaVariableStepSizeSettingsVectorTolerances<double>>,
                    tni::RungeKuttaVariableStepSizeBaseSettings<double>>(m,
                                                                         "RungeKuttaVariableStepSizeSettingsVectorTolerances",
                                                                         get_docstring(
                                                                                 "RungeKuttaVariableStepSizeSettingsVectorTolerances").c_str());

            py::class_<tni::RungeKuttaVariableStepSizeSettingsScalarTolerances<double>,
                    std::shared_ptr<tni::RungeKuttaVariableStepSizeSettingsScalarTolerances<double>>,
                    tni::RungeKuttaVariableStepSizeBaseSettings<double>>(m,
                                                                         "RungeKuttaVariableStepSizeSettingsScalarTolerances",
                                                                         get_docstring(
                                                                                 "RungeKuttaVariableStepSizeSettingsScalarTolerances").c_str());

            py::class_<tni::BulirschStoerIntegratorSettings<double>,
                    std::shared_ptr<tni::BulirschStoerIntegratorSettings<double>>,
                    tni::IntegratorSettings<double>>(m, "BulirschStoerIntegratorSettings",
                                                     get_docstring("BulirschStoerIntegratorSettings").c_str());


            py::class_<tni::AdamsBashforthMoultonSettings<double>,
                    std::shared_ptr<tni::AdamsBashforthMoultonSettings<double>>,
                    tni::IntegratorSettings<double>>(m, "AdamsBashforthMoultonSettings",
                                                     get_docstring("AdamsBashforthMoultonSettings").c_str());


// FACTORY FUNCTIONS
            m.def("euler",
                  &tni::eulerSettings<double>,
                  py::arg("initial_time"),
                  py::arg("initial_time_step"),
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false,
                  get_docstring("euler").c_str());

            m.def("runge_kutta_4",
                  &tni::rungeKutta4Settings<double>,
                  py::arg("initial_time"),
                  py::arg("initial_time_step"),
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false,
                  get_docstring("runge_kutta_4").c_str());

            m.def("runge_kutta_variable_step_size",
                  &tni::rungeKuttaVariableStepSettingsScalarTolerances<double>,
                  py::arg("initial_time"),
                  py::arg("initial_time_step"),
                  py::arg("coefficient_set"),
                  py::arg("minimum_step_size"),
                  py::arg("maximum_step_size"),
                  py::arg("relative_error_tolerance"),
                  py::arg("absolute_error_tolerance"),
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false,
                  py::arg("safety_factor") = 0.8,
                  py::arg("maximum_factor_increase") = 4.0,
                  py::arg("minimum_factor_increase") = 0.1,
                  get_docstring("runge_kutta_variable_step_size").c_str());

            m.def("runge_kutta_variable_step_size_vector_tolerances",
                  &tni::rungeKuttaVariableStepSettingsVectorTolerances<double>,
                  py::arg("initial_time"),
                  py::arg("initial_time_step"),esa
                  py::arg("coefficient_set"),
                  py::arg("minimum_step_size"),
                  py::arg("maximum_step_size"),
                  py::arg("relative_error_tolerance"),
                  py::arg("absolute_error_tolerance"),
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false,
                  py::arg("safety_factor") = 0.8,
                  py::arg("maximum_factor_increase") = 4.0,
                  py::arg("minimum_factor_increase") = 0.1,
                  get_docstring("runge_kutta_variable_step_size_vector_tolerances").c_str());

            m.def("bulirsch_stoer",
                  &tni::bulirschStoerIntegratorSettings<double>,
                  py::arg("initial_time"),
                  py::arg("initial_time_step"),
                  py::arg("extrapolation_sequence"),
                  py::arg("maximum_number_of_steps"),
                  py::arg("minimum_step_size"),
                  py::arg("maximum_step_size"),
                  py::arg("relative_error_tolerance") = 1.0E-12,
                  py::arg("absolute_error_tolerance") = 1.0E-12,
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false,
                  py::arg("safety_factor") = 0.7,
                  py::arg("maximum_factor_increase") = 10.0,
                  py::arg("minimum_factor_increase") = 0.1,
                  get_docstring("bulirsch_stoer").c_str());

            m.def("adams_bashforth_moulton",
                  &tni::adamsBashforthMoultonSettings<double>,
                  py::arg("initial_time"),
                  py::arg("initial_time_step"),
                  py::arg("minimum_step_size"),
                  py::arg("maximum_step_size"),
                  py::arg("relative_error_tolerance") = 1.0E-12,
                  py::arg("absolute_error_tolerance") = 1.0E-12,
                  py::arg("minimum_order") = 6,
                  py::arg("maximum_order") = 11,
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false,
                  py::arg("bandwidth") = 200.0,
                  get_docstring("adams_bashforth_moulton").c_str());
        }

}// namespace integrator
}// namespace propagation_setup
}// namespace numerical_simulation
}// namespace tudatpy
