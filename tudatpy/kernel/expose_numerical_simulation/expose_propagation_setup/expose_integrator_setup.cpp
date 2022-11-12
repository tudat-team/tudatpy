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
              //       .value("euler_type", tni::AvailableIntegrators::euler)
              //       .value("runge_kutta_4_type", tni::AvailableIntegrators::rungeKutta4)
                    .value("runge_kutta_fixed_step_size_type", tni::AvailableIntegrators::rungeKuttaFixedStepSize)
                    .value("runge_kutta_variable_step_size_type", tni::AvailableIntegrators::rungeKuttaVariableStepSize)
                    .value("bulirsch_stoer_type", tni::AvailableIntegrators::bulirschStoer)
                    .value("adams_bashforth_moulton_type", tni::AvailableIntegrators::adamsBashforthMoulton)
                    .export_values();

            py::enum_<tni::CoefficientSets>(m, "CoefficientSets",
                                            get_docstring("CoefficientSets").c_str())
                    .value("euler_forward", tni::forwardEuler,
                           get_docstring("CoefficientSets.euler_forward").c_str())
                    .value("rk_4", tni::rungeKutta4Classic,
                           get_docstring("CoefficientSets.rk_4").c_str())
                    .value("explicit_mid_point", tni::explicitMidPoint,
                           get_docstring("CoefficientSets.explicit_mid_point").c_str())
                    .value("explicit_trapezoid_rule", tni::explicitTrapezoidRule,
                           get_docstring("CoefficientSets.explicit_trapezoid_rule").c_str())
                    .value("ralston", tni::ralston,
                           get_docstring("CoefficientSets.ralston").c_str())
                    .value("rk_3", tni::rungeKutta3,
                           get_docstring("CoefficientSets.rk_3").c_str())
                    .value("ralston_3", tni::ralston3,
                           get_docstring("CoefficientSets.ralston_3").c_str())
                    .value("SSPRK3", tni::SSPRK3,
                           get_docstring("CoefficientSets.SSPRK3").c_str())
                    .value("ralston_4", tni::ralston4,
                           get_docstring("CoefficientSets.ralston_4").c_str())
                    .value("three_eight_rule_rk_4", tni::threeEighthRuleRK4,
                           get_docstring("CoefficientSets.three_eight_rule_rk_4").c_str())
                    .value("heun_euler", tni::heunEuler,
                           get_docstring("CoefficientSets.heun_euler").c_str())
                    .value("rkf_12", tni::rungeKuttaFehlberg12,
                           get_docstring("CoefficientSets.rkf_12").c_str())
                    .value("rkf_45", tni::rungeKuttaFehlberg45,
                           get_docstring("CoefficientSets.rkf_45").c_str())
                    .value("rkf_56", tni::rungeKuttaFehlberg56,
                           get_docstring("CoefficientSets.rkf_56").c_str())
                    .value("rkf_78", tni::rungeKuttaFehlberg78,
                           get_docstring("CoefficientSets.rkf_78").c_str())
                    .value("rkdp_87", tni::rungeKutta87DormandPrince,
                           get_docstring("CoefficientSets.rkdp_87").c_str())
                    .value("rkf_89", tni::rungeKuttaFehlberg89,
                           get_docstring("CoefficientSets.rkf_89").c_str())
                    .value("rkv_89", tni::rungeKuttaVerner89,
                           get_docstring("CoefficientSets.rkv_89").c_str())
                    .value("rkf_108", tni::rungeKuttaFeagin108,
                           get_docstring("CoefficientSets.rkf_108").c_str())
                    .value("rkf_1210", tni::rungeKuttaFeagin1210,
                           get_docstring("CoefficientSets.rkf_1210").c_str())
                    .value("rkf_1412", tni::rungeKuttaFeagin1412,
                           get_docstring("CoefficientSets.rkf_1412").c_str())
                    .export_values();

            py::enum_<tni::RungeKuttaCoefficients::OrderEstimateToIntegrate>(m, "OrderToIntegrate",
                      get_docstring("OrderToIntegrate").c_str())
                     .value("lower", tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::lower,
                           get_docstring("OrderToIntegrate.lower").c_str())
                     .value("higher", tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::higher,
                           get_docstring("OrderToIntegrate.higher").c_str())
                     .export_values();

            py::enum_<tni::ExtrapolationMethodStepSequences>(m, "ExtrapolationMethodStepSequences",
                                                             get_docstring("ExtrapolationMethodStepSequences").c_str())
                    .value("bulirsch_stoer_sequence", tni::ExtrapolationMethodStepSequences::bulirsch_stoer_sequence,
                           get_docstring("ExtrapolationMethodStepSequences.bulirsch_stoer_sequence").c_str())
                    .value("deufelhard_sequence", tni::ExtrapolationMethodStepSequences::deufelhard_sequence,
                           get_docstring("ExtrapolationMethodStepSequences.deufelhard_sequence").c_str())
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
                    .def_readwrite("initial_time", &tni::IntegratorSettings<double>::initialTimeDeprecated_);

            py::class_<tni::RungeKuttaFixedStepSizeSettings<double>,
                    std::shared_ptr<tni::RungeKuttaFixedStepSizeSettings<double>>,
                    tni::IntegratorSettings<double>>(m, "RungeKuttaFixedStepSizeSettings",
                                                     get_docstring("RungeKuttaFixedStepSizeSettings").c_str());

            py::class_<tni::RungeKuttaVariableStepSizeBaseSettings<double>,
                    std::shared_ptr<tni::RungeKuttaVariableStepSizeBaseSettings<double>>,
                    tni::IntegratorSettings<double>>(m, "RungeKuttaVariableStepSizeBaseSettings",
                                                     get_docstring("RungeKuttaVariableStepSizeBaseSettings").c_str());

            py::class_<tni::RungeKuttaVariableStepSizeSettingsVectorTolerances<double>,
                    std::shared_ptr<tni::RungeKuttaVariableStepSizeSettingsVectorTolerances<double>>,
                    tni::RungeKuttaVariableStepSizeBaseSettings<double>>(m,
                                                                         "RungeKuttaVariableStepSizeSettingsVectorTolerances",
                                                                         get_docstring("RungeKuttaVariableStepSizeSettingsVectorTolerances").c_str());

            py::class_<tni::RungeKuttaVariableStepSizeSettingsScalarTolerances<double>,
                    std::shared_ptr<tni::RungeKuttaVariableStepSizeSettingsScalarTolerances<double>>,
                    tni::RungeKuttaVariableStepSizeBaseSettings<double>>(m,
                                                                         "RungeKuttaVariableStepSizeSettingsScalarTolerances",
                                                                         get_docstring("RungeKuttaVariableStepSizeSettingsScalarTolerances").c_str());

            py::class_<tni::BulirschStoerIntegratorSettings<double>,
                    std::shared_ptr<tni::BulirschStoerIntegratorSettings<double>>,
                    tni::IntegratorSettings<double>>(m, "BulirschStoerIntegratorSettings",
                                                     get_docstring("BulirschStoerIntegratorSettings").c_str());


            py::class_<tni::AdamsBashforthMoultonSettings<double>,
                    std::shared_ptr<tni::AdamsBashforthMoultonSettings<double>>,
                    tni::IntegratorSettings<double>>(m, "AdamsBashforthMoultonSettings",
                                                     get_docstring("AdamsBashforthMoultonSettings").c_str());


// FACTORY FUNCTIONS
            m.def("print_butcher_tableau",
                  &tni::printButcherTableau,
                  get_docstring("print_butcher_tableau").c_str());

            m.def("euler",
                  &tni::eulerSettingsDeprecated<double>,
                  py::arg("initial_time"),
                  py::arg("initial_time_step"),
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false );

            m.def("euler",
                  &tni::eulerSettings<double>,
                  py::arg("initial_time_step"),
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false,
                  get_docstring("euler").c_str());


            m.def("runge_kutta_4",
                  &tni::rungeKutta4SettingsDeprecated<double>,
                  py::arg("initial_time"),
                  py::arg("initial_time_step"),
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false);

            m.def("runge_kutta_4",
                  &tni::rungeKutta4Settings<double>,
                  py::arg("initial_time_step"),
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false,
                  get_docstring("runge_kutta_4").c_str());

            m.def("runge_kutta_fixed_step_size",
                  &tni::rungeKuttaFixedStepSettingsDeprecated<double>,
                  py::arg("initial_time"),
                  py::arg("initial_time_step"),
                  py::arg("coefficient_set"),
                  py::arg("order_to_use") = tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::lower,
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false );

            m.def("runge_kutta_fixed_step_size",
                  &tni::rungeKuttaFixedStepSettings<double>,
                  py::arg("initial_time_step"),
                  py::arg("coefficient_set"),
                  py::arg("order_to_use") = tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::lower,
                  py::arg("save_frequency") = 1,
                  py::arg("assess_termination_on_minor_steps") = false,
                  get_docstring("runge_kutta_fixed_step_size").c_str());

            m.def("runge_kutta_variable_step_size",
                  &tni::rungeKuttaVariableStepSettingsScalarTolerancesDeprecated<double>,
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
                  py::arg("throw_exception_if_minimum_step_exceeded") = true);

            m.def("runge_kutta_variable_step_size",
                  &tni::rungeKuttaVariableStepSettingsScalarTolerances<double>,
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
                  py::arg("throw_exception_if_minimum_step_exceeded") = true,
                  get_docstring("runge_kutta_variable_step_size").c_str());

            m.def("runge_kutta_variable_step_size_vector_tolerances",
                  &tni::rungeKuttaVariableStepSettingsVectorTolerancesDeprecated<double>,
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
                  py::arg("throw_exception_if_minimum_step_exceeded") = true);

            m.def("runge_kutta_variable_step_size_vector_tolerances",
                  &tni::rungeKuttaVariableStepSettingsVectorTolerances<double>,
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
                  py::arg("throw_exception_if_minimum_step_exceeded") = true,
                  get_docstring("runge_kutta_variable_step_size_vector_tolerances").c_str());

            m.def("bulirsch_stoer",
                  &tni::bulirschStoerIntegratorSettingsDeprecated<double>,
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
                  py::arg("minimum_factor_increase") = 0.1);

            m.def("bulirsch_stoer",
                  &tni::bulirschStoerIntegratorSettings<double>,
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
                  &tni::adamsBashforthMoultonSettingsDeprecated<double>,
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
                  py::arg("bandwidth") = 200.0);

            m.def("adams_bashforth_moulton",
                  &tni::adamsBashforthMoultonSettings<double>,
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
