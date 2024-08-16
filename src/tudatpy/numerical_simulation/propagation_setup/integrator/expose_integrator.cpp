/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/simulation/propagation_setup.h>

#include "tudatpy/docstrings.h"
#include "tudatpy/scalarTypes.h"

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

                PYBIND11_MODULE(expose_integrator, m) {
                    // ENUMS
                    py::enum_<tni::MinimumIntegrationTimeStepHandling>(
                        m, "MinimumIntegrationTimeStepHandling",
                        get_docstring("MinimumIntegrationTimeStepHandling")
                            .c_str())
                        .value("throw_exception_below_minimum",
                               tni::MinimumIntegrationTimeStepHandling::
                                   throw_exception_below_minimum,
                               get_docstring("MinimumIntegrationTimeStepHandlin"
                                             "g.throw_exception_below_minimum")
                                   .c_str())
                        .value("set_to_minimum_step_silently",
                               tni::MinimumIntegrationTimeStepHandling::
                                   set_to_minimum_step_silently,
                               get_docstring("MinimumIntegrationTimeStepHandlin"
                                             "g.set_to_minimum_step_silently")
                                   .c_str())
                        .value(
                            "set_to_minimum_step_single_warning",
                            tni::MinimumIntegrationTimeStepHandling::
                                set_to_minimum_step_single_warning,
                            get_docstring("MinimumIntegrationTimeStepHandling."
                                          "set_to_minimum_step_single_warning")
                                .c_str())
                        .value("set_to_minimum_step_every_time_warning",
                               tni::MinimumIntegrationTimeStepHandling::
                                   set_to_minimum_step_every_time_warning,
                               get_docstring(
                                   "MinimumIntegrationTimeStepHandling.set_to_"
                                   "minimum_step_every_time_warning")
                                   .c_str())
                        .export_values();


                    py::enum_<tni::AvailableIntegrators>(
                        m, "AvailableIntegrators",
                        R"doc(Enumeration of integrators available with tudat.


	:member runge_kutta_fixed_step_size_type:
	:member runge_kutta_variable_step_size_type:
	:member bulirsch_stoer_type:
	:member adams_bashforth_moulton_type:
)doc")
                        //       .value("euler_type",
                        //       tni::AvailableIntegrators::euler)
                        //       .value("runge_kutta_4_type",
                        //       tni::AvailableIntegrators::rungeKutta4)
                        .value(
                            "runge_kutta_fixed_step_size_type",
                            tni::AvailableIntegrators::rungeKuttaFixedStepSize)
                        .value("runge_kutta_variable_step_size_type",
                               tni::AvailableIntegrators::
                                   rungeKuttaVariableStepSize)
                        .value("bulirsch_stoer_type",
                               tni::AvailableIntegrators::bulirschStoer)
                        .value("adams_bashforth_moulton_type",
                               tni::AvailableIntegrators::adamsBashforthMoulton)
                        .export_values();

                    py::enum_<tni::CoefficientSets>(
                        m, "CoefficientSets",
                        R"doc(Coefficient sets for Runge-Kutta-type integrators.

	Coefficient sets for Runge-Kutta-type integrators. The coefficients are defined
	in a Butcher Tableau, with an coefficient set yielding an x(y) method yielding an integrator
	with global truncation error of :math:`O(\Delta t^{x})`. Some of these coefficients also contain an embedded integrator of :math:`O(\Delta t^{y})`
	for step size control.


	:member euler_forward:
	:member rk_4:
	:member explicit_mid_point:
	:member explicit_trapezoid_rule:
	:member ralston:
	:member rk_3:
	:member ralston_3:
	:member SSPRK3:
	:member ralston_4:
	:member three_eight_rule_rk_4:
	:member rkf_12:
	:member heun_euler:
	:member rkf_45:
	:member rkf_56:
	:member rkf_78:
	:member rkdp_87:
	:member rkf_89:
	:member rkv_89:
	:member rkf_108:
	:member rkf_1210:
	:member rkf_1412:
)doc")
                        .value("euler_forward", tni::forwardEuler,
                               get_docstring("CoefficientSets.euler_forward")
                                   .c_str())
                        .value("rk_4", tni::rungeKutta4Classic,
                               get_docstring("CoefficientSets.rk_4").c_str())
                        .value(
                            "explicit_mid_point", tni::explicitMidPoint,
                            get_docstring("CoefficientSets.explicit_mid_point")
                                .c_str())
                        .value("explicit_trapezoid_rule",
                               tni::explicitTrapezoidRule,
                               get_docstring(
                                   "CoefficientSets.explicit_trapezoid_rule")
                                   .c_str())
                        .value("ralston", tni::ralston,
                               get_docstring("CoefficientSets.ralston").c_str())
                        .value("rk_3", tni::rungeKutta3,
                               get_docstring("CoefficientSets.rk_3").c_str())
                        .value(
                            "ralston_3", tni::ralston3,
                            get_docstring("CoefficientSets.ralston_3").c_str())
                        .value("SSPRK3", tni::SSPRK3,
                               get_docstring("CoefficientSets.SSPRK3").c_str())
                        .value(
                            "ralston_4", tni::ralston4,
                            get_docstring("CoefficientSets.ralston_4").c_str())
                        .value("three_eight_rule_rk_4", tni::threeEighthRuleRK4,
                               get_docstring(
                                   "CoefficientSets.three_eight_rule_rk_4")
                                   .c_str())
                        .value(
                            "heun_euler", tni::heunEuler,
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
                        .value(
                            "rkf_1210", tni::rungeKuttaFeagin1210,
                            get_docstring("CoefficientSets.rkf_1210").c_str())
                        .value(
                            "rkf_1412", tni::rungeKuttaFeagin1412,
                            get_docstring("CoefficientSets.rkf_1412").c_str())
                        .export_values();

                    py::enum_<
                        tni::RungeKuttaCoefficients::OrderEstimateToIntegrate>(
                        m, "OrderToIntegrate",
                        R"doc(Which integrator order needs to be integrated, only for coefficient sets with an embedded order.


	:member lower:
	:member higher:
)doc")
                        .value("lower",
                               tni::RungeKuttaCoefficients::
                                   OrderEstimateToIntegrate::lower,
                               get_docstring("OrderToIntegrate.lower").c_str())
                        .value("higher",
                               tni::RungeKuttaCoefficients::
                                   OrderEstimateToIntegrate::higher,
                               get_docstring("OrderToIntegrate.higher").c_str())
                        .export_values();

                    py::enum_<tni::ExtrapolationMethodStepSequences>(
                        m, "ExtrapolationMethodStepSequences",
                        R"doc(Enumeration of available extrapolation method step sequences.


	:member bulirsch_stoer_sequence:
	:member deufelhard_sequence:
)doc")
                        .value("bulirsch_stoer_sequence",
                               tni::ExtrapolationMethodStepSequences::
                                   bulirsch_stoer_sequence,
                               get_docstring("ExtrapolationMethodStepSequences."
                                             "bulirsch_stoer_sequence")
                                   .c_str())
                        .value("deufelhard_sequence",
                               tni::ExtrapolationMethodStepSequences::
                                   deufelhard_sequence,
                               get_docstring("ExtrapolationMethodStepSequences."
                                             "deufelhard_sequence")
                                   .c_str())
                        .export_values();

                    // CLASSES
                    py::class_<
                        tni::IntegratorSettings<TIME_TYPE>,
                        std::shared_ptr<tni::IntegratorSettings<TIME_TYPE>>>(
                        m, "IntegratorSettings",
                        R"doc(Functional base class to define settings for integrators.

	Class to define settings for numerical integrators, for instance for use in numerical integration of equations of motion/
	variational equations. This class can be used for simple integrators such as fixed step RK and Euler. Integrators that
	require more settings to define have their own derived class.

)doc")
                        .def_readwrite("initial_time",
                                       &tni::IntegratorSettings<
                                           TIME_TYPE>::initialTimeDeprecated_);

                    py::class_<
                        tni::RungeKuttaFixedStepSizeSettings<TIME_TYPE>,
                        std::shared_ptr<
                            tni::RungeKuttaFixedStepSizeSettings<TIME_TYPE>>,
                        tni::IntegratorSettings<TIME_TYPE>>(
                        m, "RungeKuttaFixedStepSizeSettings",
                        R"doc(`IntegratorSettings`-derived class to define settings for Runge Kutta integrators with a fixed step size

)doc");

                    py::class_<
                        tni::RungeKuttaVariableStepSizeBaseSettings<TIME_TYPE>,
                        std::shared_ptr<
                            tni::RungeKuttaVariableStepSizeBaseSettings<
                                TIME_TYPE>>,
                        tni::IntegratorSettings<TIME_TYPE>>(
                        m, "RungeKuttaVariableStepSizeBaseSettings",
                        get_docstring("RungeKuttaVariableStepSizeBaseSettings")
                            .c_str());

                    py::class_<
                        tni::RungeKuttaVariableStepSizeSettingsVectorTolerances<
                            TIME_TYPE>,
                        std::shared_ptr<
                            tni::
                                RungeKuttaVariableStepSizeSettingsVectorTolerances<
                                    TIME_TYPE>>,
                        tni::RungeKuttaVariableStepSizeBaseSettings<TIME_TYPE>>(
                        m, "RungeKuttaVariableStepSizeSettingsVectorTolerances",
                        R"doc(`IntegratorSettings`-derived class to define settings for Runge Kutta integrators with vector tolerances.

	This class is actually derived from an intermediate class (`RungeKuttaVariableStepSizeBaseSettings`, not documented),
	which is derived directly from `IntegratorSettings`.

)doc");

                    py::class_<
                        tni::RungeKuttaVariableStepSizeSettingsScalarTolerances<
                            TIME_TYPE>,
                        std::shared_ptr<
                            tni::
                                RungeKuttaVariableStepSizeSettingsScalarTolerances<
                                    TIME_TYPE>>,
                        tni::RungeKuttaVariableStepSizeBaseSettings<TIME_TYPE>>(
                        m, "RungeKuttaVariableStepSizeSettingsScalarTolerances",
                        R"doc(`IntegratorSettings`-derived class to define settings for Runge Kutta integrators with scalar tolerances.

	This
	class is actually derived from an intermediate class (`RungeKuttaVariableStepSizeBaseSettings`, not documented),
	which is derived directly from `IntegratorSettings`.

)doc");

                    py::class_<
                        tni::BulirschStoerIntegratorSettings<TIME_TYPE>,
                        std::shared_ptr<
                            tni::BulirschStoerIntegratorSettings<TIME_TYPE>>,
                        tni::IntegratorSettings<TIME_TYPE>>(
                        m, "BulirschStoerIntegratorSettings",
                        R"doc(`IntegratorSettings`-derived class to define settings for Bulirsch-Stoer integrator settings.

)doc");


                    py::class_<
                        tni::AdamsBashforthMoultonSettings<TIME_TYPE>,
                        std::shared_ptr<
                            tni::AdamsBashforthMoultonSettings<TIME_TYPE>>,
                        tni::IntegratorSettings<TIME_TYPE>>(
                        m, "AdamsBashforthMoultonSettings",
                        R"doc(`IntegratorSettings`-derived class to define settings for Adams-Bashforth-Moulton integrator settings.

)doc");

                    py::class_<tni::IntegratorStepSizeControlSettings,
                               std::shared_ptr<
                                   tni::IntegratorStepSizeControlSettings>>(
                        m, "IntegratorStepSizeControlSettings",
                        get_docstring("IntegratorStepSizeControlSettings")
                            .c_str())
                        .def_readwrite("safety_factor",
                                       &tni::IntegratorStepSizeControlSettings::
                                           safetyFactorForNextStepSize_)
                        .def_readwrite(
                            "minimum_step_decrease",
                            &tni::IntegratorStepSizeControlSettings::
                                minimumFactorDecreaseForNextStepSize_)
                        .def_readwrite(
                            "maximum_step_decrease",
                            &tni::IntegratorStepSizeControlSettings::
                                maximumFactorDecreaseForNextStepSize_);

                    py::class_<tni::IntegratorStepSizeValidationSettings,
                               std::shared_ptr<
                                   tni::IntegratorStepSizeValidationSettings>>(
                        m, "IntegratorStepSizeValidationSettings",
                        get_docstring("IntegratorStepSizeValidationSettings")
                            .c_str())
                        .def_readwrite(
                            "minimum_step",
                            &tni::IntegratorStepSizeValidationSettings::
                                minimumStep_)
                        .def_readwrite(
                            "maximum_step",
                            &tni::IntegratorStepSizeValidationSettings::
                                maximumStep_)
                        .def_readwrite(
                            "minimum_step_handling",
                            &tni::IntegratorStepSizeValidationSettings::
                                minimumIntegrationTimeStepHandling_);


                    // FACTORY FUNCTIONS
                    m.def(
                        "print_butcher_tableau", &tni::printButcherTableau,
                        R"doc(Print the Butcher tableau of a given coefficient set.

	:param coefficient_set:
		Coefficient set of which the Butcher tableau will be printed.
)doc");

                    m.def("step_size_validation",
                          &tni::stepSizeValidationSettings,
                          py::arg("minimum_step"), py::arg("maximum_step"),
                          py::arg("minimum_step_size_handling") =
                              tni::throw_exception_below_minimum,
                          py::arg("accept_infinity_step") = false,
                          py::arg("accept_nan_step") = false,
                          get_docstring("step_size_validation").c_str());

                    m.def("step_size_control_elementwise_scalar_tolerance",
                          &tni::perElementIntegratorStepSizeControlSettings<
                              double>,
                          py::arg("relative_error_tolerance"),
                          py::arg("absolute_error_tolerance"),
                          py::arg("safety_factor") = 0.8,
                          py::arg("minimum_factor_increase") = 0.1,
                          py::arg("maximum_factor_increase") = 4.0,
                          get_docstring(
                              "step_size_control_elementwise_scalar_tolerance")
                              .c_str());

                    m.def("step_size_control_elementwise_matrix_tolerance",
                          &tni::perElementIntegratorStepSizeControlSettings<
                              Eigen::MatrixXd>,
                          py::arg("relative_error_tolerance"),
                          py::arg("absolute_error_tolerance"),
                          py::arg("safety_factor") = 0.8,
                          py::arg("minimum_factor_increase") = 0.1,
                          py::arg("maximum_factor_increase") = 4.0,
                          get_docstring(
                              "step_size_control_elementwise_matrix_tolerance")
                              .c_str());

                    m.def(
                        "step_size_control_blockwise_scalar_tolerance",
                        &tni::perBlockIntegratorStepSizeControlSettings<double>,
                        py::arg("block_indices"),
                        py::arg("relative_error_tolerance"),
                        py::arg("absolute_error_tolerance"),
                        py::arg("safety_factor") = 0.8,
                        py::arg("minimum_factor_increase") = 0.1,
                        py::arg("maximum_factor_increase") = 4.0,
                        get_docstring(
                            "step_size_control_blockwise_scalar_tolerance")
                            .c_str());

                    m.def("step_size_control_blockwise_matrix_tolerance",
                          &tni::perBlockIntegratorStepSizeControlSettings<
                              Eigen::MatrixXd>,
                          py::arg("block_indices"),
                          py::arg("relative_error_tolerance"),
                          py::arg("absolute_error_tolerance"),
                          py::arg("safety_factor") = 0.8,
                          py::arg("minimum_factor_increase") = 0.1,
                          py::arg("maximum_factor_increase") = 4.0,
                          get_docstring(
                              "step_size_control_blockwise_matrix_tolerance")
                              .c_str());

                    m.def(
                        "standard_cartesian_state_element_blocks",
                        &tni::getStandardCartesianStatesElementsToCheck,
                        py::arg("number_of_rows"), py::arg("number_of_columns"),
                        get_docstring("standard_cartesian_state_element_blocks")
                            .c_str());

                    m.def("standard_rotational_state_element_blocks",
                          &tni::getStandardRotationalStatesElementsToCheck,
                          py::arg("number_of_rows"),
                          py::arg("number_of_columns"));

                    m.def(
                        "step_size_control_custom_blockwise_scalar_tolerance",
                        &tni::
                            perBlockFromFunctionIntegratorStepSizeControlSettings<
                                double>,
                        py::arg("block_indices_function"),
                        py::arg("relative_error_tolerance"),
                        py::arg("absolute_error_tolerance"),
                        py::arg("safety_factor") = 0.8,
                        py::arg("minimum_factor_increase") = 0.1,
                        py::arg("maximum_factor_increase") = 4.0,
                        get_docstring("step_size_control_custom_blockwise_"
                                      "scalar_tolerance")
                            .c_str());

                    m.def(
                        "step_size_control_custom_blockwise_matrix_tolerance",
                        &tni::
                            perBlockFromFunctionIntegratorStepSizeControlSettings<
                                Eigen::MatrixXd>,
                        py::arg("block_indices_function"),
                        py::arg("relative_error_tolerance"),
                        py::arg("absolute_error_tolerance"),
                        py::arg("safety_factor") = 0.8,
                        py::arg("minimum_factor_increase") = 0.1,
                        py::arg("maximum_factor_increase") = 4.0,
                        get_docstring("step_size_control_custom_blockwise_"
                                      "matrix_tolerance")
                            .c_str());

                    m.def("runge_kutta_fixed_step",
                          &tni::rungeKuttaFixedStepSettings<TIME_TYPE>,
                          py::arg("time_step"), py::arg("coefficient_set"),
                          py::arg("order_to_use") =
                              tni::RungeKuttaCoefficients::
                                  OrderEstimateToIntegrate::lower,
                          py::arg("assess_termination_on_minor_steps") = false,
                          get_docstring("runge_kutta_fixed_step").c_str());

                    m.def("runge_kutta_variable_step",
                          &tni::multiStageVariableStepSizeSettings<TIME_TYPE>,
                          py::arg("initial_time_step"),
                          py::arg("coefficient_set"),
                          py::arg("step_size_control_settings"),
                          py::arg("step_size_validation_settings"),
                          py::arg("assess_termination_on_minor_steps") = false,
                          get_docstring("runge_kutta_variable_step").c_str());

                    m.def(
                        "bulirsch_stoer_variable_step",
                        &tni::bulirschStoerVariableStepIntegratorSettings<
                            TIME_TYPE>,
                        py::arg("initial_time_step"),
                        py::arg("extrapolation_sequence"),
                        py::arg("maximum_number_of_steps"),
                        py::arg("step_size_control_settings"),
                        py::arg("step_size_validation_settings"),
                        py::arg("assess_termination_on_minor_steps") = false,
                        get_docstring("bulirsch_stoer_variable_step").c_str());

                    m.def("bulirsch_stoer_fixed_step",
                          &tni::bulirschStoerFixedStepIntegratorSettings<
                              TIME_TYPE>,
                          py::arg("time_step"),
                          py::arg("extrapolation_sequence"),
                          py::arg("maximum_number_of_steps"),
                          py::arg("assess_termination_on_minor_steps") = false,
                          get_docstring("bulirsch_stoer_fixed_step").c_str());

                    m.def(
                        "adams_bashforth_moulton",
                        &tni::adamsBashforthMoultonSettings<TIME_TYPE>,
                        py::arg("initial_time_step"),
                        py::arg("minimum_step_size"),
                        py::arg("maximum_step_size"),
                        py::arg("relative_error_tolerance") = 1.0E-12,
                        py::arg("absolute_error_tolerance") = 1.0E-12,
                        py::arg("minimum_order") = 6,
                        py::arg("maximum_order") = 11,
                        py::arg("assess_termination_on_minor_steps") = false,
                        py::arg("bandwidth") = 200.0,
                        R"doc(Creates the settings for the Adams-Bashforth-Moulton integrator.

	Factory function to create settings for the Adams-Bashforth-Moulton integrator.
	For this integrator, the step size is varied based on the tolerances and safety factor provided.
	The tolerance is composed of an absolute and a relative part.


	:param initial_time_step:
		Initial time step to be used.
	:param minimum_step_size:
		Minimum time step to be used during the integration.
	:param maximum_step_size:
		Maximum time step to be used during the integration.
	:param relative_error_tolerance:
		Relative tolerance to adjust the time step.
	:param absolute_error_tolerance:
		Relative tolerance to adjust the time step.
	:param minimum_order:
		Minimum order of the integrator.
	:param maximum_order:
		Maximum order of the integrator.
	:param save_frequency:
		Frequency at which to save the numerical integrated states (expressed per unit integration time step, with n = saveFrequency, so n = 2 means that the state is saved every two integration steps).
	:param assess_termination_on_minor_steps:
		Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
	:param bandwidth:
		Maximum error factor for doubling the step size.
	:return:
		AdamsBashforthMoultonSettings object.
)doc");

                    m.def("adams_bashforth_moulton_fixed_order",
                          &tni::adamsBashforthMoultonSettingsFixedOrder<
                              TIME_TYPE>,
                          py::arg("initial_time_step"),
                          py::arg("minimum_step_size"),
                          py::arg("maximum_step_size"),
                          py::arg("relative_error_tolerance") = 1.0E-12,
                          py::arg("absolute_error_tolerance") = 1.0E-12,
                          py::arg("order") = 6,
                          py::arg("assess_termination_on_minor_steps") = false,
                          py::arg("bandwidth") = 200.0,
                          get_docstring("adams_bashforth_moulton_fixed_order")
                              .c_str());

                    m.def(
                        "adams_bashforth_moulton_fixed_step",
                        &tni::adamsBashforthMoultonSettingsFixedStep<TIME_TYPE>,
                        py::arg("time_step"),
                        py::arg("relative_error_tolerance") = 1.0E-12,
                        py::arg("absolute_error_tolerance") = 1.0E-12,
                        py::arg("minimum_order") = 6,
                        py::arg("maximum_order") = 11,
                        py::arg("assess_termination_on_minor_steps") = false,
                        py::arg("bandwidth") = 200.0,
                        get_docstring("adams_bashforth_moulton_fixed_step")
                            .c_str());

                    m.def(
                        "adams_bashforth_moulton_fixed_step_fixed_order",
                        &tni::adamsBashforthMoultonSettingsFixedStepFixedOrder<
                            TIME_TYPE>,
                        py::arg("time_step"), py::arg("order") = 6,
                        py::arg("assess_termination_on_minor_steps") = false,
                        get_docstring(
                            "adams_bashforth_moulton_fixed_step_fixed_order")
                            .c_str());

                    /*!
                     * DEPRECATED
                     * -------------------------------------------------------------------------
                     *
                     */


                    m.def(
                        "runge_kutta_variable_step_size_vector_tolerances",
                        &tni::rungeKuttaVariableStepSettingsVectorTolerances<
                            TIME_TYPE>,
                        py::arg("initial_time_step"),
                        py::arg("coefficient_set"),
                        py::arg("minimum_step_size"),
                        py::arg("maximum_step_size"),
                        py::arg("relative_error_tolerance"),
                        py::arg("absolute_error_tolerance"),
                        py::arg("assess_termination_on_minor_steps") = false,
                        py::arg("safety_factor") = 0.8,
                        py::arg("maximum_factor_increase") = 4.0,
                        py::arg("minimum_factor_increase") = 0.1,
                        py::arg("throw_exception_if_minimum_step_exceeded") =
                            true,
                        R"doc(Creates the settings for the Runge-Kutta variable step size integrator with vector tolerances.

	Factory function to create settings for the Runge-Kutta variable step size integrator with vector tolerances.
	For this integrator, the step size is varied based on the tolerances and safety factor provided.
	The tolerance is composed of an absolute and a relative part.
	Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).


	:param initial_time_step:
		Initial time step to be used.
	:param coefficient_set:
		Coefficient set (Butcher's tableau) to be used in the integration.
	:param minimum_step_size:
		Minimum time step to be used during the integration.
	:param maximum_step_size:
		Maximum time step to be used during the integration.
	:param relative_error_tolerance:
		Relative vector tolerance to adjust the time step.
	:param absolute_error_tolerance:
		Absolute vector tolerance to adjust the time step.
	:param save_frequency:
		Frequency at which to save the numerical integrated states (expressed per unit integration time step,
		with n = saveFrequency, so n = 2 means that the state is saved every two integration steps).

	:param assess_termination_on_minor_steps:
		Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
		integrator (true) or only at the end of each integration step (false).

	:param safety_factor:
		Safety factor used in the step size control.
	:param maximum_factor_increase:
		Maximum increase between consecutive time steps, expressed as the factor between new and old step size.

	:param minimum_factor_increase:
		Minimum increase between consecutive time steps, expressed as the factor between new and old step size.

	:param throw_exception_if_minimum_step_exceeded:
		If set to false, the variable step integrator will use the minimum step size specified when the algorithm
		computes the optimum one to be lower, instead of throwing an exception.

	:return:
		RungeKuttaVariableStepSizeSettingsVectorTolerances object.
)doc");

                    m.def(
                        "runge_kutta_variable_step_size",
                        &tni::rungeKuttaVariableStepSettingsScalarTolerances<
                            TIME_TYPE>,
                        py::arg("initial_time_step"),
                        py::arg("coefficient_set"),
                        py::arg("minimum_step_size"),
                        py::arg("maximum_step_size"),
                        py::arg("relative_error_tolerance"),
                        py::arg("absolute_error_tolerance"),
                        py::arg("assess_termination_on_minor_steps") = false,
                        py::arg("safety_factor") = 0.8,
                        py::arg("maximum_factor_increase") = 4.0,
                        py::arg("minimum_factor_increase") = 0.1,
                        py::arg("throw_exception_if_minimum_step_exceeded") =
                            true,
                        R"doc(Creates the settings for the Runge-Kutta variable step size integrator with scalar tolerances.

	Factory function to create settings for the Runge-Kutta variable step size integrator with scalar tolerances.
	For this integrator, the step size is varied based on the tolerances and safety factor provided.
	The tolerance is composed of an absolute and a relative part.
	Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).


	:param initial_time_step:
		Initial time step to be used.
	:param coefficient_set:
		Coefficient set (Butcher's tableau) to be used in the integration.
	:param minimum_step_size:
		Minimum time step to be used during the integration.
	:param maximum_step_size:
		Maximum time step to be used during the integration.
	:param relative_error_tolerance:
		Relative vector tolerance to adjust the time step.
	:param absolute_error_tolerance:
		Absolute vector tolerance to adjust the time step.
	:param save_frequency:
		Frequency at which to save the numerical integrated states (expressed per unit integration time step,
		with n = saveFrequency, so n = 2 means that the state is saved every two integration steps).

	:param assess_termination_on_minor_steps:
		Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
		integrator (true) or only at the end of each integration step (false).

	:param safety_factor:
		Safety factor used in the step size control.
	:param maximum_factor_increase:
		Maximum increase between consecutive time steps, expressed as the factor between new and old step size.

	:param minimum_factor_increase:
		Minimum increase between consecutive time steps, expressed as the factor between new and old step size.

	:param throw_exception_if_minimum_step_exceeded:
		If set to false, the variable step integrator will use the minimum step size specified when the algorithm
		computes the optimum one to be lower, instead of throwing an exception.

	:return:
		RungeKuttaVariableStepSettingsScalarTolerances object.
)doc");


                    /*!
                     * DEPRECATED UNDOCUMENTED
                     * -------------------------------------------------------------------------
                     *
                     */

                    m.def("euler", &tni::eulerSettingsDeprecated<TIME_TYPE>,
                          py::arg("initial_time"), py::arg("initial_time_step"),
                          py::arg("assess_termination_on_minor_steps") = false);

                    m.def("euler", &tni::eulerSettings<TIME_TYPE>,
                          py::arg("initial_time_step"),
                          py::arg("assess_termination_on_minor_steps") = false);

                    m.def("runge_kutta_4",
                          &tni::rungeKutta4SettingsDeprecated<TIME_TYPE>,
                          py::arg("initial_time"), py::arg("initial_time_step"),
                          py::arg("assess_termination_on_minor_steps") = false);

                    m.def("runge_kutta_4", &tni::rungeKutta4Settings<TIME_TYPE>,
                          py::arg("initial_time_step"),
                          py::arg("assess_termination_on_minor_steps") = false);

                    m.def(
                        "runge_kutta_fixed_step_size",
                        &tni::rungeKuttaFixedStepSettingsDeprecated<TIME_TYPE>,
                        py::arg("initial_time"), py::arg("initial_time_step"),
                        py::arg("coefficient_set"),
                        py::arg("order_to_use") = tni::RungeKuttaCoefficients::
                            OrderEstimateToIntegrate::lower,
                        py::arg("assess_termination_on_minor_steps") = false);

                    m.def("runge_kutta_fixed_step_size",
                          &tni::rungeKuttaFixedStepSettings<TIME_TYPE>,
                          py::arg("initial_time_step"),
                          py::arg("coefficient_set"),
                          py::arg("order_to_use") =
                              tni::RungeKuttaCoefficients::
                                  OrderEstimateToIntegrate::lower,
                          py::arg("assess_termination_on_minor_steps") = false);

                    m.def(
                        "runge_kutta_variable_step_size",
                        &tni::
                            rungeKuttaVariableStepSettingsScalarTolerancesDeprecated<
                                TIME_TYPE>,
                        py::arg("initial_time"), py::arg("initial_time_step"),
                        py::arg("coefficient_set"),
                        py::arg("minimum_step_size"),
                        py::arg("maximum_step_size"),
                        py::arg("relative_error_tolerance"),
                        py::arg("absolute_error_tolerance"),
                        py::arg("assess_termination_on_minor_steps") = false,
                        py::arg("safety_factor") = 0.8,
                        py::arg("maximum_factor_increase") = 4.0,
                        py::arg("minimum_factor_increase") = 0.1,
                        py::arg("throw_exception_if_minimum_step_exceeded") =
                            true);


                    m.def(
                        "runge_kutta_variable_step_size_vector_tolerances",
                        &tni::
                            rungeKuttaVariableStepSettingsVectorTolerancesDeprecated<
                                TIME_TYPE>,
                        py::arg("initial_time"), py::arg("initial_time_step"),
                        py::arg("coefficient_set"),
                        py::arg("minimum_step_size"),
                        py::arg("maximum_step_size"),
                        py::arg("relative_error_tolerance"),
                        py::arg("absolute_error_tolerance"),
                        py::arg("assess_termination_on_minor_steps") = false,
                        py::arg("safety_factor") = 0.8,
                        py::arg("maximum_factor_increase") = 4.0,
                        py::arg("minimum_factor_increase") = 0.1,
                        py::arg("throw_exception_if_minimum_step_exceeded") =
                            true);


                    m.def("bulirsch_stoer",
                          &tni::bulirschStoerIntegratorSettingsDeprecated<
                              TIME_TYPE>,
                          py::arg("initial_time"), py::arg("initial_time_step"),
                          py::arg("extrapolation_sequence"),
                          py::arg("maximum_number_of_steps"),
                          py::arg("minimum_step_size"),
                          py::arg("maximum_step_size"),
                          py::arg("relative_error_tolerance") = 1.0E-12,
                          py::arg("absolute_error_tolerance") = 1.0E-12,
                          py::arg("assess_termination_on_minor_steps") = false,
                          py::arg("safety_factor") = 0.7,
                          py::arg("maximum_factor_increase") = 10.0,
                          py::arg("minimum_factor_increase") = 0.1);

                    m.def(
                        "bulirsch_stoer",
                        &tni::bulirschStoerIntegratorSettingsDeprecatedNew<
                            TIME_TYPE>,
                        py::arg("initial_time_step"),
                        py::arg("extrapolation_sequence"),
                        py::arg("maximum_number_of_steps"),
                        py::arg("minimum_step_size"),
                        py::arg("maximum_step_size"),
                        py::arg("relative_error_tolerance") = 1.0E-12,
                        py::arg("absolute_error_tolerance") = 1.0E-12,
                        py::arg("assess_termination_on_minor_steps") = false,
                        py::arg("safety_factor") = 0.7,
                        py::arg("maximum_factor_increase") = 10.0,
                        py::arg("minimum_factor_increase") = 0.1,
                        R"doc(Creates the settings for the Bulirsch-Stoer integrator.

	Factory function to create settings for the Bulirsch-Stoer integrator.
	For this integrator, the step size is varied based on the tolerances and safety factor provided.
	The tolerance is composed of an absolute and a relative part.
	Different extrapolation sequences can be used (see the `ExtrapolationMethodStepSequences` enum).


	:param initial_time_step:
		Initial time step to be used.
	:param extrapolation_sequence:
		Extrapolation sequence to be used in the integration.
	:param maximum_number_of_steps:
		Number of entries in the sequence (e.g., number of integrations used for a single extrapolation).
	:param minimum_step_size:
		Minimum time step to be used during the integration.
	:param maximum_step_size:
		Maximum time step to be used during the integration.
	:param relative_error_tolerance:
		Relative tolerance to adjust the time step.
	:param absolute_error_tolerance:
		Relative tolerance to adjust the time step.
	:param save_frequency:
		Frequency at which to save the numerical integrated states (expressed per unit integration time step, with n = saveFrequency, so n = 2 means that the state is saved every two integration steps).
	:param assess_termination_on_minor_steps:
		Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
	:param safety_factor:
		Safety factor used in the step size control.
	:param maximum_factor_increase:
		Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
	:param minimum_factor_increase:
		Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
	:return:
		BulirschStoerIntegratorSettings object.
)doc");

                    m.def("adams_bashforth_moulton",
                          &tni::adamsBashforthMoultonSettingsDeprecated<
                              TIME_TYPE>,
                          py::arg("initial_time"), py::arg("initial_time_step"),
                          py::arg("minimum_step_size"),
                          py::arg("maximum_step_size"),
                          py::arg("relative_error_tolerance") = 1.0E-12,
                          py::arg("absolute_error_tolerance") = 1.0E-12,
                          py::arg("minimum_order") = 6,
                          py::arg("maximum_order") = 11,
                          py::arg("assess_termination_on_minor_steps") = false,
                          py::arg("bandwidth") = 200.0);
                }


            }  // namespace integrator
        }  // namespace propagation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
