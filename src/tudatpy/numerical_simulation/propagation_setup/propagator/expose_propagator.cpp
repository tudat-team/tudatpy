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
#include <tudat/astro/propagators/getZeroProperModeRotationalInitialState.h>
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
            namespace propagator {

                PYBIND11_MODULE(expose_propagator, m) {
                    py::class_<tp::PropagationPrintSettings,
                               std::shared_ptr<tp::PropagationPrintSettings>>(
                        m, "PropagationPrintSettings",
                        R"doc(Class to save settings on what is to be written to the console during the propagation of a single arc.

)doc")
                        .def_property(
                            "print_dependent_variable_indices",
                            &tp::PropagationPrintSettings::
                                getPrintDependentVariableData,
                            &tp::PropagationPrintSettings::
                                setPrintDependentVariableData,
                            R"doc(Boolean defining whether the meaning and indices of the
entries of the dependent variable data are to be printed to
the console (before the propagation).

.. note:: The same information can be retrieved from the
          :py:attr:`SingleArcSimulationResults.dependent_variable_ids`
          attribute.

	)doc")
                        .def_property(
                            "print_state_indices",
                            &tp::PropagationPrintSettings::
                                getPrintPropagatedStateData,
                            &tp::PropagationPrintSettings::
                                setPrintPropagatedStateData,
                            R"doc(Boolean defining whether the meaning and indices of the
entries of the state vector are to be printed to
the console (before the propagation).

.. note:: The same information can be retrieved from the
          :py:attr:`SingleArcSimulationResults.propagated_state_ids`
          attribute.

	)doc")
                        .def_property(
                            "print_processed_state_indices",
                            &tp::PropagationPrintSettings::
                                getPrintProcessedStateData,
                            &tp::PropagationPrintSettings::
                                setPrintProcessedStateData,
                            R"doc(Boolean defining whether the meaning and indices of the
entries of the processed state vector are to be printed to
the console (after the propagation). The distinction between the
propagated and processed (or conventuional) state representation is described in
detail `here <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/propagator_types.html>`_.
Summarizing: the processed state is the 'typical' formulation of the state (for translational dynamics: Cartesian states).

.. note:: The same information can be retrieved from the
          :py:attr:`SingleArcSimulationResults.processed_state_ids`
          attribute.


	)doc")
                        .def_property(
                            "print_number_of_function_evaluations",
                            &tp::PropagationPrintSettings::
                                getPrintNumberOfFunctionEvaluations,
                            &tp::PropagationPrintSettings::
                                setPrintNumberOfFunctionEvaluations,
                            R"doc(Boolean defining whether the number of function evaluations that
were performed is to be printed to the console (after propagation).

	)doc")
                        .def_property(
                            "print_propagation_clock_time",
                            &tp::PropagationPrintSettings::
                                getPrintPropagationTime,
                            &tp::PropagationPrintSettings::
                                setPrintPropagationTime,
                            R"doc(Boolean defining whether the total clock time taken for the propagation
is to be printed to the console (after propagation).

	)doc")
                        .def_property(
                            "print_termination_reason",
                            &tp::PropagationPrintSettings::
                                getPrintTerminationReason,
                            &tp::PropagationPrintSettings::
                                setPrintTerminationReason,
                            R"doc(Boolean defining whether the reason for propagation termination
is to be printed to the console (after propagation).

	)doc")
                        .def_property(
                            "print_initial_and_final_conditions",
                            &tp::PropagationPrintSettings::
                                getPrintInitialAndFinalConditions,
                            &tp::PropagationPrintSettings::
                                setPrintInitialAndFinalConditions,
                            R"doc(Boolean defining whether the initial and final conditions (state and time)
are to be printed to the console (beforee and after propagation, respectively).

	)doc")
                        .def_property(
                            "results_print_frequency_in_seconds",
                            &tp::PropagationPrintSettings::
                                getResultsPrintFrequencyInSeconds,
                            &tp::PropagationPrintSettings::
                                setResultsPrintFrequencyInSeconds,
                            R"doc(Variable indicating how often (in seconds of simulation time)
the current state and time are to be printed to the console (by default, set to NaN - they are never printed).
In case this setting is active (e.g. not NaN), and the ``results_print_frequency_in_steps`` setting is active,
the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.

	)doc")
                        .def_property(
                            "results_print_frequency_in_steps",
                            &tp::PropagationPrintSettings::
                                getResultsPrintFrequencyInSteps,
                            &tp::PropagationPrintSettings::
                                setResultsPrintFrequencyInSteps,
                            R"doc(Variable indicating how often (in number of full integration steps)
the current state and time are to be printed to the console (by default, set to 0 - they are never printed).
In case this setting is active (e.g. not 0), and the ``results_print_frequency_in_seconds`` setting is active,
the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.

	)doc")
                        .def_property(
                            "print_dependent_variables_during_propagation",
                            &tp::PropagationPrintSettings::
                                getPrintDependentVariableDuringPropagation,
                            &tp::PropagationPrintSettings::
                                setPrintDependentVariableDuringPropagation,
                            R"doc(Boolean defining whether the dependent variables are to be printed during the propagation along with the state,
at steps/epochs define by the ``results_print_frequency_in_seconds`` and/or ``results_print_frequency_in_steps`` inputs.

	)doc")
                        .def(
                            "enable_all_boolean_printing",
                            py::overload_cast<>(&tp::PropagationPrintSettings::
                                                    enableAllPrinting),
                            R"doc(Function enabling all True/False printing (e.g. sets all boolean attributes to True)

)doc")
                        .def(
                            "enable_all_printing",
                            py::overload_cast<const double, const int>(
                                &tp::PropagationPrintSettings::
                                    enableAllPrinting),
                            py::arg("results_print_frequency_in_seconds"),
                            py::arg("results_print_frequency_in_steps"),
                            R"doc(Function enabling all True/False printing (e.g. sets all boolean attributes to True), and setting the non-boolean
attributes to values defined here.



	:param results_print_frequency_in_seconds:
		See ``results_print_frequency_in_seconds`` class attribute
	:param results_print_frequency_in_steps:
		See ``results_print_frequency_in_steps`` class attribute
)doc")
                        .def(
                            "disable_all_printing",
                            &tp::PropagationPrintSettings::disableAllPrinting,
                            R"doc(Function enabling all printing (e.g. sets all boolean attributes to False, and disables all other output as well)

)doc");

                    py::class_<
                        tp::PropagatorProcessingSettings,
                        std::shared_ptr<tp::PropagatorProcessingSettings>>(
                        m, "PropagatorProcessingSettings",
                        R"doc(Base class to define settings on how the numerical results are to be used

	Base class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
	Instances of this class are typically not created by the user. Settings objects for derived class of single-, multi- and hybrid arc propagation are
	instantiated through the factory functions to define propagator settings (such as :func:`~translational` or :func:`~multi_arc`) in this module


)doc")
                        .def_property(
                            "set_integrated_result",
                            &tp::PropagatorProcessingSettings::
                                getSetIntegratedResult,
                            &tp::PropagatorProcessingSettings::
                                setIntegratedResult,
                            R"doc(Boolean defining whether the propagation results are to
be used to update the environment. If this variable is set
to False, the numerical propagation results can be
retrieved from this object (provided the
:py:attr:`~clear_numerical_solution` is set to False),
but the (for instance) Ephemeris of the propagated body
is not updated with the propagation results. If this
variable is set to True, the properties of the propagated
:class:`~tudatpy.numerical_simulation.environment.Body`
object will be updated as per the numerical results
(see `here <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/console_output.html#automatic-processing>`_ for details).

	)doc")
                        .def_property(
                            "clear_numerical_solution",
                            &tp::PropagatorProcessingSettings::
                                getClearNumericalSolutions,
                            &tp::PropagatorProcessingSettings::
                                setClearNumericalSolutions,
                            R"doc(Boolean defining whether the propagation results should be
deleted after the propagation is terminated. If this is
done, the :py:attr:`~state_history`,
:py:attr:`~unprocessed_state_history` and
:py:attr:`~dependent_variable_history` will not be
available in the :py:class:`~tudatpy.numerical_simulation.propagator.SingleArcSimulationResults` class. Putting this setting to True (deleting the
results) is only sensible when the
:py:attr:`~set_integrated_result` is set to True. In that
case, the propagated states are *not* accessible directly
from this objects, but the results are used to update the
environment, *e.g.* update the ephemeris of the propagated
body with the numerical results.

	)doc")
                        .def_property(
                            "create_dependent_variable_interface",
                            &tp::PropagatorProcessingSettings::
                                getUpdateDependentVariableInterpolator,
                            &tp::PropagatorProcessingSettings::
                                setUpdateDependentVariableInterpolator,
                            get_docstring("PropagatorProcessingSettings.create_"
                                          "dependent_variable_interface")
                                .c_str());

                    py::class_<tp::SingleArcPropagatorProcessingSettings,
                               std::shared_ptr<
                                   tp::SingleArcPropagatorProcessingSettings>,
                               tp::PropagatorProcessingSettings>(
                        m, "SingleArcPropagatorProcessingSettings",
                        R"doc(Class to define settings on how the numerical results are to be used for single-arc propagations

	Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
	derived from :class:`PropagatorProcessingSettings`.
	Instances of this class are typically not created by the user. A settings object is
	instantiated through the factory functions to define single-arc propagator settings (such as :func:`~translational` or :func:`~rotational`) in this module

)doc")
                        .def_property_readonly(
                            "print_settings",
                            &tp::SingleArcPropagatorProcessingSettings::
                                getPrintSettings,
                            R"doc(Settings object defining which quantities should be printed to the console before, during and after the propagation. By default, this
object is instantiated to print nothing.

	)doc")
                        .def_property(
                            "results_save_frequency_in_steps",
                            &tp::SingleArcPropagatorProcessingSettings::
                                getResultsSaveFrequencyInSteps,
                            &tp::SingleArcPropagatorProcessingSettings::
                                setResultsSaveFrequencyInSteps,
                            R"doc(Variable indicating how often (in number of integrator steps)
the propagated time, state, dependent variables, etc. are to be saved to data structures containing the results
(by default, set to 1 - they are saved every time step). If this setting is set to 0, the data is never saved based on number of steps.
In case this setting is active (e.g. not 0), and the ``results_save_frequency_in_seconds`` setting is active,
the data is saved as soon as *one* of the two conditions (number of seconds, or number of steps) is met.

	)doc")
                        .def_property(
                            "results_save_frequency_in_seconds",
                            &tp::SingleArcPropagatorProcessingSettings::
                                getResultsSaveFrequencyInSeconds,
                            &tp::SingleArcPropagatorProcessingSettings::
                                setResultsSaveFrequencyInSeconds,
                            R"doc(Variable indicating how often (in seconds of simulation time)
the propagated time, state, dependent variables, etc. are to be saved to data structures containing the results
(by default, set to NaN - they are not saved based on propagation time; see below and ``results_save_frequency_in_steps`` attribute ).
In case this setting is active, and set to :math:`\Delta t`, the data are saved as soon as the current time step is :math:`\ge \Delta t` after the
last step at which data was saved.
In case this setting is active (e.g. not NaN), and the ``results_save_frequency_in_steps`` setting is active,
the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.

	)doc");

                    py::class_<tp::MultiArcPropagatorProcessingSettings,
                               std::shared_ptr<
                                   tp::MultiArcPropagatorProcessingSettings>,
                               tp::PropagatorProcessingSettings>(
                        m, "MultiArcPropagatorProcessingSettings",
                        R"doc(Class to define settings on how the numerical results are to be used for multi-arc propagations

	Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
	derived from :class:`PropagatorProcessingSettings`.
	Instances of this class are typically not created by the user. A settings object is
	instantiated through the factory function :func:`~multi_arc` to define multi-arc propagator settings.
	This object contains a list of :class:`SingleArcPropagatorProcessingSettings` objects, containing the processing settings for each constituent arc.

)doc")
                        .def(
                            "set_print_settings_for_all_arcs",
                            &tp::MultiArcPropagatorProcessingSettings::
                                resetAndApplyConsistentSingleArcPrintSettings,
                            py::arg("single_arc_print_settings"),
                            R"doc(Function that sets the same print settings for each arc in the multi-arc propagation.


	:param single_arc_print_settings:
		Propagation print settings that are applied to each constituent single-arc settings, overriding any existing settings.
)doc")
                        .def_property(
                            "print_output_on_first_arc_only",
                            &tp::MultiArcPropagatorProcessingSettings::
                                getPrintFirstArcOnly,
                            &tp::MultiArcPropagatorProcessingSettings::
                                resetPrintFirstArcOnly,
                            R"doc(Variable defining whether the ``set_print_settings_for_all_arcs`` function has been used to define identical print settings for each arc.

	)doc")
                        .def_property_readonly(
                            "identical_settings_per_arc",
                            &tp::MultiArcPropagatorProcessingSettings::
                                useIdenticalSettings,
                            get_docstring("MultiArcPropagatorProcessingSettings"
                                          ".identical_settings_per_arc")
                                .c_str())
                        .def_property_readonly(
                            "single_arc_settings",
                            &tp::MultiArcPropagatorProcessingSettings::
                                getSingleArcSettings,
                            R"doc(List containing the processing settings for each constituent arc

	)doc");


                    py::class_<tp::HybridArcPropagatorProcessingSettings,
                               std::shared_ptr<
                                   tp::HybridArcPropagatorProcessingSettings>,
                               tp::PropagatorProcessingSettings>(
                        m, "HybridArcPropagatorProcessingSettings",
                        R"doc(Class to define settings on how the numerical results are to be used for hybrid-arc propagations

	Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
	derived from :class:`PropagatorProcessingSettings`.
	Instances of this class are typically not created by the user. A settings object is
	instantiated through the factory function :func:`~hybrid_arc` to define hybrid-arc propagator settings.
	This object contains a :class:`SingleArcPropagatorProcessingSettings` object and a :class:`MultuArcPropagatorProcessingSettings` ,
	containing the processing settings for the constituents of the hybrid-arc propagatioon.

)doc")
                        .def_property(
                            "set_integrated_result",
                            &tp::HybridArcPropagatorProcessingSettings::
                                getSetIntegratedResult,
                            &tp::HybridArcPropagatorProcessingSettings::
                                setIntegratedResult)
                        .def_property(
                            "clear_numerical_solution",
                            &tp::HybridArcPropagatorProcessingSettings::
                                getClearNumericalSolutions,
                            &tp::HybridArcPropagatorProcessingSettings::
                                setClearNumericalSolutions)
                        .def(
                            "set_print_settings_for_all_arcs",
                            &tp::HybridArcPropagatorProcessingSettings::
                                resetAndApplyConsistentPrintSettings,
                            py::arg("print_settings"),
                            R"doc(Function that sets the same print settings for each arc in the multi-arc propagation.


	:param single_arc_print_settings:
		Propagation print settings that are applied to each constituent single-arc settings, overriding any existing settings.
)doc")
                        .def_property_readonly(
                            "single_arc_settings",
                            &tp::HybridArcPropagatorProcessingSettings::
                                getSingleArcSettings,
                            R"doc(Processing settings for the single-arc component of the hybrid-arc propagation.

	)doc")
                        .def_property_readonly(
                            "multi_arc_settings",
                            &tp::HybridArcPropagatorProcessingSettings::
                                getMultiArcSettings,
                            R"doc(Processing settings for the single-arc component of the multi-arc propagation.

	)doc");

                    ///////////////////////////////////////////////////////////////////////////////////////

                    // ENUMS
                    py::enum_<tp::TranslationalPropagatorType>(
                        m, "TranslationalPropagatorType",
                        R"doc(Enumeration of available translational propagator types.


	:member undefined_translational_propagator:
	:member cowell:
	:member encke:
	:member gauss_keplerian:
	:member gauss_modified_equinoctial:
	:member unified_state_model_quaternions:
	:member unified_state_model_modified_rodrigues_parameters:
	:member unified_state_model_exponential_map:
)doc")
                        .value(
                            "undefined_translational_propagator",
                            tp::TranslationalPropagatorType::
                                undefined_translational_propagator,
                            get_docstring("TranslationalPropagatorType."
                                          "undefined_translational_propagator")
                                .c_str())
                        .value(
                            "cowell", tp::TranslationalPropagatorType::cowell,
                            get_docstring("TranslationalPropagatorType.cowell")
                                .c_str())
                        .value(
                            "encke", tp::TranslationalPropagatorType::encke,
                            get_docstring("TranslationalPropagatorType.encke")
                                .c_str())
                        .value(
                            "gauss_keplerian",
                            tp::TranslationalPropagatorType::gauss_keplerian,
                            get_docstring(
                                "TranslationalPropagatorType.gauss_keplerian")
                                .c_str())
                        .value("gauss_modified_equinoctial",
                               tp::TranslationalPropagatorType::
                                   gauss_modified_equinoctial,
                               get_docstring("TranslationalPropagatorType."
                                             "gauss_modified_equinoctial")
                                   .c_str())
                        .value("unified_state_model_quaternions",
                               tp::TranslationalPropagatorType::
                                   unified_state_model_quaternions,
                               get_docstring("TranslationalPropagatorType."
                                             "unified_state_model_quaternions")
                                   .c_str())
                        .value(
                            "unified_state_model_modified_rodrigues_parameters",
                            tp::TranslationalPropagatorType::
                                unified_state_model_modified_rodrigues_parameters,
                            get_docstring(
                                "TranslationalPropagatorType.unified_state_"
                                "model_modified_rodrigues_parameters")
                                .c_str())
                        .value(
                            "unified_state_model_exponential_map",
                            tp::unified_state_model_exponential_map,
                            get_docstring("TranslationalPropagatorType.unified_"
                                          "state_model_exponential_map")
                                .c_str())
                        .export_values();

                    py::enum_<tp::RotationalPropagatorType>(
                        m, "RotationalPropagatorType",
                        R"doc(Enumeration of available rotational propagator types.


	:member undefined_rotational_propagator:
	:member quaternions:
	:member modified_rodrigues_parameters:
	:member exponential_map:
)doc")
                        .value("undefined_rotational_propagator",
                               tp::RotationalPropagatorType::
                                   undefined_rotational_propagator,
                               get_docstring("RotationalPropagatorType."
                                             "undefined_rotational_propagator")
                                   .c_str())
                        .value("quaternions",
                               tp::RotationalPropagatorType::quaternions,
                               get_docstring(
                                   "RotationalPropagatorType.quaternions")
                                   .c_str())
                        .value("modified_rodrigues_parameters",
                               tp::RotationalPropagatorType::
                                   modified_rodrigues_parameters,
                               get_docstring("RotationalPropagatorType."
                                             "modified_rodrigues_parameters")
                                   .c_str())
                        .value("exponential_map",
                               tp::RotationalPropagatorType::exponential_map,
                               get_docstring(
                                   "RotationalPropagatorType.exponential_map")
                                   .c_str())
                        .export_values();

                    py::enum_<tp::PropagationTerminationTypes>(
                        m, "PropagationTerminationTypes",
                        R"doc(Enumeration of possible propagation termination types


	:member time_stopping_condition:
	:member cpu_time_stopping_condition:
	:member dependent_variable_stopping_condition:
	:member hybrid_stopping_condition:
	:member custom_stopping_condition:
)doc")
                        .value("time_stopping_condition_type",
                               tp::PropagationTerminationTypes::
                                   time_stopping_condition,
                               get_docstring("PropagationTerminationTypes.time_"
                                             "stopping_condition_type")
                                   .c_str())
                        .value("cpu_time_stopping_condition_type",
                               tp::PropagationTerminationTypes::
                                   cpu_time_stopping_condition,
                               get_docstring("PropagationTerminationTypes.cpu_"
                                             "time_stopping_condition_type")
                                   .c_str())
                        .value("dependent_variable_stopping_condition_type",
                               tp::PropagationTerminationTypes::
                                   dependent_variable_stopping_condition,
                               get_docstring(
                                   "PropagationTerminationTypes.dependent_"
                                   "variable_stopping_condition_type")
                                   .c_str())
                        .value("hybrid_stopping_condition_type",
                               tp::PropagationTerminationTypes::
                                   hybrid_stopping_condition,
                               get_docstring("PropagationTerminationTypes."
                                             "hybrid_stopping_condition_type")
                                   .c_str())
                        .value("custom_stopping_condition_type",
                               tp::PropagationTerminationTypes::
                                   custom_stopping_condition,
                               get_docstring("PropagationTerminationTypes."
                                             "custom_stopping_condition_type")
                                   .c_str())
                        .export_values();

                    py::enum_<tp::IntegratedStateType>(
                        m, "StateType",
                        R"doc(Enumeration of available integrated state types.


	:member hybrid_type:
	:member translational_type:
	:member rotational_type:
	:member body_mass_type:
	:member custom_type:
)doc")
                        .value("hybrid_type", tp::IntegratedStateType::hybrid,
                               get_docstring("StateType.hybrid_type").c_str())
                        .value("translational_type",
                               tp::IntegratedStateType::translational_state,
                               get_docstring("StateType.translational_type")
                                   .c_str())
                        .value(
                            "rotational_type",
                            tp::IntegratedStateType::rotational_state,
                            get_docstring("StateType.rotational_type").c_str())
                        .value("mass_type",
                               tp::IntegratedStateType::body_mass_state,
                               get_docstring("StateType.mass_type").c_str())
                        .value("custom_type",
                               tp::IntegratedStateType::custom_state,
                               get_docstring("StateType.custom_type").c_str())
                        .export_values();


                    //        py::enum_<tp::VariableType>(m, "VariableType")
                    //                .value("independent_variable",
                    //                tp::VariableType::independentVariable)
                    //                .value("cpu_time_variable",
                    //                tp::VariableType::cpuTimeVariable)
                    //                .value("state_variable",
                    //                tp::VariableType::stateVariable)
                    //                .value("dependent_variable",
                    //                tp::VariableType::dependentVariable)
                    //                .export_values();


                    // CLASSES
                    py::class_<tp::PropagatorSettings<double>,
                               std::shared_ptr<tp::PropagatorSettings<double>>>(
                        m, "PropagatorSettings",
                        R"doc(Functional base class to define settings for propagators.

	Base class to define settings for propagators. Derived classes are split into settings for single- and multi-arc dynamics.

)doc")
                        .def_property(
                            "initial_states",
                            &tp::PropagatorSettings<double>::getInitialStates,
                            &tp::PropagatorSettings<double>::resetInitialStates,
                            get_docstring("PropagatorSettings.initial_states")
                                .c_str());

                    py::class_<
                        tp::MultiArcPropagatorSettings<double, TIME_TYPE>,
                        std::shared_ptr<
                            tp::MultiArcPropagatorSettings<double, TIME_TYPE>>,
                        tp::PropagatorSettings<double>>(
                        m, "MultiArcPropagatorSettings",
                        R"doc(`PropagatorSettings`-derived class to define settings for multi-arc dynamics.

)doc")
                        .def_property_readonly(
                            "processing_settings",
                            &tp::MultiArcPropagatorSettings<
                                double, TIME_TYPE>::getOutputSettings,
                            get_docstring(
                                "MultiArcPropagatorSettings.print_settings")
                                .c_str());

                    py::class_<
                        tp::HybridArcPropagatorSettings<double, TIME_TYPE>,
                        std::shared_ptr<
                            tp::HybridArcPropagatorSettings<double, TIME_TYPE>>,
                        tp::PropagatorSettings<double>>(
                        m, "HybridArcPropagatorSettings",
                        R"doc(`PropagatorSettings`-derived class to define settings for hybrid-arc dynamics.

)doc")
                        .def_property_readonly(
                            "processing_settings",
                            &tp::HybridArcPropagatorSettings<
                                double, TIME_TYPE>::getOutputSettings,
                            get_docstring(
                                "HybridArcPropagatorSettings.print_settings")
                                .c_str());

                    py::class_<
                        tp::SingleArcPropagatorSettings<double, TIME_TYPE>,
                        std::shared_ptr<
                            tp::SingleArcPropagatorSettings<double, TIME_TYPE>>,
                        tp::PropagatorSettings<double>>(
                        m, "SingleArcPropagatorSettings",
                        R"doc(`PropagatorSettings`-derived class to define settings for single-arc dynamics.

)doc")
                        .def_property(
                            "termination_settings",
                            &tp::SingleArcPropagatorSettings<
                                double, TIME_TYPE>::getTerminationSettings,
                            &tp::SingleArcPropagatorSettings<
                                double, TIME_TYPE>::resetTerminationSettings,
                            get_docstring("SingleArcPropagatorSettings."
                                          "termination_settings")
                                .c_str())
                        .def_property(
                            "integrator_settings",
                            &tp::SingleArcPropagatorSettings<
                                double, TIME_TYPE>::getIntegratorSettings,
                            &tp::SingleArcPropagatorSettings<
                                double, TIME_TYPE>::setIntegratorSettings,
                            get_docstring("SingleArcPropagatorSettings."
                                          "termination_settings")
                                .c_str())
                        .def_property_readonly(
                            "processing_settings",
                            &tp::SingleArcPropagatorSettings<
                                double, TIME_TYPE>::getOutputSettings,
                            get_docstring("SingleArcPropagatorSettings."
                                          "processing_settings")
                                .c_str())
                        .def_property_readonly(
                            "print_settings",
                            &tp::SingleArcPropagatorSettings<
                                double, TIME_TYPE>::getPrintSettings,
                            get_docstring(
                                "SingleArcPropagatorSettings.print_settings")
                                .c_str());


                    py::class_<
                        tp::TranslationalStatePropagatorSettings<double,
                                                                 TIME_TYPE>,
                        std::shared_ptr<
                            tp::TranslationalStatePropagatorSettings<
                                double, TIME_TYPE>>,
                        tp::SingleArcPropagatorSettings<double, TIME_TYPE>>(
                        m, "TranslationalStatePropagatorSettings",
                        R"doc(`SingleArcPropagatorSettings`-derived class to define settings for single-arc translational dynamics.

)doc")

                        .def("get_propagated_state_size",
                             &tp::TranslationalStatePropagatorSettings<
                                 double, TIME_TYPE>::getPropagatedStateSize)
                        .def("reset_and_recreate_acceleration_models",
                             &tp::TranslationalStatePropagatorSettings<
                                 double, TIME_TYPE>::resetAccelerationModelsMap,
                             py::arg("new_acceleration_settings"),
                             py::arg("bodies"));


                    //            .def_property_readonly("acceleration_settings",
                    //                                   &tp::TranslationalStatePropagatorSettings<double>::getAccelerationSettingsMap);


                    py::class_<
                        tp::MultiTypePropagatorSettings<double, TIME_TYPE>,
                        std::shared_ptr<
                            tp::MultiTypePropagatorSettings<double, TIME_TYPE>>,
                        tp::SingleArcPropagatorSettings<double, TIME_TYPE>>(
                        m, "MultiTypePropagatorSettings",
                        R"doc(`SingleArcPropagatorSettings`-derived class to define settings for propagation of multiple quantities.

)doc")
                        .def("reset_initial_states",
                             &tp::MultiTypePropagatorSettings<
                                 double, TIME_TYPE>::resetInitialStates,
                             py::arg("initial_states"))
                        .def("recreate_state_derivative_models",
                             &tp::MultiTypePropagatorSettings<
                                 double, TIME_TYPE>::resetIntegratedStateModels,
                             py::arg("bodies"))
                        .def("single_type_settings",
                             &tp::MultiTypePropagatorSettings<
                                 double,
                                 TIME_TYPE>::getSingleTypePropagatorSettings,
                             py::arg("state_type"))
                        .def_property_readonly(
                            "propagator_settings_per_type",
                            &tp::MultiTypePropagatorSettings<
                                double, TIME_TYPE>::getPropagatorSettingsMap,
                            get_docstring("MultiTypePropagatorSettings."
                                          "propagator_settings_per_type")
                                .c_str());

                    py::class_<
                        tp::RotationalStatePropagatorSettings<double,
                                                              TIME_TYPE>,
                        std::shared_ptr<tp::RotationalStatePropagatorSettings<
                            double, TIME_TYPE>>,
                        tp::SingleArcPropagatorSettings<double, TIME_TYPE>>(
                        m, "RotationalStatePropagatorSettings",
                        R"doc(`SingleArcPropagatorSettings`-derived class to define settings for single-arc rotational state propagation.

)doc");

                    py::class_<
                        tp::MassPropagatorSettings<double, TIME_TYPE>,
                        std::shared_ptr<
                            tp::MassPropagatorSettings<double, TIME_TYPE>>,
                        tp::SingleArcPropagatorSettings<double, TIME_TYPE>>(
                        m, "MassPropagatorSettings",
                        get_docstring("MassPropagatorSettings").c_str());

                    py::class_<
                        tp::CustomStatePropagatorSettings<double, TIME_TYPE>,
                        std::shared_ptr<tp::CustomStatePropagatorSettings<
                            double, TIME_TYPE>>,
                        tp::SingleArcPropagatorSettings<double, TIME_TYPE>>(
                        m, "CustomStatePropagatorSettings",
                        get_docstring("CustomStatePropagatorSettings").c_str());


                    m.def(
                        "translational",
                        py::overload_cast<
                            const std::vector<std::string> &,
                            const tba::AccelerationMap &,
                            const std::vector<std::string> &,
                            const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                            const std::shared_ptr<
                                tp::PropagationTerminationSettings>,
                            const tp::TranslationalPropagatorType,
                            const std::vector<std::shared_ptr<
                                tp::SingleDependentVariableSaveSettings>> &,
                            const double>(
                            &tp::translationalStatePropagatorSettingsDeprecated<
                                double, TIME_TYPE>),
                        py::arg("central_bodies"),
                        py::arg("acceleration_models"),
                        py::arg("bodies_to_integrate"),
                        py::arg("initial_states"),
                        py::arg("termination_settings"),
                        py::arg("propagator") = tp::cowell,
                        py::arg("output_variables") =
                            std::vector<std::shared_ptr<
                                tp::SingleDependentVariableSaveSettings>>(),
                        py::arg("print_interval") = TUDAT_NAN);

                    m.def(
                        "translational",
                        tp::translationalStatePropagatorSettings<double,
                                                                 TIME_TYPE>,
                        py::arg("central_bodies"),
                        py::arg("acceleration_models"),
                        py::arg("bodies_to_integrate"),
                        py::arg("initial_states"), py::arg("initial_time"),
                        py::arg("integrator_settings"),
                        py::arg("termination_settings"),
                        py::arg("propagator") = tp::cowell,
                        py::arg("output_variables") =
                            std::vector<std::shared_ptr<
                                tp::SingleDependentVariableSaveSettings>>(),
                        py::arg("processing_settings") = std::make_shared<
                            tp::SingleArcPropagatorProcessingSettings>(),
                        R"doc(Factory function to create translational state propagator settings with stopping condition at given final time.

	Factory function to create translational state propagator settings for N bodies.
	The propagated state vector is defined by the combination of integrated bodies, and their central body, the combination
	of which define the relative translational states for which a differential equation is to be solved. The propagator
	input defines the formulation in which the differential equations are set up
	The dynamical models are defined by an ``AccelerationMap``, as created by :func:`~create_acceleration_models` function.
	Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/dynamics_types/translational.html>`_


	:param central_bodies:
		List of central bodies with respect to which the bodies to be integrated are propagated.
	:param acceleration_models:
		Set of accelerations acting on the bodies to propagate, provided as acceleration models.
	:param bodies_to_integrate:
		List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
	:param initial_states:
		Initial states of the bodies to integrate (one initial state for each body, concatenated into a single array), provided in the same order as the bodies to integrate. The initial states must be expressed in Cartesian elements, w.r.t. the central body of each integrated body. The states must be defined with the same frame orientation as the global frame orientation of the environment (specified when creating a system of bodies, see for instance :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings` and :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`). Consequently, for N integrated bodies, this input is a vector with size size 6N.
	:param initial_time:
		Initial epoch of the numerical propagation
	:param integrator_settings:
		Settings defining the numerical integrator that is to be used for the propagation

		.. note:: The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time

	:param termination_settings:
		Generic termination settings object to check whether the propagation should be ended.
	:param propagator:
		Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
	:param output_variables:
		Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
	:return:
		Translational state propagator settings object.
)doc");

                    m.def("rotational",
                          py::overload_cast<
                              const tba::TorqueModelMap &,
                              const std::vector<std::string> &,
                              const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                              const std::shared_ptr<
                                  tp::PropagationTerminationSettings>,
                              const tp::RotationalPropagatorType,
                              const std::vector<std::shared_ptr<
                                  tp::SingleDependentVariableSaveSettings>>,
                              const double>(
                              &tp::rotationalStatePropagatorSettingsDeprecated<
                                  double, TIME_TYPE>),
                          py::arg("torque_models"),
                          py::arg("bodies_to_integrate"),
                          py::arg("initial_states"),
                          py::arg("termination_settings"),
                          py::arg("propagator") = tp::quaternions,
                          py::arg("output_variables") =
                              std::vector<std::shared_ptr<
                                  tp::SingleDependentVariableSaveSettings>>(),
                          py::arg("print_interval") = TUDAT_NAN);

                    m.def(
                        "rotational",
                        &tp::rotationalStatePropagatorSettings<double,
                                                               TIME_TYPE>,
                        py::arg("torque_models"),
                        py::arg("bodies_to_integrate"),
                        py::arg("initial_states"), py::arg("initial_time"),
                        py::arg("integrator_settings"),
                        py::arg("termination_settings"),
                        py::arg("propagator") = tp::quaternions,
                        py::arg("output_variables") =
                            std::vector<std::shared_ptr<
                                tp::SingleDependentVariableSaveSettings>>(),
                        py::arg("processing_settings") = std::make_shared<
                            tp::SingleArcPropagatorProcessingSettings>(),
                        R"doc(Factory function to create rotational state propagator settings.

	Factory function to create rotational state propagator settings for N bodies.
	The propagated state vector is defined by the integrated bodies, which defines the bodies for which the
	differential equation defining the evolution of the rotational state between an
	inertial and body-fixed frame are to be solved. The propagator input defines the
	formulation in which the differential equations are set up. The dynamical models are
	defined by an ``TorqueModelMap``, as created by ``create_torque_models`` function.
	Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/dynamics_types/rotational.html>`_


	:param torque_models:
		Set of torques acting on the bodies to propagate, provided as torque models.
	:param bodies_to_integrate:
		List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
	:param initial_states:
		Initial rotational states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
		Regardless of the propagator that is selected, the initial rotational state is always defined as four quaternion entries, and the angular velocity of the body,
		as defined in more detail `here <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/environment_setup/use_of_reference_frames.html#definition-of-rotational-state>`_.

	:param termination_settings:
		Generic termination settings object to check whether the propagation should be ended.
	:param propagator:
		Type of rotational propagator to be used (see `RotationalPropagatorType` enum).
	:param output_variables:
		List of dependent variables to be saved (by default, no dependent variables are saved).
	:param print_interval:
		Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
	:return:
		Rotational state propagator settings object.
)doc");

                    m.def("mass",
                          py::overload_cast<
                              const std::vector<std::string>,
                              const std::map<std::string,
                                             std::vector<std::shared_ptr<
                                                 tba::MassRateModel>>> &,
                              const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                              const std::shared_ptr<
                                  tp::PropagationTerminationSettings>,
                              const std::vector<std::shared_ptr<
                                  tp::SingleDependentVariableSaveSettings>>,
                              const double>(
                              &tp::massPropagatorSettingsDeprecated<double,
                                                                    TIME_TYPE>),
                          py::arg("bodies_with_mass_to_propagate"),
                          py::arg("mass_rate_models"),
                          py::arg("initial_body_masses"),
                          py::arg("termination_settings"),
                          py::arg("output_variables") =
                              std::vector<std::shared_ptr<
                                  tp::SingleDependentVariableSaveSettings>>(),
                          py::arg("print_interval") = TUDAT_NAN);

                    m.def(
                        "mass", &tp::massPropagatorSettings<double, TIME_TYPE>,
                        py::arg("bodies_with_mass_to_propagate"),
                        py::arg("mass_rate_models"),
                        py::arg("initial_body_masses"), py::arg("initial_time"),
                        py::arg("integrator_settings"),
                        py::arg("termination_settings"),
                        py::arg("output_variables") =
                            std::vector<std::shared_ptr<
                                tp::SingleDependentVariableSaveSettings>>(),
                        py::arg("processing_settings") = std::make_shared<
                            tp::SingleArcPropagatorProcessingSettings>(),
                        R"doc(Factory function to create mass propagator settings

	Factory function to create mass propagator settings
	It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
	each body. In this function, the dependent variables to save are provided
	as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
	through the termination settings object provided.
	Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/dynamics_types/mass.html>`_


	:param bodies_with_mass_to_propagate:
		List of bodies whose mass should be numerically propagated.
	:param mass_rate_settings:
		Mass rates associated to each body, provided as a mass rate settings object.
	:param initial_body_masses:
		Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
	:param termination_settings:
		Generic termination settings object to check whether the propagation should be ended.
	:param output_variables:
		List of dependent variables to be saved (by default, no dependent variables are saved).
	:param print_interval:
		Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
	:return:
		Mass propagator settings object.
)doc");

                    m.def("custom_state",
                          &tp::customStatePropagatorSettings<double, TIME_TYPE>,
                          py::arg("state_derivative_function"),
                          py::arg("initial_state"), py::arg("initial_time"),
                          py::arg("integrator_settings"),
                          py::arg("termination_settings"),
                          py::arg("output_variables") =
                              std::vector<std::shared_ptr<
                                  tp::SingleDependentVariableSaveSettings>>(),
                          py::arg("processing_settings") = std::make_shared<
                              tp::SingleArcPropagatorProcessingSettings>(),
                          get_docstring("custom_state").c_str());


                    m.def(
                        "multitype",
                        py::overload_cast<
                            const std::vector<
                                std::shared_ptr<tp::SingleArcPropagatorSettings<
                                    double, TIME_TYPE>>>,
                            const std::shared_ptr<
                                tp::PropagationTerminationSettings>,
                            const std::vector<std::shared_ptr<
                                tp::SingleDependentVariableSaveSettings>>,
                            const double>(
                            &tp::multiTypePropagatorSettingsDeprecated<
                                double, TIME_TYPE>),
                        py::arg("propagator_settings_list"),
                        py::arg("termination_settings"),
                        py::arg("output_variables") =
                            std::vector<std::shared_ptr<
                                tp::SingleDependentVariableSaveSettings>>(),
                        py::arg("print_interval") = TUDAT_NAN);

                    m.def(
                        "multitype",
                        &tp::multiTypePropagatorSettings<double, TIME_TYPE>,
                        py::arg("propagator_settings_list"),
                        py::arg("integrator_settings"), py::arg("initial_time"),
                        py::arg("termination_settings"),
                        py::arg("output_variables") =
                            std::vector<std::shared_ptr<
                                tp::SingleDependentVariableSaveSettings>>(),
                        py::arg("processing_settings") = std::make_shared<
                            tp::SingleArcPropagatorProcessingSettings>(),
                        R"doc(Factory function to create multitype propagator settings.

	Factory function to create multitype propagator settings.
	It works by providing a list of SingleArcPropagatorSettings objects. When using this function,
	only the termination and output settings provided here are used, any such settings in the
	constituent propagator settings are ignored
	Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/dynamics_types/multi_type.html>`_

	.. note:: The propagated state contains the state types in the following order: Translational ( **C** ), Rotational ( **R** ), Mass ( **M** ), and Custom ( **C** ).
	          When propagating two bodies, an example of what the output state would look like is for instance:
	          [ **T** Body 1, **T** Body 2, **R** Body 1, **R** Body 2, **M** Body 1, **M** Body 2 ]


	:param propagator_settings_list:
		List of SingleArcPropagatorSettings objects to use.
	:param termination_settings:
		Generic termination settings object to check whether the propagation should be ended.
	:param output_variables:
		List of dependent variables to be saved (by default, no dependent variables are saved).
	:param print_interval:
		Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
	:return:
		Mass propagator settings object.
)doc");


                    m.def(
                        "multi_arc",
                        &tp::multiArcPropagatorSettings<double, TIME_TYPE>,
                        py::arg("single_arc_settings"),
                        py::arg("transfer_state_to_next_arc") = false,
                        py::arg("processing_settings") = std::make_shared<
                            tp::MultiArcPropagatorProcessingSettings>(),
                        R"doc(Factory function to create multi-arc propagator settings.

	Factory function to create multi-arc propagator settings. It works by providing separate settings for
	each arc in a list.


	:param single_arc_settings:
		List of SingleArcPropagatorSettings objects to use, one for each arc.
	:param transfer_state_to_next_arc:
		It denotes whether whether the initial state of arc N+1 is to be taken from arc N (for N>0).
	:return:
		Multi-arc propagator settings object.
)doc");

                    m.def(
                        "hybrid_arc",
                        &tp::hybridArcPropagatorSettings<double, TIME_TYPE>,
                        py::arg("single_arc_settings"),
                        py::arg("multi_arc_settings"),
                        py::arg("processing_settings") = std::make_shared<
                            tp::HybridArcPropagatorProcessingSettings>(),
                        R"doc(Factory function to create hybrid-arc propagator settings.

	Factory function to create hybrid-arc propagator settings (i.e., a combination of single- and multi-arc dynamics).


	:param single_arc_settings:
		SingleArcPropagatorSettings object to use for the propagation.
	:param multi_arc_settings:
		MultiArcPropagatorSettings object to use for the propagation.
	:return:
		Hybrid-arc propagator settings object.
)doc");

                    ///////////////////////////////////////////////////////////////////////////////////////


                    py::class_<
                        tp::PropagationTerminationSettings,
                        std::shared_ptr<tp::PropagationTerminationSettings>>(
                        m, "PropagationTerminationSettings",
                        R"doc(Functional base class to define termination settings for the propagation.

)doc");

                    py::class_<
                        tp::PropagationDependentVariableTerminationSettings,
                        std::shared_ptr<
                            tp::PropagationDependentVariableTerminationSettings>,
                        tp::PropagationTerminationSettings>(
                        m, "PropagationDependentVariableTerminationSettings",
                        R"doc(`PropagationTerminationSettings`-derived class to define termination settings for the propagation from dependent variables.

)doc");

                    py::class_<
                        tp::PropagationTimeTerminationSettings,
                        std::shared_ptr<tp::PropagationTimeTerminationSettings>,
                        tp::PropagationTerminationSettings>(
                        m, "PropagationTimeTerminationSettings",
                        R"doc(`PropagationTerminationSettings`-derived class to define termination settings for the propagation from propagation time.

)doc");

                    py::class_<tp::PropagationCPUTimeTerminationSettings,
                               std::shared_ptr<
                                   tp::PropagationCPUTimeTerminationSettings>,
                               tp::PropagationTerminationSettings>(
                        m, "PropagationCPUTimeTerminationSettings",
                        R"doc(`PropagationTerminationSettings`-derived class to define termination settings for the propagation from CPU time.

)doc");

                    py::class_<tp::PropagationCustomTerminationSettings,
                               std::shared_ptr<
                                   tp::PropagationCustomTerminationSettings>,
                               tp::PropagationTerminationSettings>(
                        m, "PropagationCustomTerminationSettings",
                        R"doc(`PropagationTerminationSettings`-derived class to define custom termination settings for the propagation.

)doc");

                    py::class_<tp::PropagationHybridTerminationSettings,
                               std::shared_ptr<
                                   tp::PropagationHybridTerminationSettings>,
                               tp::PropagationTerminationSettings>(
                        m, "PropagationHybridTerminationSettings",
                        R"doc(`PropagationTerminationSettings`-derived class to define hybrid termination settings for the propagation.

)doc");

                    py::class_<
                        tp::NonSequentialPropagationTerminationSettings,
                        std::shared_ptr<
                            tp::NonSequentialPropagationTerminationSettings>,
                        tp::PropagationTerminationSettings>(
                        m, "NonSequentialPropagationTerminationSettings",
                        get_docstring(
                            "NonSequentialPropagationTerminationSettings")
                            .c_str());

                    //                .def(py::init<
                    //                             const
                    //                             std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
                    //                             const double,
                    //                             const bool,
                    //                             const bool,
                    //                             const
                    //                             std::shared_ptr<tudat::root_finders::RootFinderSettings>>(),
                    //                     py::arg("dependent_variadble_settings"),
                    //                     py::arg("limit_value"),
                    //                     py::arg("use_as_lower_limit"),
                    //                     py::arg("terminate_exactly_on_final_condition")
                    //                     = false,
                    //                     py::arg("termination_root_finder_settings")
                    //                     = nullptr);

                    m.def(
                        "time_termination",
                        &tp::propagationTimeTerminationSettings,
                        py::arg("termination_time"),
                        py::arg("terminate_exactly_on_final_condition") = false,
                        R"doc(Factory function to create time termination settings for the propagation.

	Factory function to create time termination settings for the propagation.
	The propagation is stopped when the final time provided is reached. Note that the termination time is set as the
	absolute time (in seconds since J2000), not the time since the start of the propagation.
	Depending on the sign of the time step of the numerical integrator, the termination time will be treated as an
	upper bound (for positive time step) or lower bound (for negative time step).
	The simulator will normally finish the final time-step, which may cause the termination time to be slightly exceeded.
	This behaviour can be suppressed by providing the optional input argument
	``terminate_exactly_on_final_condition=True``, in which case the final propagation step will be *exactly* on the
	specified time.


	:param termination_time:
		Final time of the propagation.
	:param terminate_exactly_on_final_condition:
		Denotes whether the propagation is to terminate exactly on the final condition, or whether it is to terminate on the first step where it is violated.
	:return:
		Time termination settings object.
)doc");

                    m.def(
                        "cpu_time_termination",
                        &tp::propagationCPUTimeTerminationSettings,
                        py::arg("cpu_termination_time"),
                        R"doc(Factory function to create CPU time termination settings for the propagation.

	Factory function to create CPU time termination settings for the propagation.
	The propagation is stopped when the final CPU time provided is reached.


	:param cpu_termination_time:
		Maximum CPU time for the propagation.
	:return:
		CPU time termination settings object.
)doc");

                    m.def(
                        "dependent_variable_termination",
                        &tp::propagationDependentVariableTerminationSettings,
                        py::arg("dependent_variable_settings"),
                        py::arg("limit_value"), py::arg("use_as_lower_limit"),
                        py::arg("terminate_exactly_on_final_condition") = false,
                        py::arg("termination_root_finder_settings") = nullptr,
                        R"doc(Factory function to create termination settings for the propagation based on a dependent variable.

	Factory function to create termination settings for the propagation based on the value of a dependent variable.
	The propagation is stopped when a provided upper or lower limit value is reached.
	The simulator will normally finish the final time-step, which may cause the dependent variable to be slightly exceeded.
	This behaviour can be suppressed by providing the optional input argument
	``terminate_exactly_on_final_condition=True``, in which case the final propagation step will be *exactly* on the
	specified dependent variable value.


	:param dependent_variable_settings:
		Dependent variable object to be used as termination setting.
	:param limit_value:
		Limit value of the dependent variable; if reached, the propagation is stopped.
	:param use_as_lower_limit:
		Denotes whether the limit value should be used as lower or upper limit.
	:param terminate_exactly_on_final_condition:
		Denotes whether the propagation is to terminate exactly on the final condition, or whether it is to terminate on the first step where it is violated.
	:param termination_root_finder_settings:
		Settings object to create root finder used to converge on exact final condition.
	:return:
		Dependent variable termination settings object.
)doc");

                    m.def(
                        "custom_termination",
                        &tp::popagationCustomTerminationSettings,
                        py::arg("custom_condition"),
                        R"doc(Factory function to create custom termination settings for the propagation.

	Factory function to create custom termination settings for the propagation.
	The propagation is stopped when the condition provided is verified.
	This custom function should take the current time as input and output a Boolean. It can use internal variables
	and calculations, for example retrieved from the environment.


	:param custom_condition:
		Function of time (independent variable) which is called during the propagation and returns a boolean value denoting whether the termination condition is verified.
	:return:
		Custom termination settings object.
)doc");

                    m.def(
                        "hybrid_termination",
                        &tp::propagationHybridTerminationSettings,
                        py::arg("termination_settings"),
                        py::arg("fulfill_single_condition"),
                        R"doc(Factory function to create hybrid termination settings for the propagation.

	Factory function to create hybrid termination settings for the propagation. This function can be used
	to define that all conditions or a single condition of the conditions provided must be met to
	stop the propagation. Each termination condition should be created according to each individual factory function
	and then added to a list of termination conditions.


	:param termination_settings:
		List of single PropagationTerminationSettings objects to be checked during the propagation.
	:param fulfill_single_condition:
		Whether only a single condition of those provided must be met to stop the propagation (true) or all of them simultaneously (false).
	:return:
		Hybrid termination settings object.
)doc");

                    m.def("non_sequential_termination",
                          &tp::nonSequentialPropagationTerminationSettings,
                          py::arg("forward_termination_settings"),
                          py::arg("backward_termination_settings"),
                          get_docstring("non_sequential_termination").c_str());

                    m.def(
                        "add_dependent_variable_settings",
                        &tp::addDepedentVariableSettings<double>,
                        py::arg("dependent_variable_settings"),
                        py::arg("propagator_settings"),
                        R"doc(Function to add dependent variables to existing propagator settings.

	Function to add dependent variables to existing :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.SingleArcPropagatorSettings`
	object. This function is added as an alternative to teh regular manner in which to defined dependent variables (use of input to factory
	functions for single-arc propagator settings :func:`~tudatpy.numerical_simulation.propagation_setup.propagator.translational`,
	:func:`~tudatpy.numerical_simulation.propagation_setup.propagator.rotational`, :func:`~tudatpy.numerical_simulation.propagation_setup.propagator.mass`,
	:func:`~tudatpy.numerical_simulation.propagation_setup.propagator.multitype`). Typically, this function is used to modify
	existing propagator settings in a loop when running multiple simulations


	:param dependent_variable_settings:
		List of dependent variable settings that are to be added to propagator settings. Note that this function adds settings, and does not replace any existing settings (nor does it check for duplicate settings).
	:param propagator_settings:
		Propagator settings to which the additional dependent variables settings are to be added.
	:return:
)doc");
                }

            }  // namespace propagator
        }  // namespace propagation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
