/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_propagator.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/propagators/getZeroProperModeRotationalInitialState.h>
#include <tudat/simulation/propagation_setup.h>

#include "scalarTypes.h"

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudatpy
{
namespace dynamics
{
namespace propagation_setup
{
namespace propagator
{

std::shared_ptr< tudat::propagators::MultiArcPropagatorProcessingSettings > multiArcProcessingSettings( )
{
    return std::make_shared< tudat::propagators::MultiArcPropagatorProcessingSettings >( );
}

void expose_propagator_setup( py::module &m )
{
    ///////////////////////////////////////////////////////////////////////////////////////

    // ENUMS
    py::enum_< tp::TranslationalPropagatorType >( m,
                                                  "TranslationalPropagatorType",
                                                  R"doc(

Enumeration of available translational propagator types.





      )doc" )
            .value( "undefined_translational_propagator",
                    tp::TranslationalPropagatorType::undefined_translational_propagator,
                    R"doc(
      )doc" )
            .value( "cowell",
                    tp::TranslationalPropagatorType::cowell,
                    R"doc(

Propagation of Cartesian elements (state vector size 6), without any transformations

)doc" )
            .value( "encke",
                    tp::TranslationalPropagatorType::encke,
                    R"doc(

Propagation of the difference in Cartesian elements of the orbit w.r.t. an unperturbed reference orbit. The reference orbit is generated from the initial state/central body, and not updated during the propagation (see :cite:t:`wakker2015`)

)doc" )
            .value( "gauss_keplerian",
                    tp::TranslationalPropagatorType::gauss_keplerian,
                    R"doc(

Propagation of Keplerian elements (state vector size 6), with true anomaly as the 'fast' element  (see :cite:t:`vallado2001`)

)doc" )
            .value( "gauss_modified_equinoctial",
                    tp::TranslationalPropagatorType::gauss_modified_equinoctial,
                    R"doc(

Propagation of Modified equinoctial elements (state vector size 6), with the element :math:`I` defining the location of the singularity based on the initial condition (see :cite:t:`hintz2008`)

)doc" )
            .value( "unified_state_model_quaternions",
                    tp::TranslationalPropagatorType::unified_state_model_quaternions,
                    R"doc(

Propagation of Unified state model using quaternions (state vector size 7, see :cite:t:`vittaldev2012`)

)doc" )
            .value( "unified_state_model_modified_rodrigues_parameters",
                    tp::TranslationalPropagatorType::unified_state_model_modified_rodrigues_parameters,
                    R"doc(

Propagation of Unified state model using modified Rodrigues parameters (state vector size 7, last element represents shadow parameter, see :cite:t:`vittaldev2012`)

)doc" )
            .value( "unified_state_model_exponential_map",
                    tp::unified_state_model_exponential_map,
                    R"doc(

Propagation of Unified state model using exponential map (state vector size 7, last element represents shadow parameter, see :cite:t:`vittaldev2012`)

)doc" )
            .export_values( );

    py::enum_< tp::RotationalPropagatorType >( m,
                                               "RotationalPropagatorType",
                                               R"doc(

Enumeration of available rotational propagator types.



      )doc" )
            .value( "undefined_rotational_propagator",
                    tp::RotationalPropagatorType::undefined_rotational_propagator,
                    R"doc(
      )doc" )
            .value( "quaternions",
                    tp::RotationalPropagatorType::quaternions,
                    R"doc(

Entries 1-4: The quaternion defining the rotation from inertial to body-fixed frame (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html#definition-of-rotational-state>`__)

Entries 5-7: The body's angular velocity vector, expressed in its body-fixed frame.

)doc" )
            .value( "modified_rodrigues_parameters",
                    tp::RotationalPropagatorType::modified_rodrigues_parameters,
                    R"doc(

Entries 1-4: The modified Rodrigues parameters defining the rotation from inertial to body-fixed frame (with entry four the shadow parameter)

Entries 5-7: The body's angular velocity vector, expressed in its body-fixed frame.

)doc" )
            .value( "exponential_map",
                    tp::RotationalPropagatorType::exponential_map,
                    R"doc(

Entries 1-4: The exponential map defining the rotation from inertial to body-fixed frame (with entry four the shadow parameter)

Entries 5-7: The body's angular velocity vector, expressed in its body-fixed frame.

)doc" )
            .export_values( );

    py::enum_< tp::PropagationTerminationTypes >( m,
                                                  "PropagationTerminationTypes",
                                                  R"doc(

Enumeration of possible propagation termination types





      )doc" )
            .value( "time_stopping_condition_type",
                    tp::PropagationTerminationTypes::time_stopping_condition,
                    R"doc(No propagator documentation found.)doc" )
            .value( "cpu_time_stopping_condition_type",
                    tp::PropagationTerminationTypes::cpu_time_stopping_condition,
                    R"doc(No propagator documentation found.)doc" )
            .value( "dependent_variable_stopping_condition_type",
                    tp::PropagationTerminationTypes::dependent_variable_stopping_condition,
                    R"doc(No propagator documentation found.)doc" )
            .value( "hybrid_stopping_condition_type",
                    tp::PropagationTerminationTypes::hybrid_stopping_condition,
                    R"doc(No propagator documentation found.)doc" )
            .value( "custom_stopping_condition_type",
                    tp::PropagationTerminationTypes::custom_stopping_condition,
                    R"doc(No propagator documentation found.)doc" )
            .export_values( );

    py::enum_< tp::IntegratedStateType >( m,
                                          "StateType",
                                          R"doc(

Enumeration of available integrated state types.





      )doc" )
            .value( "hybrid_type",
                    tp::IntegratedStateType::hybrid,
                    R"doc(
      )doc" )
            .value( "translational_type",
                    tp::IntegratedStateType::translational_state,
                    R"doc(
      )doc" )
            .value( "rotational_type",
                    tp::IntegratedStateType::rotational_state,
                    R"doc(
      )doc" )
            .value( "mass_type",
                    tp::IntegratedStateType::body_mass_state,
                    R"doc(
      )doc" )
            .value( "custom_type",
                    tp::IntegratedStateType::custom_state,
                    R"doc(
      )doc" )
            .export_values( );

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

    py::class_< tp::PropagationPrintSettings, std::shared_ptr< tp::PropagationPrintSettings > >( m,
                                                                                                 "PropagationPrintSettings",
                                                                                                 R"doc(

         Class to save settings on what is to be written to the console during the propagation of a single arc.

         Upon creation, this object has default settings such that no data is printed to the console.

      )doc" )
            .def_property( "print_dependent_variable_indices",
                           &tp::PropagationPrintSettings::getPrintDependentVariableData,
                           &tp::PropagationPrintSettings::setPrintDependentVariableData,
                           R"doc(

         Boolean defining whether the meaning and indices of the
         entries of the dependent variable data are to be printed to
         the console (before the propagation).

         .. note:: The same information can be retrieved from the
                   :py:attr:`SingleArcSimulationResults.dependent_variable_ids`
                   attribute.


         :type: bool
      )doc" )
            .def_property( "print_state_indices",
                           &tp::PropagationPrintSettings::getPrintPropagatedStateData,
                           &tp::PropagationPrintSettings::setPrintPropagatedStateData,
                           R"doc(

         Boolean defining whether the meaning and indices of the
         entries of the state vector are to be printed to
         the console (before the propagation).

         .. note:: The same information can be retrieved from the
                   :py:attr:`SingleArcSimulationResults.propagated_state_ids`
                   attribute.


         :type: bool
      )doc" )
            .def_property( "print_processed_state_indices",
                           &tp::PropagationPrintSettings::getPrintProcessedStateData,
                           &tp::PropagationPrintSettings::setPrintProcessedStateData,
                           R"doc(

         Boolean defining whether the meaning and indices of the
         entries of the processed state vector are to be printed to
         the console (after the propagation). The distinction between the
         propagated and processed (or conventional) state representation is described in
         detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/processed_propagated_elements.html>`__.
         Summarizing: the processed state is the 'typical' formulation of the state (for translational dynamics: Cartesian states).

         .. note::

             The same information can be retrieved from the :py:attr:`SingleArcSimulationResults.processed_state_ids` attribute.

         :type: bool
      )doc" )
            .def_property( "print_number_of_function_evaluations",
                           &tp::PropagationPrintSettings::getPrintNumberOfFunctionEvaluations,
                           &tp::PropagationPrintSettings::setPrintNumberOfFunctionEvaluations,
                           R"doc(

         Boolean defining whether the number of function evaluations that
         were performed is to be printed to the console (after propagation).


         :type: bool
      )doc" )
            .def_property( "print_propagation_clock_time",
                           &tp::PropagationPrintSettings::getPrintPropagationTime,
                           &tp::PropagationPrintSettings::setPrintPropagationTime,
                           R"doc(

         Boolean defining whether the total clock time taken for the propagation
         is to be printed to the console (after propagation).


         :type: bool
      )doc" )
            .def_property( "print_termination_reason",
                           &tp::PropagationPrintSettings::getPrintTerminationReason,
                           &tp::PropagationPrintSettings::setPrintTerminationReason,
                           R"doc(

         Boolean defining whether the reason for propagation termination
         is to be printed to the console (after propagation).


         :type: bool
      )doc" )
            .def_property( "print_initial_and_final_conditions",
                           &tp::PropagationPrintSettings::getPrintInitialAndFinalConditions,
                           &tp::PropagationPrintSettings::setPrintInitialAndFinalConditions,
                           R"doc(

         Boolean defining whether the initial and final conditions (state and time)
         are to be printed to the console (before and after propagation, respectively).


         :type: bool
      )doc" )
            .def_property( "results_print_frequency_in_seconds",
                           &tp::PropagationPrintSettings::getResultsPrintFrequencyInSeconds,
                           &tp::PropagationPrintSettings::setResultsPrintFrequencyInSeconds,
                           R"doc(

         Variable indicating how often (in seconds of simulation time)
         the current state and time are to be printed to the console (by default, set to NaN - they are never printed).
         In case this setting is active (e.g. not NaN), and the ``results_print_frequency_in_steps`` setting is active,
         the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.


         :type: Float
      )doc" )
            .def_property( "results_print_frequency_in_steps",
                           &tp::PropagationPrintSettings::getResultsPrintFrequencyInSteps,
                           &tp::PropagationPrintSettings::setResultsPrintFrequencyInSteps,
                           R"doc(

         Variable indicating how often (in number of full integration steps)
         the current state and time are to be printed to the console (by default, set to 0 - they are never printed).
         In case this setting is active (e.g. not 0), and the ``results_print_frequency_in_seconds`` setting is active,
         the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.


         :type: int
      )doc" )
            .def_property( "print_dependent_variables_during_propagation",
                           &tp::PropagationPrintSettings::getPrintDependentVariableDuringPropagation,
                           &tp::PropagationPrintSettings::setPrintDependentVariableDuringPropagation,
                           R"doc(

         Boolean defining whether the dependent variables are to be printed during the propagation along with the state,
         at steps/epochs define by the ``results_print_frequency_in_seconds`` and/or ``results_print_frequency_in_steps`` inputs.


         :type: float
      )doc" )
            .def( "enable_all_boolean_printing",
                  py::overload_cast<>( &tp::PropagationPrintSettings::enableAllPrinting ),
                  R"doc(

         Function enabling all True/False printing (e.g. sets all boolean attributes to True)






     )doc" )
            .def( "enable_all_printing",
                  py::overload_cast< const double, const int >( &tp::PropagationPrintSettings::enableAllPrinting ),
                  py::arg( "results_print_frequency_in_seconds" ),
                  py::arg( "results_print_frequency_in_steps" ),
                  R"doc(

         Function enabling all True/False printing (e.g. sets all boolean attributes to True), and setting the non-boolean
         attributes to values defined here.




         Parameters
         ----------
         results_print_frequency_in_seconds : float
             See ``results_print_frequency_in_seconds`` class attribute
         results_print_frequency_in_steps : int
             See ``results_print_frequency_in_steps`` class attribute




     )doc" )
            .def( "disable_all_printing",
                  &tp::PropagationPrintSettings::disableAllPrinting,
                  R"doc(

         Function enabling all printing (e.g. sets all boolean attributes to False, and disables all other output as well)






     )doc" );

    py::class_< tp::PropagatorProcessingSettings, std::shared_ptr< tp::PropagatorProcessingSettings > >( m,
                                                                                                         "PropagatorProcessingSettings",
                                                                                                         R"doc(

         Base class to define settings on how the numerical results are to be used

         Base class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
         Instances of this class are typically not created by the user. Settings objects for derived class of single-, multi- and hybrid arc propagation are
         instantiated through the functions to define propagator settings (such as :func:`~translational` or :func:`~multi_arc`) in this module.

         Upon creation, this object has default settings, ``clear_numerical_solution: false`` and ``set_integrated_result: false``.

      )doc" )
            .def_property( "set_integrated_result",
                           &tp::PropagatorProcessingSettings::getSetIntegratedResult,
                           &tp::PropagatorProcessingSettings::setIntegratedResult,
                           R"doc(

         Boolean defining whether the propagation results are to
         be used to update the environment. If this variable is set
         to False, the numerical propagation results can be
         retrieved from this object (provided the
         :py:attr:`~clear_numerical_solution` is set to False),
         but the (for instance) Ephemeris of the propagated body
         is not updated with the propagation results. If this
         variable is set to True, the properties of the propagated
         :class:`~tudatpy.dynamics.environment.Body`
         object will be updated as per the numerical results
         (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/printing_processing_results.html#automatic-processing>`__ for details).


         :type: bool
      )doc" )
            .def_property( "clear_numerical_solution",
                           &tp::PropagatorProcessingSettings::getClearNumericalSolutions,
                           &tp::PropagatorProcessingSettings::setClearNumericalSolutions,
                           R"doc(

         Boolean defining whether the propagation results should be
         deleted after the propagation is terminated. If this is
         done, the :py:attr:`~state_history`,
         :py:attr:`~unprocessed_state_history` and
         :py:attr:`~dependent_variable_history` will not be
         available in the :class:`~tudatpy.dynamics.propagation.SingleArcSimulationResults` class. Putting this setting to True (deleting the
         results) is only sensible when the
         :py:attr:`~set_integrated_result` is set to True. In that
         case, the propagated states are *not* accessible directly
         from this objects, but the results are used to update the
         environment, *e.g.* update the ephemeris of the propagated
         body with the numerical results.


         :type: bool
      )doc" )
            .def_property( "create_dependent_variable_interface",
                           &tp::PropagatorProcessingSettings::getUpdateDependentVariableInterpolator,
                           &tp::PropagatorProcessingSettings::setUpdateDependentVariableInterpolator,
                           R"doc(No propagator documentation found.)doc" );

    py::class_< tp::SingleArcPropagatorProcessingSettings,
                std::shared_ptr< tp::SingleArcPropagatorProcessingSettings >,
                tp::PropagatorProcessingSettings >( m,
                                                    "SingleArcPropagatorProcessingSettings",
                                                    R"doc(

         Class to define settings on how the numerical results are to be used for single-arc propagations

         Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
         derived from :class:`PropagatorProcessingSettings`.
         Instances of this class are typically not created by the user. A settings object is
         instantiated through the functions to define single-arc propagator settings (such as :func:`~translational` or :func:`~rotational`) in this module

         Upon creation, this object has default settings, ``print_settings: `` (see :class:`~PropagationPrintSettings` for default), ``results_save_frequency_in_steps:`` 1,
         and ``results_save_frequency_in_seconds:`` NaN (e.g. not used).


      )doc" )
            .def_property_readonly( "print_settings",
                                    &tp::SingleArcPropagatorProcessingSettings::getPrintSettings,
                                    R"doc(

         Settings object defining which quantities should be printed to the console before, during and after the propagation. By default, this
         object is instantiated to print nothing.


         :type: PropagationPrintSettings
      )doc" )
            .def_property( "results_save_frequency_in_steps",
                           &tp::SingleArcPropagatorProcessingSettings::getResultsSaveFrequencyInSteps,
                           &tp::SingleArcPropagatorProcessingSettings::setResultsSaveFrequencyInSteps,
                           R"doc(

         Variable indicating how often (in number of integrator steps)
         the propagated time, state, dependent variables, etc. are to be saved to data structures containing the results
         (by default, set to 1 - they are saved every time step). If this setting is set to 0, the data is never saved based on number of steps.
         In case this setting is active (e.g. not 0), and the ``results_save_frequency_in_seconds`` setting is active,
         the data is saved as soon as *one* of the two conditions (number of seconds, or number of steps) is met.


         :type: int
      )doc" )
            .def_property( "results_save_frequency_in_seconds",
                           &tp::SingleArcPropagatorProcessingSettings::getResultsSaveFrequencyInSeconds,
                           &tp::SingleArcPropagatorProcessingSettings::setResultsSaveFrequencyInSeconds,
                           R"doc(

         Variable indicating how often (in seconds of simulation time)
         the propagated time, state, dependent variables, etc. are to be saved to data structures containing the results
         (by default, set to NaN - they are not saved based on propagation time; see below and ``results_save_frequency_in_steps`` attribute ).
         In case this setting is active, and set to :math:`\Delta t`, the data are saved as soon as the current time step is :math:`\ge \Delta t` after the
         last step at which data was saved.
         In case this setting is active (e.g. not NaN), and the ``results_save_frequency_in_steps`` setting is active,
         the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.


         :type: float
      )doc" );

    py::class_< tp::MultiArcPropagatorProcessingSettings,
                std::shared_ptr< tp::MultiArcPropagatorProcessingSettings >,
                tp::PropagatorProcessingSettings >( m,
                                                    "MultiArcPropagatorProcessingSettings",
                                                    R"doc(

         Class to define settings on how the numerical results are to be used for multi-arc propagations

         Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
         derived from :class:`PropagatorProcessingSettings`.
         Instances of this class are typically not created by the user. A settings object is
         instantiated through the function :func:`~multi_arc` to define multi-arc propagator settings.
         This object contains a list of :class:`SingleArcPropagatorProcessingSettings` objects, containing the processing settings for each constituent arc.





      )doc" )
            .def( "set_print_settings_for_all_arcs",
                  &tp::MultiArcPropagatorProcessingSettings::resetAndApplyConsistentSingleArcPrintSettings,
                  py::arg( "single_arc_print_settings" ),
                  R"doc(

         Function that sets the same print settings for each arc in the multi-arc propagation.



         Parameters
         ----------
         single_arc_print_settings : PropagationPrintSettings
             Propagation print settings that are applied to each constituent single-arc settings, overriding any existing settings.




     )doc" )
            .def_property( "print_output_on_first_arc_only",
                           &tp::MultiArcPropagatorProcessingSettings::getPrintFirstArcOnly,
                           &tp::MultiArcPropagatorProcessingSettings::resetPrintFirstArcOnly,
                           R"doc(

         **read-only**

         Variable defining whether the ``set_print_settings_for_all_arcs`` function has been used to define identical print settings for each arc.


         :type: bool
      )doc" )
            .def_property_readonly( "identical_settings_per_arc",
                                    &tp::MultiArcPropagatorProcessingSettings::useIdenticalSettings,
                                    R"doc(No propagator documentation found.)doc" )
            .def_property_readonly( "single_arc_settings",
                                    &tp::MultiArcPropagatorProcessingSettings::getSingleArcSettings,
                                    R"doc(

         **read-only**

         List containing the processing settings for each constituent arc


         :type: list[SingleArcPropagatorProcessingSettings]
      )doc" );

    py::class_< tp::HybridArcPropagatorProcessingSettings,
                std::shared_ptr< tp::HybridArcPropagatorProcessingSettings >,
                tp::PropagatorProcessingSettings >( m,
                                                    "HybridArcPropagatorProcessingSettings",
                                                    R"doc(

         Class to define settings on how the numerical results are to be used for hybrid-arc propagations

         Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
         derived from :class:`PropagatorProcessingSettings`.
         Instances of this class are typically not created by the user. A settings object is
         instantiated through the function :func:`~hybrid_arc` to define hybrid-arc propagator settings.
         This object contains a :class:`SingleArcPropagatorProcessingSettings` object and a :class:`MultiArcPropagatorProcessingSettings` ,
         containing the processing settings for the constituents of the hybrid-arc propagation.





      )doc" )
            .def_property( "set_integrated_result",
                           &tp::HybridArcPropagatorProcessingSettings::getSetIntegratedResult,
                           &tp::HybridArcPropagatorProcessingSettings::setIntegratedResult )
            .def_property( "clear_numerical_solution",
                           &tp::HybridArcPropagatorProcessingSettings::getClearNumericalSolutions,
                           &tp::HybridArcPropagatorProcessingSettings::setClearNumericalSolutions )
            .def( "set_print_settings_for_all_arcs",
                  &tp::HybridArcPropagatorProcessingSettings::resetAndApplyConsistentPrintSettings,
                  py::arg( "print_settings" ),
                  R"doc(

         Function that sets the same print settings for each arc in the multi-arc propagation.



         Parameters
         ----------
         single_arc_print_settings : PropagationPrintSettings
             Propagation print settings that are applied to each constituent single-arc settings, overriding any existing settings.




     )doc" )
            .def_property_readonly( "single_arc_settings",
                                    &tp::HybridArcPropagatorProcessingSettings::getSingleArcSettings,
                                    R"doc(

         Processing settings for the single-arc component of the hybrid-arc propagation.


         :type: SingleArcPropagatorProcessingSettings
      )doc" )
            .def_property_readonly( "multi_arc_settings",
                                    &tp::HybridArcPropagatorProcessingSettings::getMultiArcSettings,
                                    R"doc(

         Processing settings for the single-arc component of the multi-arc propagation.


         :type: MultiArcPropagatorProcessingSettings
      )doc" );

    // CLASSES
    py::class_< tp::PropagatorSettings< STATE_SCALAR_TYPE >, std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > > >(
            m,
            "PropagatorSettings",
            R"doc(

         Base class to define settings for propagators.

         Base class to define settings for propagators. Derived classes are split into settings for single-, multi-arc and hybrid-arc dynamics.
         Single-arc settings are in turn split into specific derived classes per dynamics type


      )doc" )
            .def_property( "initial_states",
                           &tp::PropagatorSettings< STATE_SCALAR_TYPE >::getInitialStates,
                           &tp::PropagatorSettings< STATE_SCALAR_TYPE >::resetInitialStates,
                           R"doc(

            Vector of initial state values for the numerical propagation. For the multi-arc, this contains the
            concatenated single-arc initial states. For the hybrid-arc, this contains first the list of single-arc component initial states,
            followed by the multi-arc component initial states.

            :type: np.array
)doc" );

    py::class_< tp::MultiArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::MultiArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::PropagatorSettings< STATE_SCALAR_TYPE > >( m,
                                                               "MultiArcPropagatorSettings",
                                                               R"doc(

            Class derived from :class:`PropagatorSettings` to define settings for multi-arc dynamics

            Class derived from :class:`PropagatorSettings` to define settings for multi-arc dynamics
            An object of this type is typically created using the dedicated creation function :func:`~multi_arc`


      )doc" )
            .def_property_readonly( "processing_settings",
                                    &tp::MultiArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getOutputSettings,
                                    R"doc(
            **read-only**

            Settings that determine how the multi-arc propagation results are processed (e.g. if the results are used to update the body ephemeris, if
            data from each epoch is saved, etc.), and which data is printed to the console
            See `user guide <https://docs.tudat.space/en/latest/user-guide/state-propagation/propagation-setup/printing-processing-results.html#automatic-processing>`_ for more details.

            :type: MultiArcPropagatorProcessingSettings

)doc" )
        .def_property_readonly( "initial_state_list",
                                &tp::MultiArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getInitialStateList,
                                R"doc(
            **read-only**

            List of initial states per arc (e.g. entry j of this list is the initial state for arc j).

            :type: list[np.array]

)doc" )

        .def_property_readonly( "single_arc_settings",
                                &tp::MultiArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getInitialStateList,
                                R"doc(
            **read-only**

            List of single arc settings (e.g. entry j of this list is the single-arc propagator setting for arc j).

            :type: list[SingleArcPropagatorSettings]

)doc" );





    py::class_< tp::HybridArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::HybridArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::PropagatorSettings< STATE_SCALAR_TYPE > >( m,
                                                               "HybridArcPropagatorSettings",
                                                               R"doc(

            Class derived from :class:`PropagatorSettings` to define settings for hybrid-arc dynamics

            Class derived from :class:`PropagatorSettings` to define settings for hybrid-arc dynamics
            An object of this type is typically created using the dedicated creation function :func:`~hybrid_arc`

      )doc" )
            .def_property_readonly( "processing_settings",
                                    &tp::HybridArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getOutputSettings,
                                    R"doc(
            **read-only**

            Settings that determine how the hybrid-arc propagation results are processed (e.g. if the results are used to update the body ephemeris, if
            data from each epoch is saved, etc.), and which data is printed to the console
            See `user guide <https://docs.tudat.space/en/latest/user-guide/state-propagation/propagation-setup/printing-processing-results.html#automatic-processing>`_ for more details.

            :type: HybridArcPropagatorProcessingSettings

.)doc" );

    py::class_< tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::PropagatorSettings< STATE_SCALAR_TYPE > >( m,
                                                               "SingleArcPropagatorSettings",
                                                               R"doc(

            Class derived from :class:`PropagatorSettings` to define settings for single-arc dynamics (of any type, including translational, rotational, etc.)
            An object of this type is typically created using the specific propagator settings creation function, such as :func:`~translational`,
            :func:`~rotational` of :func:`~multitype`

      )doc" )
            .def_property( "termination_settings",
                           &tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getTerminationSettings,
                           &tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::resetTerminationSettings,
                           R"doc(
            Settings for creating the object that checks whether the propagation is finished (e.g. on a final time, a final dependent variable, etc.)

            :type: PropagationTerminationSettings
                           )doc" )
            .def_property( "integrator_settings",
                           &tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getIntegratorSettings,
                           &tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::setIntegratorSettings,
                           R"doc(

            Settings for creating the numerical integrator object that is used to solve the equations of motion

            :type: IntegratorSettings
                           )doc" )
            .def_property_readonly( "processing_settings",
                                    &tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getOutputSettings,
                                    R"doc(
            **read-only**

            Settings that determine how the propagation results are processed (e.g. if the results are used to update the body ephemeris, if
            data from each epoch is saved, etc.), and which data is printed to the console
            See `user guide <https://docs.tudat.space/en/latest/user-guide/state-propagation/propagation-setup/printing-processing-results.html#automatic-processing>`_ for more details.

            :type: SingleArcPropagatorProcessingSettings

                            )doc" )
            .def_property_readonly( "print_settings",
                                    &tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getPrintSettings,
                                    R"doc(
            **read-only**

            Settings that determine which information is printed to the console during the propagation. NOTE: this object can also be retrieved from the ``processing_settings`` attribute
            See `user guide <https://docs.tudat.space/en/latest/user-guide/state-propagation/propagation-setup/printing-processing-results.html#console-output>`_ for more details.

            :type: PropagationPrintSettings

                            )doc" );

    py::class_< tp::TranslationalStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::TranslationalStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                   "TranslationalStatePropagatorSettings",
                                                                                   R"doc(

         `SingleArcPropagatorSettings`-derived class to define settings for single-arc translational dynamics.

      )doc" )

            .def( "get_propagated_state_size",
                  &tp::TranslationalStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getPropagatedStateSize )
            .def( "reset_and_recreate_acceleration_models",
                  &tp::TranslationalStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::resetAccelerationModelsMap,
                  py::arg( "new_acceleration_settings" ),
                  py::arg( "bodies" ) );

    //            .def_property_readonly("acceleration_settings",
    //                                   &tp::TranslationalStatePropagatorSettings<STATE_SCALAR_TYPE>::getAccelerationSettingsMap);

    py::class_< tp::MultiTypePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::MultiTypePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                   "MultiTypePropagatorSettings",
                                                                                   R"doc(

         `SingleArcPropagatorSettings`-derived class to define settings for propagation of multiple quantities.





      )doc" )
            .def( "reset_initial_states",
                  &tp::MultiTypePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::resetInitialStates,
                  py::arg( "initial_states" ) )
            .def( "recreate_state_derivative_models",
                  &tp::MultiTypePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::resetIntegratedStateModels,
                  py::arg( "bodies" ) )
            .def( "single_type_settings",
                  &tp::MultiTypePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleTypePropagatorSettings,
                  py::arg( "state_type" ) )
            .def_property_readonly( "propagator_settings_per_type",
                                    &tp::MultiTypePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >::getPropagatorSettingsMap,
                                    R"doc(

         None

         :type: dict[IntegratedStateType, list[SingleArcPropagatorSettings]]
      )doc" );

    py::class_< tp::RotationalStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::RotationalStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                   "RotationalStatePropagatorSettings",
                                                                                   R"doc(

         `SingleArcPropagatorSettings`-derived class to define settings for single-arc rotational state propagation.





      )doc" );

    py::class_< tp::MassPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::MassPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >(
            m, "MassPropagatorSettings", R"doc(No propagator documentation found.)doc" );

    py::class_< tp::CustomStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::CustomStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >(
            m, "CustomStatePropagatorSettings", R"doc(No propagator documentation found.)doc" );

    ///////////////////////////////////////////////////////////////////////////////////////

    py::class_< tp::PropagationTerminationSettings, std::shared_ptr< tp::PropagationTerminationSettings > >(
            m,
            "PropagationTerminationSettings",
            R"doc(

         Functional base class to define termination settings for the propagation.





      )doc" );

    py::class_< tp::PropagationDependentVariableTerminationSettings,
                std::shared_ptr< tp::PropagationDependentVariableTerminationSettings >,
                tp::PropagationTerminationSettings >( m,
                                                      "PropagationDependentVariableTerminationSettings",
                                                      R"doc(

         `PropagationTerminationSettings`-derived class to define termination settings for the propagation from dependent variables.





      )doc" );

    py::class_< tp::PropagationTimeTerminationSettings,
                std::shared_ptr< tp::PropagationTimeTerminationSettings >,
                tp::PropagationTerminationSettings >( m,
                                                      "PropagationTimeTerminationSettings",
                                                      R"doc(

         `PropagationTerminationSettings`-derived class to define termination settings for the propagation from propagation time.





      )doc" );

    py::class_< tp::PropagationCPUTimeTerminationSettings,
                std::shared_ptr< tp::PropagationCPUTimeTerminationSettings >,
                tp::PropagationTerminationSettings >( m,
                                                      "PropagationCPUTimeTerminationSettings",
                                                      R"doc(

         `PropagationTerminationSettings`-derived class to define termination settings for the propagation from CPU time.





      )doc" );

    py::class_< tp::PropagationCustomTerminationSettings,
                std::shared_ptr< tp::PropagationCustomTerminationSettings >,
                tp::PropagationTerminationSettings >( m,
                                                      "PropagationCustomTerminationSettings",
                                                      R"doc(

         `PropagationTerminationSettings`-derived class to define custom termination settings for the propagation.





      )doc" );

    py::class_< tp::PropagationHybridTerminationSettings,
                std::shared_ptr< tp::PropagationHybridTerminationSettings >,
                tp::PropagationTerminationSettings >( m,
                                                      "PropagationHybridTerminationSettings",
                                                      R"doc(

         `PropagationTerminationSettings`-derived class to define hybrid termination settings for the propagation.





      )doc" );

    py::class_< tp::NonSequentialPropagationTerminationSettings,
                std::shared_ptr< tp::NonSequentialPropagationTerminationSettings >,
                tp::PropagationTerminationSettings >(
            m, "NonSequentialPropagationTerminationSettings", R"doc(No propagator documentation found.)doc" );

    ///////////////////////////////////////////////////////////////////////////////////////

    m.def( "translational",
           py::overload_cast< const std::vector< std::string > &,
                              const tba::AccelerationMap &,
                              const std::vector< std::string > &,
                              const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > &,
                              const std::shared_ptr< tp::PropagationTerminationSettings >,
                              const tp::TranslationalPropagatorType,
                              const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > > &,
                              const double >( &tp::translationalStatePropagatorSettingsDeprecated< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "central_bodies" ),
           py::arg( "acceleration_models" ),
           py::arg( "bodies_to_integrate" ),
           py::arg( "initial_states" ),
           py::arg( "termination_settings" ),
           py::arg( "propagator" ) = tp::cowell,
           py::arg( "output_variables" ) = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
           py::arg( "print_interval" ) = TUDAT_NAN );

    m.def( "translational",
           tp::translationalStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "central_bodies" ),
           py::arg( "acceleration_models" ),
           py::arg( "bodies_to_integrate" ),
           py::arg( "initial_states" ),
           py::arg( "initial_time" ),
           py::arg( "integrator_settings" ),
           py::arg( "termination_settings" ),
           py::arg( "propagator" ) = tp::cowell,
           py::arg( "output_variables" ) = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
           py::arg( "processing_settings" ) = nullptr,
           R"doc(

Function to create translational state propagator settings with stopping condition at given final time.

Function to create translational state propagator settings for N bodies.
The propagated state vector is defined by the combination of integrated bodies, and their central body, the combination
of which define the relative translational states for which a differential equation is to be solved. The propagator
input defines the formulation in which the differential equations are set up
The dynamical models are defined by an ``AccelerationMap`` (dict[str, list[AccelerationModel]]), as created by :func:`~tudatpy.dynamics.propagation_setup.create_acceleration_models` function.
Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational.html>`__.

Parameters
----------
central_bodies : list[str]
    List of central bodies with respect to which the bodies to be integrated are propagated.
acceleration_models : dict[str, list[AccelerationModel]]
    Set of accelerations acting on the bodies to propagate, provided as acceleration models.
bodies_to_integrate : list[str]
    List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
initial_states : numpy.ndarray
    Initial states of the bodies to integrate (one initial state for each body, concatenated into a single array), provided in the same order as the bodies to integrate. The initial states must be expressed in Cartesian elements, w.r.t. the central body of each integrated body. The states must be defined with the same frame orientation as the global frame orientation of the environment (specified when creating a system of bodies, see for instance :func:`~tudatpy.dynamics.environment_setup.get_default_body_settings` and :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies`). Consequently, for N integrated bodies, this input is a vector with size size 6N.
initial_time : float
    Initial epoch of the numerical propagation
integrator_settings : IntegratorSettings
    Settings defining the numerical integrator that is to be used for the propagation

    .. note:: The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time

termination_settings : PropagationTerminationSettings
    Generic termination settings object to check whether the propagation should be ended.
propagator : TranslationalPropagatorType, default=cowell
    Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
output_variables : list[SingleDependentVariableSaveSettings], default=[]
    Object to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
processing_settings: SingleArcPropagatorProcessingSettings, default=[]
    Object to define how the numerical results are to be processed after the propagation terminates, and which information to print to the console during the propagation. See `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/printing_processing_results.html>`__ for details on all options.
    If this object is left empty default settings of the :class:`~SingleArcPropagatorProcessingSettings` class are used.

Returns
-------
TranslationalStatePropagatorSettings
    Translational state propagator settings object.






     )doc" );

    m.def( "rotational",
           py::overload_cast< const tba::TorqueModelMap &,
                              const std::vector< std::string > &,
                              const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > &,
                              const std::shared_ptr< tp::PropagationTerminationSettings >,
                              const tp::RotationalPropagatorType,
                              const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
                              const double >( &tp::rotationalStatePropagatorSettingsDeprecated< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "torque_models" ),
           py::arg( "bodies_to_integrate" ),
           py::arg( "initial_states" ),
           py::arg( "termination_settings" ),
           py::arg( "propagator" ) = tp::quaternions,
           py::arg( "output_variables" ) = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
           py::arg( "print_interval" ) = TUDAT_NAN );

    m.def( "rotational",
           &tp::rotationalStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "torque_models" ),
           py::arg( "bodies_to_integrate" ),
           py::arg( "initial_states" ),
           py::arg( "initial_time" ),
           py::arg( "integrator_settings" ),
           py::arg( "termination_settings" ),
           py::arg( "propagator" ) = tp::quaternions,
           py::arg( "output_variables" ) = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
           py::arg( "processing_settings" ) = nullptr,
           R"doc(

Function to create rotational state propagator settings.

Function to create rotational state propagator settings for N bodies.
The propagated state vector is defined by the integrated bodies, which defines the bodies for which the
differential equation defining the evolution of the rotational state between an
inertial and body-fixed frame are to be solved. The propagator input defines the
formulation in which the differential equations are set up. The dynamical models are
defined by a ``TorqueModelMap``, as created by :func:`~tudatpy.dynamics.propagation_setup.create_torque_models` function.
Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/rotational.html>`__


Parameters
----------
torque_models : TorqueModelMap
    Set of torques acting on the bodies to propagate, provided as torque models.
bodies_to_integrate : list[str]
    List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
initial_states : numpy.ndarray
    Initial rotational states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
    Regardless of the propagator that is selected, the initial rotational state is always defined as four quaternion entries, and the angular velocity of the body,
    as defined in more detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html#definition-of-rotational-state>`__.
initial_time : float
    Initial epoch of the numerical propagation
integrator_settings : IntegratorSettings
    Settings defining the numerical integrator that is to be used for the propagation

    .. note:: 
    
        The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time
    
termination_settings : PropagationTerminationSettings
    Generic termination settings object to check whether the propagation should be ended.
propagator : RotationalPropagatorType, default=quaternions
    Type of rotational propagator to be used (see `RotationalPropagatorType` enum).
output_variables : list[SingleDependentVariableSaveSettings], default=[]
    Object to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
processing_settings: SingleArcPropagatorProcessingSettings, default=[]
    Object to define how the numerical results are to be processed after the propagation terminates, and which information to print to the console during the propagation. See `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/printing_processing_results.html>`__ for details on all options.
    If this object is left empty default settings of the :class:`~SingleArcPropagatorProcessingSettings` class are used.

Returns
-------
RotationalStatePropagatorSettings
    Rotational state propagator settings object.






     )doc" );

    m.def( "mass",
           py::overload_cast< const std::vector< std::string >,
                              const std::map< std::string, std::vector< std::shared_ptr< tba::MassRateModel > > > &,
                              const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > &,
                              const std::shared_ptr< tp::PropagationTerminationSettings >,
                              const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
                              const double >( &tp::massPropagatorSettingsDeprecated< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "bodies_with_mass_to_propagate" ),
           py::arg( "mass_rate_models" ),
           py::arg( "initial_body_masses" ),
           py::arg( "termination_settings" ),
           py::arg( "output_variables" ) = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
           py::arg( "print_interval" ) = TUDAT_NAN );

    m.def( "mass",
           &tp::massPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies_with_mass_to_propagate" ),
           py::arg( "mass_rate_models" ),
           py::arg( "initial_body_masses" ),
           py::arg( "initial_time" ),
           py::arg( "integrator_settings" ),
           py::arg( "termination_settings" ),
           py::arg( "output_variables" ) = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
           py::arg( "processing_settings" ) = nullptr,
           R"doc(

Function to create mass propagator settings

Function to create mass propagator settings
It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
each body. In this function, the dependent variables to save are provided
as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
through the termination settings object provided.
Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/mass.html>`__


Parameters
----------
bodies_with_mass_to_propagate : list[str]
    List of bodies whose mass should be numerically propagated.
mass_rate_models : SelectedMassRateModelMap
    Mass rates associated to each body, provided as a mass rate settings object.
initial_body_masses : numpy.ndarray
    Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
initial_time : float
    Initial epoch of the numerical propagation
integrator_settings : IntegratorSettings
    Settings defining the numerical integrator that is to be used for the propagation

    .. note:: The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time

termination_settings : PropagationTerminationSettings
    Generic termination settings object to check whether the propagation should be ended.
output_variables : list[SingleDependentVariableSaveSettings], default=[]
    Object to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
processing_settings: SingleArcPropagatorProcessingSettings, default=[]
    Object to define how the numerical results are to be processed after the propagation terminates, and which information to print to the console during the propagation. See `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/printing_processing_results.html>`__ for details on all options.
    If this object is left empty default settings of the :class:`~SingleArcPropagatorProcessingSettings` class are used.

Returns
-------
SingleArcPropagatorSettings
    Mass propagator settings object.






     )doc" );

    m.def( "custom_state",
           &tp::customStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "state_derivative_function" ),
           py::arg( "initial_state" ),
           py::arg( "initial_time" ),
           py::arg( "integrator_settings" ),
           py::arg( "termination_settings" ),
           py::arg( "output_variables" ) = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
           py::arg( "processing_settings" ) = nullptr,
           R"doc(

Function to create custom propagator settings.

Function to create custom propagator settings.
By using this propagator, the user can define their own differential equation to be solved, rather than have the differential equation be
defined by accelerations (as is the case for translational state) or torques (as is the case for rotational state). This permits the user
additional flexibility to define their own model by adding (for instance) the integration of co-states to the propagation without having to
implement the governing dynamics into Tudat. This propagator requires a function of the form :math:`\frac{d\mathbf{x}}{dt}=\mathbf{f}(t,\mathbf{x})`,
with :math:`t` the current time, :math:`\mathbf{x}` the current state, and :math:`\mathbf{f}` the state derivative function. This function can
depend on any quantities of the user's choosing, for details on how to link the properties of the environment to this function, see `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/custom_models.html>`__.

Parameters
----------
state_derivative_function : callable[[float, numpy.ndarray[numpy.float64[m, 1]]], numpy.ndarray[numpy.float64[m, 1]]])
    Function :math:`\mathbf{f}` (ser above) to compute the derivative of the current custom state
initial_state : numpy.ndarray
    Initial value of the propagated custom state
initial_time : float
    Initial epoch of the numerical propagation
integrator_settings : IntegratorSettings
    Settings defining the numerical integrator that is to be used for the propagation

    .. note:: The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time

termination_settings : PropagationTerminationSettings
    Generic termination settings object to check whether the propagation should be ended.
output_variables : list[SingleDependentVariableSaveSettings], default=[]
    Object to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
processing_settings: SingleArcPropagatorProcessingSettings, default=[]
    Object to define how the numerical results are to be processed after the propagation terminates, and which information to print to the console during the propagation. See `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/printing_processing_results.html>`__ for details on all options.
    If this object is left empty default settings of the :class:`~SingleArcPropagatorProcessingSettings` class are used.

Returns
-------
SingleArcPropagatorSettings
    Custom propagator settings object.
     )doc" );

    m.def( "multitype",
           py::overload_cast< const std::vector< std::shared_ptr< tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > > >,
                              const std::shared_ptr< tp::PropagationTerminationSettings >,
                              const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
                              const double >( &tp::multiTypePropagatorSettingsDeprecated< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "propagator_settings_list" ),
           py::arg( "termination_settings" ),
           py::arg( "output_variables" ) = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
           py::arg( "print_interval" ) = TUDAT_NAN );

    m.def( "multitype",
           &tp::multiTypePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "propagator_settings_list" ),
           py::arg( "integrator_settings" ),
           py::arg( "initial_time" ),
           py::arg( "termination_settings" ),
           py::arg( "output_variables" ) = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
           py::arg( "processing_settings" ) = nullptr,
           R"doc(

Function to create multitype propagator settings.

Function to create multitype propagator settings.
It works by providing a list of SingleArcPropagatorSettings objects. When using this function,
only the termination and output settings provided here are used, any such settings in the
constituent propagator settings are ignored
Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/multi_type.html>`__

.. note:: The propagated state contains the state types in the following order: Translational ( **C** ), Rotational ( **R** ), Mass ( **M** ), and Custom ( **C** ).
        When propagating two bodies, an example of what the output state would look like is for instance:
        [ **T** Body 1, **T** Body 2, **R** Body 1, **R** Body 2, **M** Body 1, **M** Body 2 ]


Parameters
----------
propagator_settings_list : list[SingleArcPropagatorSettings]
    List of SingleArcPropagatorSettings objects to use.
integrator_settings : IntegratorSettings
    Settings defining the numerical integrator that is to be used for the propagation

    .. note:: The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time

initial_time : float
    Initial epoch of the numerical propagation
termination_settings : PropagationTerminationSettings
    Generic termination settings object to check whether the propagation should be ended.
output_variables : list[SingleDependentVariableSaveSettings], default=[]
    Object to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
processing_settings: SingleArcPropagatorProcessingSettings, default=[]
    Object to define how the numerical results are to be processed after the propagation terminates, and which information to print to the console during the propagation. See `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/printing_processing_results.html>`__ for details on all options.
    If this object is left empty default settings of the :class:`~SingleArcPropagatorProcessingSettings` class are used.

Returns
-------
MultiTypePropagatorSettings
    Multi-type propagator settings object.

     )doc" );

    m.def( "multi_arc_processing_settings", &multiArcProcessingSettings );

    m.def( "multi_arc",
           &tp::multiArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "single_arc_settings" ),
           py::arg( "transfer_state_to_next_arc" ) = false,
           py::arg( "processing_settings" ) = nullptr,
           R"doc(

Function to create multi-arc propagator settings.

Function to create multi-arc propagator settings. It works by providing separate settings for
each arc in a list.


Parameters
----------
single_arc_settings : list[SingleArcPropagatorSettings]
    List of SingleArcPropagatorSettings objects to use, one for each arc.
transfer_state_to_next_arc : bool, default=False
    Object to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
Returns
-------
MultiArcPropagatorSettings
    Multi-arc propagator settings object.






     )doc" );

    m.def( "hybrid_arc",
           &tp::hybridArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "single_arc_settings" ),
           py::arg( "multi_arc_settings" ),
           py::arg( "processing_settings" ) = nullptr,
           R"doc(

Function to create hybrid-arc propagator settings.

Function to create hybrid-arc propagator settings (i.e., a combination of single- and multi-arc dynamics).


Parameters
----------
single_arc_settings : SingleArcPropagatorSettings
    SingleArcPropagatorSettings object to use for the propagation.
multi_arc_settings : MultiArcPropagatorSettings
    Object to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
Returns
-------
HybridArcPropagatorSettings
    Hybrid-arc propagator settings object.






     )doc" );

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

    m.def( "time_termination",
           &tp::propagationTimeTerminationSettings,
           py::arg( "termination_time" ),
           py::arg( "terminate_exactly_on_final_condition" ) = false,
           R"doc(

 Function to create time termination settings for the propagation.

 Function to create time termination settings for the propagation.
 The propagation is stopped when the final time provided is reached. Note that the termination time is set as the
 absolute time (in seconds since J2000), not the time since the start of the propagation.
 Depending on the sign of the time step of the numerical integrator, the termination time will be treated as an
 upper bound (for positive time step) or lower bound (for negative time step).
 The simulator will normally finish the final time-step, which may cause the termination time to be slightly exceeded.
 This behaviour can be suppressed by providing the optional input argument
 ``terminate_exactly_on_final_condition=True``, in which case the final propagation step will be *exactly* on the
 specified time.

 Parameters
 ----------
 termination_time : float
     Final time of the propagation.
 terminate_exactly_on_final_condition : bool, default=False
     Denotes whether the propagation is to terminate exactly on the final condition, or whether it is to terminate on the first step where it is violated.
 Returns
 -------
 PropagationTimeTerminationSettings
     Time termination settings object.



 Notes
 -----
 To reach *exactly* the final time, state derivative function evaluations beyond the final
 time may be required by the propagator. Reaching the final condition exactly is an iterative process and
 very minor deviations from the specified final condition can occur.



 Examples
 --------
 In this example, we set the termination time of the propagation equal to one day (86400 s).

 .. code-block:: python

   # Set termination time (in seconds since J2000)
   termination_time = simulation_start_epoch + 86400.0
   # Create time termination settings
   termination_settings = propagation_setup.propagator.time_termination( termination_time )


     )doc" );

    m.def( "cpu_time_termination",
           &tp::propagationCPUTimeTerminationSettings,
           py::arg( "cpu_termination_time" ),
           R"doc(

 Function to create CPU time termination settings for the propagation.

 Function to create CPU time termination settings for the propagation.
 The propagation is stopped when the final CPU time provided is reached.


 Parameters
 ----------
 cpu_termination_time : float
     Maximum CPU time for the propagation.
 Returns
 -------
 PropagationCPUTimeTerminationSettings
     CPU time termination settings object.





 Examples
 --------
 In this case, we set a CPU time termination setting so that the propagation stops once your computer has run it
 for 120 seconds.

 .. code-block:: python

   # Set CPU time to 120 seconds
   cpu_termination_time = 120.0
   # Create termination settings
   termination_settings = propagation_setup.propagator.cpu_time_termination( cpu_termination_time )


     )doc" );

    m.def( "dependent_variable_termination",
           &tp::propagationDependentVariableTerminationSettings,
           py::arg( "dependent_variable_settings" ),
           py::arg( "limit_value" ),
           py::arg( "use_as_lower_limit" ),
           py::arg( "terminate_exactly_on_final_condition" ) = false,
           py::arg( "termination_root_finder_settings" ) = nullptr,
           R"doc(

 Function to create termination settings for the propagation based on a dependent variable.

 Function to create termination settings for the propagation based on the value of a dependent variable.
 The propagation is stopped when a provided upper or lower limit value is reached.
 The simulator will normally finish the final time-step, which may cause the dependent variable to be slightly exceeded.
 This behaviour can be suppressed by providing the optional input argument
 ``terminate_exactly_on_final_condition=True``, in which case the final propagation step will be *exactly* on the
 specified dependent variable value.


 Parameters
 ----------
 dependent_variable_settings : SingleDependentVariableSaveSettings
     Dependent variable object to be used as termination setting.
 limit_value : float
     Limit value of the dependent variable; if reached, the propagation is stopped.
 use_as_lower_limit : bool, default=False
     Denotes whether the limit value should be used as lower or upper limit.
 terminate_exactly_on_final_condition : bool, default=False
     Denotes whether the propagation is to terminate exactly on the final condition, or whether it is to terminate on the first step where it is violated.
 termination_root_finder_settings : RootFinderSettings, default=None
     Settings object to create root finder used to converge on exact final condition.
 Returns
 -------
 PropagationDependentVariableTerminationSettings
     Dependent variable termination settings object.



 Notes
 -----
 To reach *exactly* the final dependent variable value, state derivative function evaluations beyond the final
 time may be required by the propagator. Reaching the final condition exactly is an iterative process and
 very minor deviations from the specified final condition can occur.



 Examples
 --------
 Below, an example is shown for termination on a given vehicle altitude. The exact termination condition is defined
 in the ``termination_settings``. The propagation is terminated once the *lower* limit of 25 km in altitude is
 reached (as the ``use_as_lower_limit`` is set to ``True``). To use the above settings to terminate when an
 *upper* limit of 25 km is reached, set this boolean to ``False``. In this example, we also want to stop exactly
 at 25 km, so we set ``terminate_exactly_on_final_condition`` to ``True``, and we specify ``termination_root_finder_settings``.

 .. code-block:: python

     # Set dependent variable to be checked as termination setting
     termination_variable = propagation_setup.dependent_variable.altitude(
         "Spacecraft", "Earth"
     )
     # Create termination settings
     termination_settings = propagation_setup.propagator.dependent_variable_termination(
         dependent_variable_settings=termination_variable,
         limit_value=25.0e3,
         use_as_lower_limit=True,
         terminate_exactly_on_final_condition=True,
         termination_root_finder_settings=root_finders.secant(
             maximum_iteration=5,
             maximum_iteration_handling=root_finders.MaximumIterationHandling.accept_result,
         ),
     )


     )doc" );

    m.def( "custom_termination",
           &tp::popagationCustomTerminationSettings,
           py::arg( "custom_condition" ),
           R"doc(

 Function to create custom termination settings for the propagation.

 Function to create custom termination settings for the propagation.
 The propagation is stopped when the condition provided is verified.
 This custom function should take the current time as input and output a Boolean. It can use internal variables
 and calculations, for example retrieved from the environment.


 Parameters
 ----------
 custom_condition : callable[[float], bool]
     Function of time (independent variable) which is called during the propagation and returns a boolean value denoting whether the termination condition is verified.
 Returns
 -------
 PropagationCustomTerminationSettings
     Custom termination settings object.





 Examples
 --------

 .. code-block:: python

   # Create custom function returning a bool
   def custom_termination_function(time: float):
       # Do something
       set_condition = ...
       # Return bool
       return set_condition

   # Create termination settings
   termination_settings = propagation_setup.propagator.custom_termination(
     custom_termination_function)


     )doc" );

    m.def( "custom_termination_with_state_input",
           &tp::popagationCustomTerminationSettingsFromFullState,
           py::arg( "custom_condition" ),
           R"doc(No propagator documentation found.)doc" );

    m.def( "hybrid_termination",
           &tp::propagationHybridTerminationSettings,
           py::arg( "termination_settings" ),
           py::arg( "fulfill_single_condition" ),
           R"doc(

 Function to create hybrid termination settings for the propagation.

 Function to create hybrid termination settings for the propagation. This function can be used
 to define that all conditions or a single condition of the conditions provided must be met to
 stop the propagation. Each termination condition should be created according to each individual function
 and then added to a list of termination conditions.

 .. note::

     When using this option, the :attr:`~tudatpy.dynamics.propagation.SingleArcSimulationResults.termination_details` of
     the simulation results object (obtained from here after a propagation: :attr:`~tudatpy.dynamics.simulator.SingleArcSimulator.propagation_results`)
     is of derived type :class:`~tudatpy.dynamics.propagation.PropagationTerminationDetailsFromHybridCondition`.

     See the :attr:`~tudatpy.dynamics.propagation.PropagationTerminationDetailsFromHybridCondition.was_condition_met_when_stopping` attribute for an example of how to retrieve which condition was met when the propagation was terminated.


 Parameters
 ----------
 termination_settings : list[PropagationTerminationSettings]
     List of single PropagationTerminationSettings objects to be checked during the propagation.
 fulfill_single_condition : bool
     Whether only a single condition of those provided must be met to stop the propagation (true) or all of them simultaneously (false).
 Returns
 -------
 PropagationHybridTerminationSettings
     Hybrid termination settings object.





 Examples
 --------
 In the following example, the propagation will terminate once *one of the three* termination settings (simulated time, cpu time, altitude)
 has reached the imposed limit value. The ``fulfill_single_condition`` variable determines whether the propagation
 terminates once a *single* condition is met (if True, as above) or once *all* conditions are met (False).

 .. code-block:: python

     # Set simulation termination time
     termination_time = simulation_start_epoch + 86400.0
     # Create simulation time termination setting
     time_termination_settings = propagation_setup.propagator.time_termination(
         termination_time
     )

     # Set dependent variable termination setting
     termination_variable = propagation_setup.dependent_variable.altitude(
         "Spacecraft", "Earth"
     )
     # Create altitude-based termination setting
     altitude_termination_settings = (
         propagation_setup.propagator.dependent_variable_termination(
             dependent_variable_settings=termination_variable,
             limit_value=25.0e3,
             use_as_lower_limit=True,
         )
     )

     # Set cpu termination time
     cpu_termination_time = 120.0
     # Create cpu time termination setting
     cpu_termination_settings = propagation_setup.propagator.cpu_time_termination(
         cpu_termination_time
     )

     # Store termination setting objects in a list
     termination_settings_list = [
         time_termination_settings,
         altitude_termination_settings,
         cpu_termination_settings,
     ]

     # Create hybrid termination settings
     termination_settings = propagation_setup.propagator.hybrid_termination(
         termination_settings_list, fulfill_single_condition=True
     )


     )doc" );

    m.def( "non_sequential_termination",
           &tp::nonSequentialPropagationTerminationSettings,
           py::arg( "forward_termination_settings" ),
           py::arg( "backward_termination_settings" ),
           R"doc(

 Function to create non-sequential termination settings for the propagation.

 Function to create non-sequential termination settings for the propagation. By using this setting,
 the propagation of the dynamics along an arc is propagated starting from some point (initial time and state) along the arc, and then
 propagating both forwards and backwards in time. This termination condition allows the user to specify
 termination conditions for the propagations forwards and backwards in time. These two propagations are then
 internally performed separately, but the resulting propagation results provide the concatenated results
 from the two and effectively constitute
 By using this function, the propagation controlled by the ``forward_termination_settings`` automatically has a positive
 time step, while the ``backward_termination_settings`` automatically has a negative time step. By definition,
 both will start from the same initial time and state, which are provided in the propagation settings


 Parameters
 ----------
 forward_termination_settings : PropagationTerminationSettings
     Propagation termination setting for the forward-in-time propagation
 backward_termination_settings : PropagationTerminationSettings
     Propagation termination setting for the backwards-in-time propagation
 Returns
 -------
 PropagationTerminationSettings
     Termination settings object for forward- and backwards-in-time propagation.






     )doc" );

    m.def( "add_dependent_variable_settings",
           &tp::addDepedentVariableSettings< double >,
           py::arg( "dependent_variable_settings" ),
           py::arg( "propagator_settings" ),
           R"doc(

 Function to add dependent variables to existing propagator settings.

 Function to add dependent variables to existing :class:`~tudatpy.dynamics.propagation_setup.propagator.SingleArcPropagatorSettings`
 object. This function is added as an alternative to teh regular manner in which to defined dependent variables (use of input to
 functions for single-arc propagator settings :func:`~tudatpy.dynamics.propagation_setup.propagator.translational`,
 :func:`~tudatpy.dynamics.propagation_setup.propagator.rotational`, :func:`~tudatpy.dynamics.propagation_setup.propagator.mass`,
 :func:`~tudatpy.dynamics.propagation_setup.propagator.multitype`). Typically, this function is used to modify
 existing propagator settings in a loop when running multiple simulations


 Parameters
 ----------
 dependent_variable_settings : List[ SingleDependentVariableSaveSettings ]
     List of dependent variable settings that are to be added to propagator settings. Note that this function adds settings, and does not replace any existing settings (nor does it check for duplicate settings).
 propagator_settings : SingleArcPropagatorSettings
     Propagator settings to which the additional dependent variables settings are to be added.
 Returns
 -------
 None
     None






     )doc" );

     m.def( "get_integrated_type_and_body_list",
           &tp::getIntegratedTypeAndBodyList< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "propagator_settings" ) );

    m.def( "get_single_integration_size", &tp::getSingleIntegrationSize, py::arg( "state_type" ) );
}

}  // namespace propagator
}  // namespace propagation_setup
}  // namespace dynamics
}  // namespace tudatpy
