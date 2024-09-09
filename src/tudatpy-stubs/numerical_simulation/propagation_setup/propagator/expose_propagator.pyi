import numpy
import typing
__all__ = ['CustomStatePropagatorSettings', 'HybridArcPropagatorProcessingSettings', 'HybridArcPropagatorSettings', 'MassPropagatorSettings', 'MultiArcPropagatorProcessingSettings', 'MultiArcPropagatorSettings', 'MultiTypePropagatorSettings', 'NonSequentialPropagationTerminationSettings', 'PropagationCPUTimeTerminationSettings', 'PropagationCustomTerminationSettings', 'PropagationDependentVariableTerminationSettings', 'PropagationHybridTerminationSettings', 'PropagationPrintSettings', 'PropagationTerminationSettings', 'PropagationTerminationTypes', 'PropagationTimeTerminationSettings', 'PropagatorProcessingSettings', 'PropagatorSettings', 'RotationalPropagatorType', 'RotationalStatePropagatorSettings', 'SingleArcPropagatorProcessingSettings', 'SingleArcPropagatorSettings', 'StateType', 'TranslationalPropagatorType', 'TranslationalStatePropagatorSettings', 'add_dependent_variable_settings', 'cowell', 'cpu_time_stopping_condition_type', 'cpu_time_termination', 'custom_state', 'custom_stopping_condition_type', 'custom_termination', 'custom_type', 'dependent_variable_stopping_condition_type', 'dependent_variable_termination', 'encke', 'exponential_map', 'gauss_keplerian', 'gauss_modified_equinoctial', 'hybrid_arc', 'hybrid_stopping_condition_type', 'hybrid_termination', 'hybrid_type', 'mass', 'mass_type', 'modified_rodrigues_parameters', 'multi_arc', 'multitype', 'non_sequential_termination', 'quaternions', 'rotational', 'rotational_type', 'time_stopping_condition_type', 'time_termination', 'translational', 'translational_type', 'undefined_rotational_propagator', 'undefined_translational_propagator', 'unified_state_model_exponential_map', 'unified_state_model_modified_rodrigues_parameters', 'unified_state_model_quaternions']

class CustomStatePropagatorSettings(SingleArcPropagatorSettings):
    """
		"""

class HybridArcPropagatorProcessingSettings(PropagatorProcessingSettings):
    """Class to define settings on how the numerical results are to be used for hybrid-arc propagations
	
	Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
	derived from :class:`PropagatorProcessingSettings`.
	Instances of this class are typically not created by the user. A settings object is
	instantiated through the factory function :func:`~hybrid_arc` to define hybrid-arc propagator settings.
	This object contains a :class:`SingleArcPropagatorProcessingSettings` object and a :class:`MultuArcPropagatorProcessingSettings` ,
	containing the processing settings for the constituents of the hybrid-arc propagatioon.
	"""
    clear_numerical_solution: bool
    set_integrated_result: bool

    def set_print_settings_for_all_arcs(self, print_settings: PropagationPrintSettings) -> None:
        """
        Function that sets the same print settings for each arc in the multi-arc propagation.
        
        
        	:param single_arc_print_settings:
        		Propagation print settings that are applied to each constituent single-arc settings, overriding any existing settings.
        """

    @property
    def multi_arc_settings(self) -> MultiArcPropagatorProcessingSettings:
        """
        Processing settings for the single-arc component of the multi-arc propagation.
        """

    @property
    def single_arc_settings(self) -> SingleArcPropagatorProcessingSettings:
        """
        Processing settings for the single-arc component of the hybrid-arc propagation.
        """

class HybridArcPropagatorSettings(PropagatorSettings):
    """`PropagatorSettings`-derived class to define settings for hybrid-arc dynamics.
	"""

    @property
    def processing_settings(self) -> HybridArcPropagatorProcessingSettings:
        ...

class MassPropagatorSettings(SingleArcPropagatorSettings):
    """
		"""

class MultiArcPropagatorProcessingSettings(PropagatorProcessingSettings):
    """Class to define settings on how the numerical results are to be used for multi-arc propagations
	
	Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
	derived from :class:`PropagatorProcessingSettings`.
	Instances of this class are typically not created by the user. A settings object is
	instantiated through the factory function :func:`~multi_arc` to define multi-arc propagator settings.
	This object contains a list of :class:`SingleArcPropagatorProcessingSettings` objects, containing the processing settings for each constituent arc.
	"""

    def set_print_settings_for_all_arcs(self, single_arc_print_settings: PropagationPrintSettings) -> None:
        """
        Function that sets the same print settings for each arc in the multi-arc propagation.
        
        
        	:param single_arc_print_settings:
        		Propagation print settings that are applied to each constituent single-arc settings, overriding any existing settings.
        """

    @property
    def identical_settings_per_arc(self) -> bool:
        ...

    @property
    def print_output_on_first_arc_only(self) -> bool:
        """
        Variable defining whether the ``set_print_settings_for_all_arcs`` function has been used to define identical print settings for each arc.
        """

    @print_output_on_first_arc_only.setter
    def print_output_on_first_arc_only(self, arg1: bool) -> None:
        ...

    @property
    def single_arc_settings(self) -> list[SingleArcPropagatorProcessingSettings]:
        """
        List containing the processing settings for each constituent arc
        """

class MultiArcPropagatorSettings(PropagatorSettings):
    """`PropagatorSettings`-derived class to define settings for multi-arc dynamics.
	"""

    @property
    def processing_settings(self) -> MultiArcPropagatorProcessingSettings:
        ...

class MultiTypePropagatorSettings(SingleArcPropagatorSettings):
    """`SingleArcPropagatorSettings`-derived class to define settings for propagation of multiple quantities.
	"""

    def recreate_state_derivative_models(self, bodies: ...) -> None:
        ...

    def reset_initial_states(self, initial_states: numpy.ndarray) -> None:
        ...

    def single_type_settings(self, state_type: StateType) -> SingleArcPropagatorSettings:
        ...

    @property
    def propagator_settings_per_type(self) -> dict[StateType, list[SingleArcPropagatorSettings]]:
        ...

class NonSequentialPropagationTerminationSettings(PropagationTerminationSettings):
    """
		"""

class PropagationCPUTimeTerminationSettings(PropagationTerminationSettings):
    """`PropagationTerminationSettings`-derived class to define termination settings for the propagation from CPU time.
	"""

class PropagationCustomTerminationSettings(PropagationTerminationSettings):
    """`PropagationTerminationSettings`-derived class to define custom termination settings for the propagation.
	"""

class PropagationDependentVariableTerminationSettings(PropagationTerminationSettings):
    """`PropagationTerminationSettings`-derived class to define termination settings for the propagation from dependent variables.
	"""

class PropagationHybridTerminationSettings(PropagationTerminationSettings):
    """`PropagationTerminationSettings`-derived class to define hybrid termination settings for the propagation.
	"""

class PropagationPrintSettings:
    """Class to save settings on what is to be written to the console during the propagation of a single arc.
	"""

    def disable_all_printing(self) -> None:
        """
        Function enabling all printing (e.g. sets all boolean attributes to False, and disables all other output as well)
        """

    def enable_all_boolean_printing(self) -> None:
        """
        Function enabling all True/False printing (e.g. sets all boolean attributes to True)
        """

    def enable_all_printing(self, results_print_frequency_in_seconds: float, results_print_frequency_in_steps: int) -> None:
        """
        Function enabling all True/False printing (e.g. sets all boolean attributes to True), and setting the non-boolean
        attributes to values defined here.
        
        
        
        	:param results_print_frequency_in_seconds:
        		See ``results_print_frequency_in_seconds`` class attribute
        	:param results_print_frequency_in_steps:
        		See ``results_print_frequency_in_steps`` class attribute
        """

    @property
    def print_dependent_variable_indices(self) -> bool:
        """
        Boolean defining whether the meaning and indices of the
        entries of the dependent variable data are to be printed to
        the console (before the propagation).
        
        .. note:: The same information can be retrieved from the
                  :py:attr:`SingleArcSimulationResults.dependent_variable_ids`
                  attribute.
        """

    @print_dependent_variable_indices.setter
    def print_dependent_variable_indices(self, arg1: bool) -> None:
        ...

    @property
    def print_dependent_variables_during_propagation(self) -> bool:
        """
        Boolean defining whether the dependent variables are to be printed during the propagation along with the state,
        at steps/epochs define by the ``results_print_frequency_in_seconds`` and/or ``results_print_frequency_in_steps`` inputs.
        """

    @print_dependent_variables_during_propagation.setter
    def print_dependent_variables_during_propagation(self, arg1: bool) -> None:
        ...

    @property
    def print_initial_and_final_conditions(self) -> bool:
        """
        Boolean defining whether the initial and final conditions (state and time)
        are to be printed to the console (beforee and after propagation, respectively).
        """

    @print_initial_and_final_conditions.setter
    def print_initial_and_final_conditions(self, arg1: bool) -> None:
        ...

    @property
    def print_number_of_function_evaluations(self) -> bool:
        """
        Boolean defining whether the number of function evaluations that
        were performed is to be printed to the console (after propagation).
        """

    @print_number_of_function_evaluations.setter
    def print_number_of_function_evaluations(self, arg1: bool) -> None:
        ...

    @property
    def print_processed_state_indices(self) -> bool:
        """
        Boolean defining whether the meaning and indices of the
        entries of the processed state vector are to be printed to
        the console (after the propagation). The distinction between the
        propagated and processed (or conventuional) state representation is described in
        detail `here <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/propagator_types.html>`_.
        Summarizing: the processed state is the 'typical' formulation of the state (for translational dynamics: Cartesian states).
        
        .. note:: The same information can be retrieved from the
                  :py:attr:`SingleArcSimulationResults.processed_state_ids`
                  attribute.
        """

    @print_processed_state_indices.setter
    def print_processed_state_indices(self, arg1: bool) -> None:
        ...

    @property
    def print_propagation_clock_time(self) -> bool:
        """
        Boolean defining whether the total clock time taken for the propagation
        is to be printed to the console (after propagation).
        """

    @print_propagation_clock_time.setter
    def print_propagation_clock_time(self, arg1: bool) -> None:
        ...

    @property
    def print_state_indices(self) -> bool:
        """
        Boolean defining whether the meaning and indices of the
        entries of the state vector are to be printed to
        the console (before the propagation).
        
        .. note:: The same information can be retrieved from the
                  :py:attr:`SingleArcSimulationResults.propagated_state_ids`
                  attribute.
        """

    @print_state_indices.setter
    def print_state_indices(self, arg1: bool) -> None:
        ...

    @property
    def print_termination_reason(self) -> bool:
        """
        Boolean defining whether the reason for propagation termination
        is to be printed to the console (after propagation).
        """

    @print_termination_reason.setter
    def print_termination_reason(self, arg1: bool) -> None:
        ...

    @property
    def results_print_frequency_in_seconds(self) -> float:
        """
        Variable indicating how often (in seconds of simulation time)
        the current state and time are to be printed to the console (by default, set to NaN - they are never printed).
        In case this setting is active (e.g. not NaN), and the ``results_print_frequency_in_steps`` setting is active,
        the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.
        """

    @results_print_frequency_in_seconds.setter
    def results_print_frequency_in_seconds(self, arg1: float) -> None:
        ...

    @property
    def results_print_frequency_in_steps(self) -> float:
        """
        Variable indicating how often (in number of full integration steps)
        the current state and time are to be printed to the console (by default, set to 0 - they are never printed).
        In case this setting is active (e.g. not 0), and the ``results_print_frequency_in_seconds`` setting is active,
        the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.
        """

    @results_print_frequency_in_steps.setter
    def results_print_frequency_in_steps(self, arg1: float) -> None:
        ...

class PropagationTerminationSettings:
    """Functional base class to define termination settings for the propagation.
	"""

class PropagationTerminationTypes:
    """Enumeration of possible propagation termination types
	
	
	:member time_stopping_condition:
	:member cpu_time_stopping_condition:
	:member dependent_variable_stopping_condition:
	:member hybrid_stopping_condition:
	:member custom_stopping_condition:
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, PropagationTerminationTypes]]
    cpu_time_stopping_condition_type: typing.ClassVar[PropagationTerminationTypes]
    custom_stopping_condition_type: typing.ClassVar[PropagationTerminationTypes]
    dependent_variable_stopping_condition_type: typing.ClassVar[PropagationTerminationTypes]
    hybrid_stopping_condition_type: typing.ClassVar[PropagationTerminationTypes]
    time_stopping_condition_type: typing.ClassVar[PropagationTerminationTypes]

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

class PropagationTimeTerminationSettings(PropagationTerminationSettings):
    """`PropagationTerminationSettings`-derived class to define termination settings for the propagation from propagation time.
	"""

class PropagatorProcessingSettings:
    """Base class to define settings on how the numerical results are to be used
	
	Base class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
	Instances of this class are typically not created by the user. Settings objects for derived class of single-, multi- and hybrid arc propagation are
	instantiated through the factory functions to define propagator settings (such as :func:`~translational` or :func:`~multi_arc`) in this module
	"""
    create_dependent_variable_interface: bool

    @property
    def clear_numerical_solution(self) -> bool:
        """
        Boolean defining whether the propagation results should be
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
        """

    @clear_numerical_solution.setter
    def clear_numerical_solution(self, arg1: bool) -> None:
        ...

    @property
    def set_integrated_result(self) -> bool:
        """
        Boolean defining whether the propagation results are to
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
        """

    @set_integrated_result.setter
    def set_integrated_result(self, arg1: bool) -> None:
        ...

class PropagatorSettings:
    """Functional base class to define settings for propagators.
	
	Base class to define settings for propagators. Derived classes are split into settings for single- and multi-arc dynamics.
	"""
    initial_states: numpy.ndarray

class RotationalPropagatorType:
    """Enumeration of available rotational propagator types.
	
	
	:member undefined_rotational_propagator:
	:member quaternions:
	:member modified_rodrigues_parameters:
	:member exponential_map:
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, RotationalPropagatorType]]
    exponential_map: typing.ClassVar[RotationalPropagatorType]
    modified_rodrigues_parameters: typing.ClassVar[RotationalPropagatorType]
    quaternions: typing.ClassVar[RotationalPropagatorType]
    undefined_rotational_propagator: typing.ClassVar[RotationalPropagatorType]

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

class RotationalStatePropagatorSettings(SingleArcPropagatorSettings):
    """`SingleArcPropagatorSettings`-derived class to define settings for single-arc rotational state propagation.
	"""

class SingleArcPropagatorProcessingSettings(PropagatorProcessingSettings):
    """Class to define settings on how the numerical results are to be used for single-arc propagations
	
	Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
	derived from :class:`PropagatorProcessingSettings`.
	Instances of this class are typically not created by the user. A settings object is
	instantiated through the factory functions to define single-arc propagator settings (such as :func:`~translational` or :func:`~rotational`) in this module
	"""

    @property
    def print_settings(self) -> PropagationPrintSettings:
        """
        Settings object defining which quantities should be printed to the console before, during and after the propagation. By default, this
        object is instantiated to print nothing.
        """

    @property
    def results_save_frequency_in_seconds(self) -> float:
        """
        Variable indicating how often (in seconds of simulation time)
        the propagated time, state, dependent variables, etc. are to be saved to data structures containing the results
        (by default, set to NaN - they are not saved based on propagation time; see below and ``results_save_frequency_in_steps`` attribute ).
        In case this setting is active, and set to :math:`\\Delta t`, the data are saved as soon as the current time step is :math:`\\ge \\Delta t` after the
        last step at which data was saved.
        In case this setting is active (e.g. not NaN), and the ``results_save_frequency_in_steps`` setting is active,
        the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.
        """

    @results_save_frequency_in_seconds.setter
    def results_save_frequency_in_seconds(self, arg1: float) -> None:
        ...

    @property
    def results_save_frequency_in_steps(self) -> int:
        """
        Variable indicating how often (in number of integrator steps)
        the propagated time, state, dependent variables, etc. are to be saved to data structures containing the results
        (by default, set to 1 - they are saved every time step). If this setting is set to 0, the data is never saved based on number of steps.
        In case this setting is active (e.g. not 0), and the ``results_save_frequency_in_seconds`` setting is active,
        the data is saved as soon as *one* of the two conditions (number of seconds, or number of steps) is met.
        """

    @results_save_frequency_in_steps.setter
    def results_save_frequency_in_steps(self, arg1: int) -> None:
        ...

class SingleArcPropagatorSettings(PropagatorSettings):
    """`PropagatorSettings`-derived class to define settings for single-arc dynamics.
	"""
    integrator_settings: ...
    termination_settings: ...

    @property
    def print_settings(self) -> PropagationPrintSettings:
        ...

    @property
    def processing_settings(self) -> SingleArcPropagatorProcessingSettings:
        ...

class StateType:
    """Enumeration of available integrated state types.
	
	
	:member hybrid_type:
	:member translational_type:
	:member rotational_type:
	:member body_mass_type:
	:member custom_type:
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, StateType]]
    custom_type: typing.ClassVar[StateType]
    hybrid_type: typing.ClassVar[StateType]
    mass_type: typing.ClassVar[StateType]
    rotational_type: typing.ClassVar[StateType]
    translational_type: typing.ClassVar[StateType]

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

class TranslationalPropagatorType:
    """Enumeration of available translational propagator types.
	
	
	:member undefined_translational_propagator:
	:member cowell:
	:member encke:
	:member gauss_keplerian:
	:member gauss_modified_equinoctial:
	:member unified_state_model_quaternions:
	:member unified_state_model_modified_rodrigues_parameters:
	:member unified_state_model_exponential_map:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, TranslationalPropagatorType]]
    cowell: typing.ClassVar[TranslationalPropagatorType]
    encke: typing.ClassVar[TranslationalPropagatorType]
    gauss_keplerian: typing.ClassVar[TranslationalPropagatorType]
    gauss_modified_equinoctial: typing.ClassVar[TranslationalPropagatorType]
    undefined_translational_propagator: typing.ClassVar[TranslationalPropagatorType]
    unified_state_model_exponential_map: typing.ClassVar[TranslationalPropagatorType]
    unified_state_model_modified_rodrigues_parameters: typing.ClassVar[TranslationalPropagatorType]
    unified_state_model_quaternions: typing.ClassVar[TranslationalPropagatorType]

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

class TranslationalStatePropagatorSettings(SingleArcPropagatorSettings):
    """`SingleArcPropagatorSettings`-derived class to define settings for single-arc translational dynamics.
	"""

    def get_propagated_state_size(self) -> int:
        ...

    def reset_and_recreate_acceleration_models(self, new_acceleration_settings: dict[str, dict[str, list[...]]], bodies: ...) -> None:
        ...

def add_dependent_variable_settings(dependent_variable_settings: list[...], propagator_settings: SingleArcPropagatorSettings) -> None:
    """Function to add dependent variables to existing propagator settings.
	
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
	"""

def cpu_time_termination(cpu_termination_time: float) -> PropagationTerminationSettings:
    """Factory function to create CPU time termination settings for the propagation.
	
	Factory function to create CPU time termination settings for the propagation.
	The propagation is stopped when the final CPU time provided is reached.
	
	
	:param cpu_termination_time:
			Maximum CPU time for the propagation.
	:return:
			CPU time termination settings object.
	"""

def custom_state(state_derivative_function: typing.Callable[[float, numpy.ndarray], numpy.ndarray], initial_state: numpy.ndarray, initial_time: float, integrator_settings: ..., termination_settings: ..., output_variables: list[...]=[], processing_settings: SingleArcPropagatorProcessingSettings=...) -> CustomStatePropagatorSettings:
    ...

def custom_termination(custom_condition: typing.Callable[[float], bool]) -> PropagationTerminationSettings:
    """Factory function to create custom termination settings for the propagation.
	
	Factory function to create custom termination settings for the propagation.
	The propagation is stopped when the condition provided is verified.
	This custom function should take the current time as input and output a Boolean. It can use internal variables
	and calculations, for example retrieved from the environment.
	
	
	:param custom_condition:
			Function of time (independent variable) which is called during the propagation and returns a boolean value denoting whether the termination condition is verified.
	:return:
			Custom termination settings object.
	"""

def dependent_variable_termination(dependent_variable_settings: ..., limit_value: float, use_as_lower_limit: bool, terminate_exactly_on_final_condition: bool=False, termination_root_finder_settings: ...=None) -> PropagationTerminationSettings:
    """Factory function to create termination settings for the propagation based on a dependent variable.
	
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
	"""

def hybrid_arc(single_arc_settings: SingleArcPropagatorSettings, multi_arc_settings: MultiArcPropagatorSettings, processing_settings: HybridArcPropagatorProcessingSettings=...) -> HybridArcPropagatorSettings:
    """Factory function to create hybrid-arc propagator settings.
	
	Factory function to create hybrid-arc propagator settings (i.e., a combination of single- and multi-arc dynamics).
	
	
	:param single_arc_settings:
			SingleArcPropagatorSettings object to use for the propagation.
	:param multi_arc_settings:
			MultiArcPropagatorSettings object to use for the propagation.
	:return:
			Hybrid-arc propagator settings object.
	"""

def hybrid_termination(termination_settings: list[PropagationTerminationSettings], fulfill_single_condition: bool) -> PropagationTerminationSettings:
    """Factory function to create hybrid termination settings for the propagation.
	
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
	"""

@typing.overload
def mass(bodies_with_mass_to_propagate: list[str], mass_rate_models: dict[str, list[...]], initial_body_masses: numpy.ndarray, termination_settings: ..., output_variables: list[...]=[], print_interval: float=...) -> MassPropagatorSettings:
    ...

@typing.overload
def mass(bodies_with_mass_to_propagate: list[str], mass_rate_models: dict[str, list[...]], initial_body_masses: numpy.ndarray, initial_time: float, integrator_settings: ..., termination_settings: ..., output_variables: list[...]=[], processing_settings: SingleArcPropagatorProcessingSettings=...) -> MassPropagatorSettings:
    """Factory function to create mass propagator settings
	
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
	"""

def multi_arc(single_arc_settings: list[SingleArcPropagatorSettings], transfer_state_to_next_arc: bool=False, processing_settings: MultiArcPropagatorProcessingSettings=...) -> MultiArcPropagatorSettings:
    """Factory function to create multi-arc propagator settings.
	
	Factory function to create multi-arc propagator settings. It works by providing separate settings for
	each arc in a list.
	
	
	:param single_arc_settings:
			List of SingleArcPropagatorSettings objects to use, one for each arc.
	:param transfer_state_to_next_arc:
			It denotes whether whether the initial state of arc N+1 is to be taken from arc N (for N>0).
	:return:
			Multi-arc propagator settings object.
	"""

@typing.overload
def multitype(propagator_settings_list: list[SingleArcPropagatorSettings], termination_settings: ..., output_variables: list[...]=[], print_interval: float=...) -> MultiTypePropagatorSettings:
    ...

@typing.overload
def multitype(propagator_settings_list: list[SingleArcPropagatorSettings], integrator_settings: ..., initial_time: float, termination_settings: ..., output_variables: list[...]=[], processing_settings: SingleArcPropagatorProcessingSettings=...) -> MultiTypePropagatorSettings:
    """Factory function to create multitype propagator settings.
	
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
	"""

def non_sequential_termination(forward_termination_settings: PropagationTerminationSettings, backward_termination_settings: PropagationTerminationSettings) -> PropagationTerminationSettings:
    ...

@typing.overload
def rotational(torque_models: dict[str, dict[str, list[...]]], bodies_to_integrate: list[str], initial_states: numpy.ndarray, termination_settings: ..., propagator: RotationalPropagatorType=..., output_variables: list[...]=[], print_interval: float=...) -> RotationalStatePropagatorSettings:
    ...

@typing.overload
def rotational(torque_models: dict[str, dict[str, list[...]]], bodies_to_integrate: list[str], initial_states: numpy.ndarray, initial_time: float, integrator_settings: ..., termination_settings: ..., propagator: RotationalPropagatorType=..., output_variables: list[...]=[], processing_settings: SingleArcPropagatorProcessingSettings=...) -> RotationalStatePropagatorSettings:
    """Factory function to create rotational state propagator settings.
	
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
	"""

def time_termination(termination_time: float, terminate_exactly_on_final_condition: bool=False) -> PropagationTerminationSettings:
    """Factory function to create time termination settings for the propagation.
	
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
	"""

@typing.overload
def translational(central_bodies: list[str], acceleration_models: dict[str, dict[str, list[..., 3, 1, 0, 3, ...]]], bodies_to_integrate: list[str], initial_states: numpy.ndarray, termination_settings: ..., propagator: TranslationalPropagatorType=..., output_variables: list[...]=[], print_interval: float=...) -> TranslationalStatePropagatorSettings:
    ...

@typing.overload
def translational(central_bodies: list[str], acceleration_models: dict[str, dict[str, list[..., 3, 1, 0, 3, ...]]], bodies_to_integrate: list[str], initial_states: numpy.ndarray, initial_time: float, integrator_settings: ..., termination_settings: ..., propagator: TranslationalPropagatorType=..., output_variables: list[...]=[], processing_settings: SingleArcPropagatorProcessingSettings=...) -> TranslationalStatePropagatorSettings:
    """Factory function to create translational state propagator settings with stopping condition at given final time.
	
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
	"""
cowell: TranslationalPropagatorType
cpu_time_stopping_condition_type: PropagationTerminationTypes
custom_stopping_condition_type: PropagationTerminationTypes
custom_type: StateType
dependent_variable_stopping_condition_type: PropagationTerminationTypes
encke: TranslationalPropagatorType
exponential_map: RotationalPropagatorType
gauss_keplerian: TranslationalPropagatorType
gauss_modified_equinoctial: TranslationalPropagatorType
hybrid_stopping_condition_type: PropagationTerminationTypes
hybrid_type: StateType
mass_type: StateType
modified_rodrigues_parameters: RotationalPropagatorType
quaternions: RotationalPropagatorType
rotational_type: StateType
time_stopping_condition_type: PropagationTerminationTypes
translational_type: StateType
undefined_rotational_propagator: RotationalPropagatorType
undefined_translational_propagator: TranslationalPropagatorType
unified_state_model_exponential_map: TranslationalPropagatorType
unified_state_model_modified_rodrigues_parameters: TranslationalPropagatorType
unified_state_model_quaternions: TranslationalPropagatorType