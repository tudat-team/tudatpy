import numpy
import pybind11_stubgen.typing_ext
from ...astro import time_representation
from ...dynamics import environment
from ...dynamics.propagation_setup import dependent_variable
from ...dynamics.propagation_setup import propagator
import typing
__all__ = ['AccelerationModel', 'AerodynamicGuidance', 'ConstantThrustMagnitudeWrapper', 'CustomThrustMagnitudeWrapper', 'DependentVariablesInterface', 'HybridArcSimulationResults', 'HybridArcVariationalSimulationResults', 'MassRateModel', 'MultiArcSimulationResults', 'MultiArcVariationalSimulationResults', 'PropagationTerminationDetails', 'PropagationTerminationDetailsFromHybridCondition', 'PropagationTerminationReason', 'RotationalProperModeDampingResults', 'SimulationResults', 'SingleArcSimulationResults', 'SingleArcVariationalSimulationResults', 'ThrustMagnitudeWrapper', 'TorqueModel', 'combine_initial_states', 'get_damped_proper_mode_initial_rotational_state', 'get_generalized_acceleration_size', 'get_initial_rotational_state_of_body', 'get_initial_state_of_bodies', 'get_initial_state_of_body', 'get_single_integration_differential_equation_order', 'get_single_integration_size', 'get_state_of_bodies', 'nan_or_inf_detected_in_state', 'propagation_never_run', 'runtime_error_caught_in_propagation', 'termination_condition_reached', 'unknown_reason']

class AccelerationModel:
    pass

class AerodynamicGuidance:
    angle_of_attack: float
    bank_angle: float
    sideslip_angle: float

    def __init__(self) -> None:
        ...

    def updateGuidance(self, current_time: float) -> None:
        ...

class ConstantThrustMagnitudeWrapper(ThrustMagnitudeWrapper):
    constant_thrust_magnitude: float

class CustomThrustMagnitudeWrapper(ThrustMagnitudeWrapper):
    pass

class DependentVariablesInterface:
    """No documentation found."""

class HybridArcSimulationResults(SimulationResults):
    """Class that stores all the results (including logging data) of a hybrid-arc propagation of states.
    
    Class that stores all the results (including logging data) of a hybrid-arc propagation of states.
    The results of the constituent arcs are accessed through the ``single_arc_results`` and ``multi_arc_results`` attribute1."""

    @property
    def multi_arc_results(self) -> MultiArcSimulationResults:
        """
                 **read-only**
        
                 Object with results from the multi-arc propagation component of the hybrid-arc propagation. This object stores
                 all numerical and logging results from the multi-arc component.
        
                 :type: MultiArcSimulationResults
        """

    @property
    def single_arc_results(self) -> SingleArcSimulationResults:
        """
                 **read-only**
        
                 Object with results from the single-arc propagation component of the hybrid-arc propagation. This object stores
                 all numerical and logging results from the single-arc component.
        
                 :type: SingleArcSimulationResults
        """

class HybridArcVariationalSimulationResults(SimulationResults):
    """Class that stores all the results (including logging data) of a hybrid-arc propagation of states and variational equations.
    
    Class that stores all the results (including logging data) of a hybrid-arc propagation of states and variational equations.
    The results of the constituent arcs are accessed through the ``single_arc_results`` and ``multi_arc_results`` attribute1."""

    @property
    def multi_arc_results(self) -> MultiArcVariationalSimulationResults:
        """
                 **read-only**
        
                 Object with results from the multi-arc propagation component of the hybrid-arc propagation. This object stores
                 all numerical and logging results from the multi-arc component.
        
                 :type: MultiArcVariationalSimulationResults
        """

    @property
    def single_arc_results(self) -> SingleArcVariationalSimulationResults:
        """
                 **read-only**
        
                 Object with results from the single-arc propagation component of the hybrid-arc propagation. This object stores
                 all numerical and logging results from the single-arc component.
        
                 :type: SingleArcVariationalSimulationResults
        """

class MassRateModel:
    pass

class MultiArcSimulationResults(SimulationResults):
    """Class that stores all the results (including logging data) of a multi-arc propagation of states.
    
    Class that stores all the results (including logging data) of a multi-arc propagation of states.
    The results of the constituent arcs are accessed through the ``single_arc_results`` attribute."""

    @property
    def arc_end_times(self) -> list[float]:
        """
                 **read-only**
        
                 List of epochs at which each of the arcs end (e.g. latest epoch in state history per arc)
        
                 :type: list[float]
        """

    @property
    def arc_start_times(self) -> list[float]:
        """
                 **read-only**
        
                 List of epochs at which each of the arcs start (e.g. earliest epoch in state history per arc)
        
                 :type: list[float]
        """

    @property
    def propagation_is_performed(self) -> bool:
        """
                 **read-only**
        
                 Boolean indicating whether the propagation for which this object stores the results has been performed or not
        
                 :type: int
        """

    @property
    def single_arc_results(self) -> list[SingleArcSimulationResults]:
        """
                 **read-only**
        
                 List of results from the single-arc propagations. The list of single-arc results objects store
                 all numerical and logging results from the constituent arc propagations. These single-arc result
                 objects are where the state history, etc. can be retrieved
        
                 :type: list[SingleArcSimulationResults]
        """

class MultiArcVariationalSimulationResults(SimulationResults):
    """Class that stores all the results (including logging data) of a multi-arc propagation of states and variational equations
    
    Class that stores all the results (including logging data) of a multi-arc propagation of states and variational equations. The
    results of the constituent arcs are accessed through the ``single_arc_results`` attribute."""

    @property
    def propagation_is_performed(self) -> bool:
        """
                 **read-only**
        
                 Boolean indicating whether the propagation for which this object stores the results has been performed or not
        
                 :type: int
        """

    @property
    def single_arc_results(self) -> list[SingleArcVariationalSimulationResults]:
        """
         R"doc(
                 **read-only**
        
                 List of results from the single-arc propagations. The list of single-arc results objects store
                 all numerical and logging results from the constituent arc propagations. These single-arc result
                 objects are where the state history, etc. can be retrieved
        
                 :type: list[SingleArcVariationalSimulationResults]
        """

class PropagationTerminationDetails:
    """Class that provides information on the reason for the
    termination of the propagation."""

    @property
    def terminated_on_exact_condition(self) -> bool:
        """
                 Boolean defining whether the propagation was terminated on an *exact* final condition,
                 or once the propagation went *past* the determined final condition. The choice of behaviour is
                 defined by the termination settings provided as input to the Simulator object. This variable only
                 has a meaningful definition if the ``termination_reason`` has value ``termination_condition_reached``
        
        
                 :type: bool
        """

    @property
    def termination_reason(self) -> PropagationTerminationReason:
        """
                 Enum defining the reason the propagation was terminated
        
        
                 :type: PropagationTerminationReason
        """

class PropagationTerminationDetailsFromHybridCondition(PropagationTerminationDetails):
    """Class that provides information on the reason for the termination of the propagation, for hybrid termination conditions
    
    
    Derived class from :class:`PropagationTerminationDetails` that provides information on the reason for the termination of the propagation,
    for the case of hybrid termination conditions (defined using the :func:`~tudatpy.dynamics.propagation_setup.propagator.hybrid_termination`)
    function"""

    @property
    def was_condition_met_when_stopping(self) -> list[bool]:
        """
                 List of booleans defining, per entry in ``termination_settings`` when calling :func:`~tudatpy.dynamics.propagation_setup.propagator.hybrid_termination`,
                 whether the corresponding entry of the hybrid termination settings was met or not.
        
        
                 :type: list[bool]
        
                 Examples
                 --------
        
                 Assuming the hybrid termination settings were defined similar to the example in the :func:`~tudatpy.dynamics.propagation_setup.propagator.hybrid_termination`:
        
                 .. code-block:: python
        
                     # Store termination setting objects in a list
                     termination_settings_list = [
                         time_termination_settings,
                         altitude_termination_settings,
                         cpu_termination_settings,
                     ]
                     # Define string representations for output
                     termination_conditions_repr = [
                         "Time Termination",
                         "Altitude Termination",
                         "CPU Time Termination",
                     ]
        
                     # Create hybrid termination settings
                     termination_settings = propagation_setup.propagator.hybrid_termination(
                         termination_settings_list, fulfill_single_condition=True
                     )
        
                 The ``was_condition_met_when_stopping`` attribute contains a list of booleans, with the same length as the
                 ``termination_settings_list``. The following code can be used to print if the termination condition was met for each
                 of the termination settings:
        
                 .. code-block:: python
        
                     # perform propagation
                     ...
        
                     # post-process propagation results
                     termination_details = dynamics_simulator.propagation_results.termination_details
                     condition_met_flags = termination_details.was_condition_met_when_stopping
        
                     condition_fulfilled = [
                         f"{condition:<35}: {met}"
                         for condition, met in zip(termination_conditions_repr, condition_met_flags)
                     ]
        
                     print("Termination Conditions fulfilled:")
                     print("\\n".join(condition_fulfilled))
        """

class PropagationTerminationReason:
    """Enumeration of types of termination of propagation.
    
          
    
    Members:
    
      propagation_never_run
    
      unknown_reason
    
      termination_condition_reached
    
      runtime_error_caught_in_propagation
    
      nan_or_inf_detected_in_state"""
    __members__: typing.ClassVar[dict[str, PropagationTerminationReason]]
    nan_or_inf_detected_in_state: typing.ClassVar[PropagationTerminationReason]
    propagation_never_run: typing.ClassVar[PropagationTerminationReason]
    runtime_error_caught_in_propagation: typing.ClassVar[PropagationTerminationReason]
    termination_condition_reached: typing.ClassVar[PropagationTerminationReason]
    unknown_reason: typing.ClassVar[PropagationTerminationReason]

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

class RotationalProperModeDampingResults:
    """Object that stores the results of the algorithm to damp the proper mode of rotational dynamics for an initial state,
    as computed by the :func:`~get_damped_proper_mode_initial_rotational_state` function"""

    @property
    def damped_initial_state(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                 Initial state produced by the damping algorithm, for which the signature of the proper mode should be
                 removed (or at least, substantially reduced). Note that this initial state corresponds to the *full* state vector
                 that is provided to the ``get_damped_proper_mode_initial_rotational_state`` function (e.g. is size 7
                 for rotational dynamics of a single body, size 13 for coupled orbital-rotational dynamics of a single body, etc.)
        
        
                 :type: numpy.ndarray
        """

    @damped_initial_state.setter
    def damped_initial_state(self, arg0: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        ...

    @property
    def forward_backward_dependent_variables(self) -> list[tuple[dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]]:
        """
                 As ``forward_backward_states``, but for the dependent variables.
        
        
                 :type: list[tuple[dict[float,numpy.ndarray],dict[float,numpy.ndarray]]]
        """

    @forward_backward_dependent_variables.setter
    def forward_backward_dependent_variables(self, arg0: list[tuple[dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]]) -> None:
        ...

    @property
    def forward_backward_states(self) -> list[tuple[dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]]:
        """
                 Data structure that contains the full state histories used by the damping algorithm. The contents are are as follows:
        
                 * The :math:`i^{th}` entry of the list corresponds to the :math:`i^{th}` iteration of the forward-backward propagation
                 * Each tuple in the list contains two dictionaries, the first one corresponding to the forward propagation results, the seconds one to the backward propagation results
        
        
                 :type: list[tuple[dict[float,numpy.ndarray],dict[float,numpy.ndarray]]]
        """

    @forward_backward_states.setter
    def forward_backward_states(self, arg0: list[tuple[dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]]) -> None:
        ...

class SimulationResults:
    """Base class for objects that store all results of a numerical propagation.
    
    Base class for objects that store all results of a numerical propagation.
    Derived class are implemented for single-, multi- and hybrid-arc propagation of both dynamics and variational equations"""

    @property
    def dependent_variable_interface(self) -> DependentVariablesInterface:
        """
                 **read-only**
        
                 Attrribute that allows for automatic interpolation and retrieval of dependent variables
        
                 :type: DependentVariablesInterface
        """

class SingleArcSimulationResults(SimulationResults):
    """Class that stores all the results (including logging data) of a single-arc propagation of states"""

    def clear_data(self) -> None:
        """
                 Function to delete the contents of this object.
        
                 Function to delete the contents of this object. This function is typically called when wanting to manually reduce the
                 memory usage in large analyses by clearing data. It does not need to be used manually before repropagating.
        """

    @property
    def cumulative_computation_time_history(self) -> dict[float, float]:
        """
                 **read-only**
        
                 History of cumulative computation time in seconds needed during the propagation as key-value
                 pairs. At each epoch (key, as a float) the computation time (value) in seconds is the total computation time
                 used up to and including that time step. This includes the total time up to and including the current time step,
                 since the beginning of the (single-arc) propagation.
        
        
                 :type: dict[float, float]
        """

    @property
    def cumulative_number_of_function_evaluations_history(self) -> dict[float, int]:
        """
                 **read-only**
        
                 History of cumulative number of function evaluations (e.g. number of times the state derivative model is evaluated)
                 as key-value pairs. At each epoch (key, as a float) the cumulative number of function evaluations (value) including the current one is included.
        
                 .. note:: At present, this dictionary is updated during each function evaluation, including the minor stages of (for instance) a multi-stage integrator. This is unlike the rest of the data in this object, which is provided only after each full time step. For variable-time step integration, this also means that the values need not be monotonically increasing with propagation time (if a step size is rejected and recomputed with smaller step size).
        
        
        
                 :type: dict[float, float]
        """

    @property
    def dependent_variable_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Dependent variables computed during the propagation as key-value pairs. The key denotes the epoch as a float. If the output
                 with the higher-resolution :func:`~Time` object is required, use the :func:`~dependent_variable_history_time_object` function.
        
                 The vector of all dependent variables concatenated into a single vector as value, with the epoch as key.
                 They order of the concatenated dependent variables in a single value is provided by the ``dependent_variable_ids`` attribute of this object.
        
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def dependent_variable_history_time_object(self) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Same as :func:`~dependent_variable_history`, but using the high-resolution :func:`~Time` type used as independent variable in the propagation as key
        
                 :type: dict[Time, numpy.ndarray]
        """

    @property
    def dependent_variable_ids(self) -> dict[tuple[int, int], str]:
        """
                 **read-only**
        
                 Key-value container with the starting entry of the dependent variables saved (key), along with associated ID (value).
        
        
                 :type: dict[[int,int], str]
        """

    @property
    def initial_and_final_times(self) -> tuple[time_representation.Time, time_representation.Time]:
        """
                 **read-only**
        
                 Initial and final time of the propagation
        
                 :type: tuple[float,float]
        """

    @property
    def integration_completed_successfully(self) -> bool:
        """
                 **read-only**
        
                 Boolean defining whether the last propagation was finished
                 successfully, as defined by the termination conditions, or if
                 it was terminated prematurely (for instance due to an
                 exception, or an Inf/NaN state entry being detected).
        
        
                 :type: bool
        """

    @property
    def ordered_dependent_variable_settings(self) -> dict[tuple[int, int], dependent_variable.SingleDependentVariableSaveSettings]:
        ...

    @property
    def processed_state_ids(self) -> dict[tuple[int, int], str]:
        """
                 **read-only**
        
                 Key-value container with the starting entry of the states (key), along with associated ID (value).
        
        
                 :type: dict[[int,int] str]
        """

    @property
    def propagated_state_ids(self) -> dict[tuple[int, int], str]:
        """
                 **read-only**
        
                 Key-value container with the starting entry of the states (key), along with associated ID (value).
        
        
                 :type: dict[[int,int] str]
        """

    @property
    def propagated_state_vector_length(self) -> int:
        """
                 **read-only**
        
                 Length of the propagated state vector
        
                 :type: int
        """

    @property
    def propagation_is_performed(self) -> bool:
        """
                 **read-only**
        
                 Boolean indicating whether the propagation for which this object stores the results has been performed or not
        
                 :type: int
        """

    @property
    def state_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Numerical solution of the equations of motion as key-value pairs. The key denotes the epoch as a float. If the output
                 with the higher-resolution :func:`~Time` object is required, use the :func:`~state_history_time_object` function.
        
                 The values of this dictionary contains the
                 numerically calculated state at this epoch. For this function, the states are always converted to so-called
                 'processed' formulations (e.g. Cartesian states for translational dynamics), see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/processed_propagated_elements.html>`_
                 for details. For the history of the states that were actually propagated, use the ``unprocessed_state_history``.
        
                 .. note:: The propagated state at each epoch contains the state types in the following order: Translational ( **T** ), Rotational ( **R** ), Mass ( **M** ), and Custom ( **C** ).
                           When propagating two bodies, an example of what the output state would look like is for instance:
                           [ **T** Body 1, **T** Body 2, **R** Body 1, **R** Body 2, **M** Body 1, **M** Body 2 ] The specifics can be retrieved using the :attr:`state_ids` attribute of this class
        
                 .. note:: For propagation of translational dynamics using cowell
                           propagator, the conventional and propagated
                           coordinates are identical.
        
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def state_history_float(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        ...

    @property
    def state_history_float_split(self) -> tuple[list[float], list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]:
        """
                 **read-only**
        
                 Same as ``state_history``, but with a return type of two lists (one with times and one with states)
        
                 :type: tuple[ list[float], list[numpy.ndarray] ]
        """

    @property
    def state_history_time_object(self) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Same as :func:`~state_history`, but using the high-resolution :func:`~Time` type used as independent variable in the propagation as key
        
        
                 :type: dict[Time, numpy.ndarray]
        """

    @property
    def termination_details(self) -> PropagationTerminationDetails:
        """
                 **read-only**
        
                 Object describing the details of the event that triggered the termination of the last propagation.
        
        
                 :type: PropagationTerminationDetails
        """

    @property
    def total_computation_time(self) -> float:
        """
                 **read-only**
        
                 Total computation time (in seconds) that was required for the propagation.
        
                 :type: float
        """

    @property
    def total_number_of_function_evaluations(self) -> float:
        """
                 **read-only**
        
                 Total number of function evaluations that were computed for the propagation
        
                 :type: float
        """

    @property
    def unordered_dependent_variable_settings(self) -> list[dependent_variable.SingleDependentVariableSaveSettings]:
        ...

    @property
    def unprocessed_state_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Numerical solution of the equations of motion as key-value pairs, without any processing applied.
                 The key denotes the epoch as a float. If the output
                 with the higher-resolution :func:`~Time` object is required, use the :func:`~unprocessed_state_history_time_object` function.
        
                 The values contain the
                 numerically calculated state at this epoch. This attribute contains the states of the propagated bodies expressed in the
                 "raw" form in which the propagation took place. For instance, when using a Gauss-Kepler propagation scheme, this
                 attribute will contain the numerically propagated Keplerian elements at each time epoch
        
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def unprocessed_state_history_time_object(self) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Same as :func:`~unprocessed_state_history`, but using the high-resolution :func:`~Time` type used as independent variable in the propagation as key
        
        
                 :type: dict[Time, numpy.ndarray]
        """

class SingleArcVariationalSimulationResults(SimulationResults):
    """Class that stores all the results (including logging data) of a single-arc propagation of states and variational equations
    
    Class that stores all the results (including logging data) of a single-arc propagation of states and variational equations. The
    propagation results of the states are accessed through the ``dynamics_results`` attribute."""

    @property
    def dynamics_results(self) -> SingleArcSimulationResults:
        """
               **read-only**
        
                 Object with all results of the propagation results of the states, as well as details on propagation termination reason, runtime, function evaluations, etc.
        
        
                 :type: SingleArcSimulationResults
        """

    @property
    def sensitivity_matrix_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
        """
                 **read-only**
        
                 Numerical solution of the sensitivity matrix :math:`\\partial \\mathbf{x}/\\partial\\mathbf{p}` from the variational equations as key-value pairs. The key denotes the epoch. The value contains the
                 numerically calculated sensitivity matrix at this epoch. The state :math:`\\mathbf{x}` in the definition is always in the so-called
                 'processed' formulations (e.g. Cartesian states for translational dynamics), see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/processed_propagated_elements.html>`_
                 for details. See :attr:`~SingleArcSimulationResults.state_history` for the definition of the order of the elements of the states.
        
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def state_transition_matrix_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
        """
                 **read-only**
        
                 Numerical solution of the state transition matrix :math:`\\partial \\mathbf{x}/\\partial\\mathbf{x}_{0}` from the variational equations as key-value pairs. The key denotes the epoch. The value contains the
                 numerically calculated state transition matrix at this epoch. The states :math:`\\mathbf{x}` and :math:`\\mathbf{x}_{0}` in the definition are always in the so-called
                 'processed' formulations (e.g. Cartesian states for translational dynamics), see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/processed_propagated_elements.html>`_
                 for details. See :attr:`~SingleArcSimulationResults.state_history` for the definition of the order of the elements of the states.
        
        
                 :type: dict[float, numpy.ndarray]
        """

class ThrustMagnitudeWrapper:
    pass

class TorqueModel:
    pass

def combine_initial_states(propagator_settings_per_type: dict[propagator.StateType, list[propagator.SingleArcPropagatorSettings]]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
    """Function to retrieve the initial state for a list of propagator settings.
    
    Function to retrieve the initial state for a list of propagator settings. This way, the initial state for
    different quantities to be propagated (e.g., translational state, rotational state, mass) are retrieved and
    organized in a single container.
    
    
    Parameters
    ----------
    propagator_settings_per_type : dict
        Propagator settings where the type of propagation is reported as key and the respective list of propagator settings as value.
    Returns
    -------
    numpy.ndarray
        Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the vector of SingleArcPropagatorSettings of given type."""

def get_damped_proper_mode_initial_rotational_state(bodies: environment.SystemOfBodies, propagator_settings: propagator.SingleArcPropagatorSettings, body_mean_rotational_rate: float, dissipation_times: list[float], propagate_undamped: bool=True) -> RotationalProperModeDampingResults:
    """Function to compute an initial rotational state for which the proper mode of rotation is damped.
    
    Function to compute an initial rotational state for which the proper mode of rotation is damped, using the algorithm
    used by Rambaux et al. (2010) to compute an initial rotational state for Phobos. This algorithm propagates the
    dynamics of the system a number of times, with the settings specified by the user and a specific modification to
    damp the proper mode. Since a number of propagations are performed by this function, it may take some time to run.
    Specifically, the algorithm works as follows:
    
    * Introduce a damping torque (see below) to damp the proper mode, with damping time :math:`\\tau_{d}`
    * Propagate the dynamics forward in time for a period of :math:`10\\tau_{d}`
    * Remove the virtual torque, and propagate the dynamics back to the initial time :math:`t_{0}`
    * Repeat the above for the list of damping times provided by the user
    
    The state after the final backwards propagation to :math:`t_{0}` is provided as output by this function, to be
    used as damped initial state. The output from this function also provides the user access to the full state history
    and dependent variable history of the forward and backward propagations, to allow a user to track and validate
    the progress of the algorithm.
    
    The damping torque :math:`\\Gamma` is defined as follows:
    
    .. math::
       \\boldsymbol{\\Gamma}= -\\frac{1}{\\tau_{d}}\\mathbf{I}\\begin{pmatrix}\\omega_{x}\\\\ \\omega_{y}\\\\ \\omega_{x}-\\omega_{p} \\end{pmatrix}
    
    where :math:\\mathbf{I}` is the body's inertia tensor (in its body-fixed frame), :math:`\\tau_{d}` the damping time of the
    current propagation, and :math:`\\omega_{x}, \\omega_{y}, \\omega_{z}` the body's current rotation about its
    body-fixed, x-, y- and z-axes, respectively. The damping torque is implemented to damp out all rotations along
    the body-fixed x- and y-axes, and any deviations from constant rotation with frequency :\\omega_{p}: about the body-fixed z-axis.
    
    .. note:: The mean rotation rate of the body :math:`\\omega_{p}` is a user-defined input, and must be tuned to the dynamics of the system.
    
    
    Parameters
    ----------
    bodies : SystemOfBodies
        Set of body objects that defines the environment
    propagator_settings : SingleArcPropagatorSettings
        Propagator settings for the dynamics of which the initial rotational state is to be damped. These propagator
        settings must be for rotational dynamics only, or for multi-type rotational dynamics that contains rotational
        dynamics for a single body (e.g. translational-rotational dynamics for a single body)
    
    body_mean_rotational_rate : float
        Mean rotational rate :math:`\\omega_{p}` to which the damping algorithm will force the body-fixed rotation about its z-axis.
    dissipation_times : list[ float ]
        List of damping times :math:`\\tau_{d}` for which the algorithm is to be run. Note that this list should be organized in ascending order for the algorithm to perform properly
    propagate_undamped : bool, default = True
        Boolean defining whether the first forward/backward propagation performed by the damping algorithm has damping turned off (damping turned off if True, damping turned on if False).
        Propagating without any damping before starting the damping algorithm is useful for verification purposes, but not required for the algorithm itself.
    
    Returns
    -------
    DampedInitialRotationalStateResults
        Object that contains the results of the damping algorithm (final damped rotational state, and forward/backward propagation results)."""

def get_generalized_acceleration_size(state_type: propagator.StateType) -> int:
    ...

def get_initial_rotational_state_of_body(body_to_propagate: str, base_orientation: str, bodies: environment.SystemOfBodies, initial_time: time_representation.Time) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
    ...

def get_initial_state_of_bodies(bodies_to_propagate: list[str], central_bodies: list[str], body_system: environment.SystemOfBodies, initial_time: time_representation.Time) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
    ...

def get_initial_state_of_body(body_to_propagate: str, central_body: str, bodies: environment.SystemOfBodies, initial_time: time_representation.Time) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
    ...

def get_single_integration_differential_equation_order(state_type: propagator.StateType) -> int:
    ...

def get_single_integration_size(state_type: propagator.StateType) -> int:
    ...

def get_state_of_bodies(bodies_to_propagate: list[str], central_bodies: list[str], body_system: environment.SystemOfBodies, initial_time: time_representation.Time) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
    """Function to get the translational states of a set of bodies, with respect to some set of central bodies, at the requested time.
    
    Function to get the translational states of a set of bodies, with respect to some set of central bodies, at the requested time. This function
    is typically used to extract an initial state for a propagation of a set of bodies, for which the initial state is extracted from the
    existing ephemerides of the bodies.
    
    
    Parameters
    ----------
    bodies_to_propagate : list[str]
        List of names of bodies for which the state is to be extracted
    central_bodies : list[str]
        List of central bodies, w.r.t. which the states are to be computed (in the same order as ``bodies_to_propagate``)
    bodies_to_propagate : SystemOfBodies
        System of bodies that define the environment
    initial_time : float
        Time at which the states are to be extracted from the environment
    Returns
    -------
    numpy.ndarray
        Vector of size :math:`6\\times N`, with the translational states of each entry of body from
        ``bodies_to_propagate`` w.r.t. the corresponding central body in ``central_bodies``."""
nan_or_inf_detected_in_state: PropagationTerminationReason
propagation_never_run: PropagationTerminationReason
runtime_error_caught_in_propagation: PropagationTerminationReason
termination_condition_reached: PropagationTerminationReason
unknown_reason: PropagationTerminationReason