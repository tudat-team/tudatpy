import typing
import numpy as np
from _typeshed import Incomplete
from tudatpy.dynamics.propagation_setup.dependent_variable import VariableSettings as VariableSettings, get_dependent_variable_id as get_dependent_variable_id, get_dependent_variable_shape as get_dependent_variable_shape
from tudatpy.dynamics.simulator import SingleArcSimulator as SingleArcSimulator
from tudatpy.util import result2array as result2array
from typing import Any

class DependentVariableDictionary(dict):
    """Dictionary-like class designed to store and retrieve dependent variable time histories
    using dependent variable settings objects (instances of `VariableSettings`-derived
    classes, check the Tudat API docs for more information about this class).
    
    The goal of this class is to make dependent variable time history retrieval semantic and
    straight-forward:
    
    .. code-block:: python
    
        total_acceleration_history = dep_vars_dict[
            propagation_setup.dependent_variable.total_acceleration("Delfi-C3")
        ]
    
    Example usage:
    --------------
    
    See the example below. As you can see, we can use as keys either the
    
    * Original dependent variable settings object (`dependent_variables_to_save[0]`),
    * Or a newly created dependent variable settings object (`propagation_setup.dependent_variable.total_acceleration("Delfi-C3")`)
    
    .. code-block:: python
    
        # Create simulation object and propagate the dynamics
        dynamics_simulator = dynamics.simulator.create_dynamics_simulator(
            bodies, propagator_settings
        )
    
        # Create DependentVariableDictionary
        dep_vars_dict = result2dict(dynamics_simulator)
    
        # Retrieve the time history (in `dict[epoch: value]` form) of the total acceleration experienced by Delft-C3
        total_acceleration_history = dep_vars_dict[
            # This can be done using either the `SingleAccelerationDependentVariableSaveSettings`
            # corresponding to this dependent variable
            dependent_variables_to_save[0]
        ]
        total_acceleration_history = dep_vars_dict[
            # Or a newly created one
            propagation_setup.dependent_variable.total_acceleration("Delfi-C3")
        ]
    
    How are time histories saved in a `DependentVariableDictionary`?
    ----------------------------------------------------------------
    
    A `DependentVariableDictionary` maps which maps dependent variables, identified by either their
    corresponding dependent variable settings object (an instance of a `VariableSettings`-derived
    class) or its string ID, to their time histories.
    
    The time history of each dependent variable is a Python `dict` which maps epochs (`float`) to
    NumPy arrays (`np.ndarray`) of shape `(A, B)`: dict[epoch: np.ndarray[A, B]].
    
    **Important**: in `(A, B)`, we remove singleton/trivial dimensions (dimensions, `A` or `B`, of size 1).
    In the case of scalar dependent variables, the value associated to each epoch is a `np.ndarray` of shape `(1,)`.
    In the case of vectorial dependent variables, it is a **row** vector of size `(A,)`. This is done by using
    `np.squeeze` to remove any dimensions of size 1. Practical examples:
    
    Dimensions of dependent variable values associated to each epoch based on their type:
    
    +-----------+-------------+
    | Data Type | Shape       |
    +===========+=============+
    | Scalar    | `(1,)`      |
    +-----------+-------------+
    | Vectorial | `(3,)`      |
    +-----------+-------------+
    | Matrix    | `(A, B)`    |
    +-----------+-------------+
    | Tensor    | `(A, B, C)` |
    +-----------+-------------+"""
    time_history: Incomplete

    def __init__(self, mapping: Incomplete | None=None, /, **kwargs) -> None:
        """
        Create a `DependentVariableDictionary` from either a dictionary (`mapping`), or a series of
        keyword-value pairs (`kwargs`).
        """

    def __setitem__(self, __key: VariableSettings, __value: Any) -> None:
        """
        Set the time history corresponding to a dependent variable, identified either by
        the dependent variable settings object corresponding to the dependent variable
        or its string ID.

        Check the documentation of `DependentVariableDictionary.__read_key` for further details.
        """

    def __getitem__(self, __key: VariableSettings):
        """
        Retrieve the time history corresponding to a dependent variable, identified either by
        the dependent variable settings object corresponding to the dependent variable
        or its string ID.

        Check the documentation of `DependentVariableDictionary.__read_key` for further details.

        Output
        ------
            Time history of the dependent variable, returned as a `dict` mapping epochs (`float`)
            to `np.ndarray`s containing the value of the dependent variable at each given epoch.

        """

    def asarray(self, key: VariableSettings) -> np.ndarray:
        """
        Return the time history of a given dependent variable as a NumPy array.

        Arguments
        ---------
        key : VariableSettings
            dependent variable settings object or string ID of the dependent variable

        Returns
        ------
        time_history : np.ndarray
            time history of the dependent variable, returned as a NumPy array
        """

def create_dependent_variable_dictionary(dynamics_simulator: SingleArcSimulator) -> DependentVariableDictionary:
    """Construct a dictionary-like object (`DependentVariableDictionary`) which maps which maps dependent variables
    to their time histories. See the documentation of `DependentVariableDictionary` to learn more about how
    time histories are saved, and how the time history of a given dependent variable can be retrieved.
    
    Arguments
    ---------
    dynamics_simulator : SingleArcSimulator
        `SingleArcSimulator` object containing the results of the numerical propagation
    
    Returns
    ------
    dependent_variable_dictionary : DependentVariableDictionary
        `DependentVariableDictionary` of propagation"""