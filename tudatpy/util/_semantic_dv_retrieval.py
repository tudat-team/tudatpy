''' 
Copyright (c) 2010-2023, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

# General imports
import numpy as np
from inspect import getmro

# Tudat imports
from tudatpy.util import result2array
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation.propagation_setup.dependent_variable import \
    VariableSettings, \
    get_dependent_variable_id, \
    get_dependent_variable_shape


class SimulationVariableDict(dict):
    """
    Dictionary-like class designed to store and retrieve dependent variable time histories
    using dependent variable settings objects (instances of `VariableSettings`-derived
    classes, check the Tudat API docs for more information about this class).

    The goal of this class is to make dependent variable time history retrieval semantic and
    straight-forward:

    ```
    total_acceleration_history = dep_vars_dict[
        propagation_setup.dependent_variable.total_acceleration("Delfi-C3")
    ]
    ```
        
    Example usage:
    --------------

    See the example below. As you can see, we can use as keys either the

    * Original dependent variable settings object (`dependent_variables_to_save[0]`),

    * Or a newly created dependent variable settings object (`propagation_setup.dependent_variable.total_acceleration("Delfi-C3")`)

    ```
    # Create simulation object and propagate the dynamics
    dynamics_simulator = numerical_simulation.create_dynamics_simulator(
        bodies, propagator_settings
    )

    # Create SimulationVariableDict
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
    ```

    How are time histories saved in a `SimulationVariableDict`?
    -----------------------------------------------------------

    A `SimulationVariableDict` maps which maps dependent variables, identified by either their
    corresponding dependent variable settings object (an instance of a `VariableSettings`-derived 
    class) or its string ID, to their time histories.

    The time history of each dependent variable is a Python `dict` which maps epochs (`float`) to
    NumPy arrays (`np.ndarray`) of shape `(A, B)`:

        dict[epoch: np.ndarray[A, B]]
    
    Practical examples:
    
    * `[1, 1]` (a number, most likely a `float`) for a scalar dependent variable
    * `[3, 1]` for a vectorial dependent variable
    * `[3, 3]` for a matrix dependent variable, etc.
    """

    def __read_key(self, key: VariableSettings):
        """
        Read a `SimulationVariableDict` key, as follows:

        * If the `key` is an instance of a `VariableSettings`-derived class, return the
          object's string ID, obtained using `get_dependent_variable_id`.

        * If the `key` is a string, use it directly, assuming it is the string ID of a 
          dependent variable settings object

        * If the `key` is neither type, raise a `TypeError`
        """

        # If the key is an instance of a `VariableSettings`-derived class, return its string ID
        if VariableSettings in getmro(key.__class__):
            key = get_dependent_variable_id(key)
        
        # Else, ensure it is a string, or else raise a TypeError
        if not isinstance(key, str): raise TypeError(
            'SimulationVariableDict keys must be either instances of `VariableSettings`-derived classes, '
            'or dependent variable string IDs (see `numerical_simulation.propagation_setup.dependent_variable.get_dependent_variable_id`).'
        )

        return key

    def __init__(self, mapping=None, /, **kwargs):
        """
        Create a `SimulationVariableDict` from either a dictionary (`mapping`), or a series of 
        keyword-value pairs (`kwargs`).
        """
        if mapping is not None:
            mapping = {
                self.__read_key(key): value for key, value in mapping.items()
            }
        else:
            mapping = {}
        if kwargs:
            mapping.update(
                {self.__read_key(key): value for key, value in kwargs.items()}
            )
        super().__init__(mapping)

        # Create a time history attribute to easily retrieve the time history
        # The time history is obtained from the keys of the first dependent
        # variable's time history dictionary
        self.time_history = np.array(
            # List of keys (epochs) of the time history of the
            list(self[
                # first dependent variable stored in the `SimulationVariableDict`
                list(self.keys())[0]
            ].keys()))

    def __setitem__(self, __key: VariableSettings, __value: any) -> None:
        """
        Set the time history corresponding to a dependent variable, identified either by
        the dependent variable settings object corresponding to the dependent variable
        or its string ID.

        Check the documentation of `SimulationVariableDict.__read_key` for further details.
        """
        return super().__setitem__(self.__read_key(__key), __value)

    def __getitem__(self, __key: VariableSettings) -> dict[float: np.ndarray]:
        """
        Retrieve the time history corresponding to a dependent variable, identified either by
        the dependent variable settings object corresponding to the dependent variable
        or its string ID.

        Check the documentation of `SimulationVariableDict.__read_key` for further details.

        Output
        ------
        
        * Time history of the dependent variable, returned as a `dict` mapping epochs (`float`)
          to `np.ndarray`s containing the value of the dependent variable at each given epoch.
        """
        return super().__getitem__(self.__read_key(__key))

    def __repr__(self) -> str:
        """
        Return a string summary of the contents of a `SimulationVariableDict` for print.
        """
        
        width = max([len(ID) for ID in self.keys()])+5
        title = f'{"Depent Variable Dictionary Summary":^{width}}'

        representation_string = f'\n{"="*width}\n' + title + f'\n{"="*width}\n' + \
            ''.join([
                f'{f"[{i}]":<4} {ID}\n' for i, (ID, value) in enumerate(self.items())
            ]) + \
            f'{"="*width}\n'

        return representation_string


def result2dict(
        dynamics_simulator: numerical_simulation.SingleArcSimulator
    ) -> SimulationVariableDict:
    """
    Construct a dictionary-like object (`SimulationVariableDict`) which maps which maps dependent variables
    to their time histories. See the documentation of `SimulationVariableDict` to learn more about how
    time histories are saved, and how the time history of a given dependent variable can be retrieved.

    Arguments
    ---------
    - dynamics_simulator: `SingleArcSimulator` object containing the results of the numerical propagation

    Output
    ------
    - dependent_variable_dictionary: `SimulationVariableDict` of propagation
    """

    #--------------------------------------------------------------------
    #%% RETRIEVE DEPENDENT VARIABLE DATA
    #--------------------------------------------------------------------

    # Retrieve dependent variable settings objects
    dependent_variable_settings = dynamics_simulator.propagation_results.ordered_dependent_variable_settings

    # Retrieve /transposed/ time and dependent variable histories
    time_history = result2array(
        dynamics_simulator.dependent_variable_history).T[0,  :]
    dependent_variable_history = result2array(
        dynamics_simulator.dependent_variable_history).T[1:, :]

    # Calculate total number of epochs of propagation
    n = len(time_history)

    #--------------------------------------------------------------------
    #%% CONSTRUCT DEPENDENT VARIABLE DICTIONARY
    #--------------------------------------------------------------------
    
    # Construct dependent variable matrices
    dependent_variable_matrices = []
    for ((i, m), dependent_variable) in dependent_variable_settings.items():
        
        # Retrieve dependent variable shape
        A, B = get_dependent_variable_shape(dependent_variable)
        
        # Save dependent variable history as a tensor of (A, B)-sized 
        # matrices with `n` entries, where `n` is the number of epochs
        dependent_variable_matrices.append(
            dependent_variable_history[
                # From index i to index i+m (the flattened dimension of the dependent variable)
                i:i+m, :
            ].T.reshape((n, A, B))
        )

    # Construct dependent variable dictionary
    dependent_variable_dictionary = SimulationVariableDict({
        dependent_variable: {
            epoch: dependent_variable_matrices[i_depv][i_epoch] 
            for i_epoch, epoch in enumerate(time_history)
        }
        for i_depv, dependent_variable in enumerate(dependent_variable_settings.values())
    })

    return dependent_variable_dictionary
