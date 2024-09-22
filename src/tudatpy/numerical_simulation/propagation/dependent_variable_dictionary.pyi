"""

Copyright (c) 2010-2023, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
"""
from inspect import getmro
import numpy as np
import numpy
import tudatpy.numerical_simulation.expose_numerical_simulation
from tudatpy.numerical_simulation.expose_numerical_simulation import SingleArcSimulator
import tudatpy.numerical_simulation.propagation_setup.dependent_variable.expose_dependent_variable
from tudatpy.numerical_simulation.propagation_setup.dependent_variable.expose_dependent_variable import VariableSettings
from tudatpy.numerical_simulation.propagation_setup.dependent_variable.expose_dependent_variable import get_dependent_variable_id
from tudatpy.numerical_simulation.propagation_setup.dependent_variable.expose_dependent_variable import get_dependent_variable_shape
from tudatpy.util._support import result2array
import typing
from typing import Any
__all__ = ['Any', 'DependentVariableDictionary', 'SingleArcSimulator', 'VariableSettings', 'create_dependent_variable_dictionary', 'get_dependent_variable_id', 'get_dependent_variable_shape', 'getmro', 'np', 'result2array']

class DependentVariableDictionary(dict):
    """Dictionary-like class designed to store and retrieve dependent variable time histories
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
	```
	
	How are time histories saved in a `DependentVariableDictionary`?
	-----------------------------------------------------------
	
	A `DependentVariableDictionary` maps which maps dependent variables, identified by either their
	corresponding dependent variable settings object (an instance of a `VariableSettings`-derived
	class) or its string ID, to their time histories.
	
	The time history of each dependent variable is a Python `dict` which maps epochs (`float`) to
	NumPy arrays (`np.ndarray`) of shape `(A, B)`:
	
		dict[epoch: np.ndarray[A, B]]
	
	**Important**: in `(A, B)`, we remove singleton/trivial dimensions (dimensions, `A` or `B`, of size 1).
	In the case of scalar dependent variables, the value associated to each epoch is a `np.ndarray` of shape `(1,)`.
	In the case of vectorial dependent variables, it is a **row** vector of size `(A,)`. This is done by using
	`np.squeeze` to remove any dimensions of size 1. Practical examples:
	
	Dimensions of dependent variable values associated to each epoch based on their type:
	
	+-----------+-------------+
	| Data Type | Shape	   |
	+===========+=============+
	| Scalar	| `(1,)`	  |
	| Vectorial | `(3,)`	  |
	| Matrix	| `(A, B)`	|
	| Tensor	| `(A, B, C)` |
	+===========+=============+
	"""

    def _DependentVariableDictionary__read_key(self, key: tudatpy.numerical_simulation.propagation_setup.dependent_variable.expose_dependent_variable.VariableSettings):
        """
        
                Read a `DependentVariableDictionary` key, as follows:
        
                * If the `key` is an instance of a `VariableSettings`-derived class, return the
                  object's string ID, obtained using `get_dependent_variable_id`.
        
                * If the `key` is a string, use it directly, assuming it is the string ID of a
                  dependent variable settings object
        
                * If the `key` is neither type, raise a `TypeError`
                
        """

    def __getitem__(self, _DependentVariableDictionary__key: tudatpy.numerical_simulation.propagation_setup.dependent_variable.expose_dependent_variable.VariableSettings):
        """
        
                Retrieve the time history corresponding to a dependent variable, identified either by
                the dependent variable settings object corresponding to the dependent variable
                or its string ID.
        
                Check the documentation of `DependentVariableDictionary.__read_key` for further details.
        
                Output
                ------
        
                * Time history of the dependent variable, returned as a `dict` mapping epochs (`float`)
                  to `np.ndarray`s containing the value of the dependent variable at each given epoch.
                
        """

    def __init__(self, mapping=None, **kwargs):
        """
        
                Create a `DependentVariableDictionary` from either a dictionary (`mapping`), or a series of
                keyword-value pairs (`kwargs`).
                
        """

    def __repr__(self) -> str:
        """
        
                Return a string summary of the contents of a `DependentVariableDictionary` for print.
                
        """

    def __setitem__(self, _DependentVariableDictionary__key: tudatpy.numerical_simulation.propagation_setup.dependent_variable.expose_dependent_variable.VariableSettings, _DependentVariableDictionary__value: typing.Any) -> None:
        """
        
                Set the time history corresponding to a dependent variable, identified either by
                the dependent variable settings object corresponding to the dependent variable
                or its string ID.
        
                Check the documentation of `DependentVariableDictionary.__read_key` for further details.
                
        """

    def asarray(self, key: tudatpy.numerical_simulation.propagation_setup.dependent_variable.expose_dependent_variable.VariableSettings) -> numpy.ndarray:
        """
        
                Return the time history of a given dependent variable as a NumPy array.
        
                Arguments
                ---------
                - key: dependent variable settings object or string ID of the dependent variable
        
                Output
                ------
                - time_history: time history of the dependent variable, returned as a NumPy array
                
        """

def create_dependent_variable_dictionary(dynamics_simulator: tudatpy.numerical_simulation.expose_numerical_simulation.SingleArcSimulator) -> DependentVariableDictionary:
    """Construct a dictionary-like object (`DependentVariableDictionary`) which maps which maps dependent variables
	to their time histories. See the documentation of `DependentVariableDictionary` to learn more about how
	time histories are saved, and how the time history of a given dependent variable can be retrieved.
	
	Arguments
	---------
	- dynamics_simulator: `SingleArcSimulator` object containing the results of the numerical propagation
	
	Output
	------
	- dependent_variable_dictionary: `DependentVariableDictionary` of propagation
	"""