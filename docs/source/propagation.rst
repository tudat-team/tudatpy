``propagation``
===============
Functionalities and utilities of propagation objects.

This module provides functionalities for propagation settings
objects. It also contains some utility functions that extract specific quantities from propagation settings and body
objects. Note that the classes in this module are rarely created manually,
but are instead created by the functionality in the :ref:`\`\`propagation_setup\`\``  module.












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.propagation

.. autosummary::

   get_state_of_bodies

   get_damped_proper_mode_initial_rotational_state

   combine_initial_states

   dependent_variable_dictionary.dependent_variable_dictionary.create_dependent_variable_dictionary


.. autofunction:: tudatpy.numerical_simulation.propagation.get_state_of_bodies

.. autofunction:: tudatpy.numerical_simulation.propagation.get_damped_proper_mode_initial_rotational_state

.. autofunction:: tudatpy.numerical_simulation.propagation.combine_initial_states

.. autofunction:: tudatpy.numerical_simulation.propagation.dependent_variable_dictionary.dependent_variable_dictionary.create_dependent_variable_dictionary



Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.propagation

.. autosummary::

   PropagationTerminationReason



.. autoclass:: tudatpy.numerical_simulation.propagation.PropagationTerminationReason
   :members:




Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.propagation

.. autosummary::

   SimulationResults

   SingleArcSimulationResults

   PropagationTerminationDetails

   PropagationTerminationDetailsFromHybridCondition

   RotationalProperModeDampingResults

   dependent_variable_dictionary.dependent_variable_dictionary.DependentVariableDictionary


.. autoclass:: tudatpy.numerical_simulation.propagation.SimulationResults
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation.SingleArcSimulationResults
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation.PropagationTerminationDetails
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation.PropagationTerminationDetailsFromHybridCondition
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation.RotationalProperModeDampingResults
   :members:

.. autoclass:: tudatpy.numerical_simulation.propagation.dependent_variable_dictionary.dependent_variable_dictionary.DependentVariableDictionary
   :members:

