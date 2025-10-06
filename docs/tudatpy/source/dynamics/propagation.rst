.. _propagation:

``propagation``
===============
Functionalities and utilities of propagation objects.

This module provides functionalities for propagation settings
objects. It also contains some utility functions that extract specific quantities from propagation settings and body
objects. Note that the classes in this module are rarely created manually,
but are instead created by the functionality in the :ref:`propagation_setup`  module.












Functions
---------
.. currentmodule:: tudatpy.dynamics.propagation

.. autosummary::

   get_state_of_bodies

   get_damped_proper_mode_initial_rotational_state

   combine_initial_states

   dependent_variable_dictionary.create_dependent_variable_dictionary

.. autofunction:: tudatpy.dynamics.propagation.get_state_of_bodies

.. autofunction:: tudatpy.dynamics.propagation.get_damped_proper_mode_initial_rotational_state

.. autofunction:: tudatpy.dynamics.propagation.combine_initial_states

.. autofunction:: tudatpy.dynamics.propagation.dependent_variable_dictionary.create_dependent_variable_dictionary


Enumerations
------------
.. currentmodule:: tudatpy.dynamics.propagation

.. autosummary::

   PropagationTerminationReason



.. autoclass:: tudatpy.dynamics.propagation.PropagationTerminationReason
   :members:




Classes
-------
.. currentmodule:: tudatpy.dynamics.propagation

.. autosummary::

   SimulationResults

   SingleArcSimulationResults

   SingleArcVariationalSimulationResults

   MultiArcSimulationResults

   MultiArcVariationalSimulationResults

   HybridArcSimulationResults

   HybridArcVariationalSimulationResults

   PropagationTerminationDetails

   PropagationTerminationDetailsFromHybridCondition

   RotationalProperModeDampingResults

   dependent_variable_dictionary.DependentVariableDictionary


.. autoclass:: tudatpy.dynamics.propagation.SimulationResults
   :members:

.. autoclass:: tudatpy.dynamics.propagation.SingleArcSimulationResults
   :members:

.. autoclass:: tudatpy.dynamics.propagation.SingleArcVariationalSimulationResults
   :members:

.. autoclass:: tudatpy.dynamics.propagation.MultiArcSimulationResults
   :members:

.. autoclass:: tudatpy.dynamics.propagation.MultiArcVariationalSimulationResults
   :members:

.. autoclass:: tudatpy.dynamics.propagation.HybridArcSimulationResults
   :members:

.. autoclass:: tudatpy.dynamics.propagation.HybridArcVariationalSimulationResults
   :members:

.. autoclass:: tudatpy.dynamics.propagation.PropagationTerminationDetails
   :members:

.. autoclass:: tudatpy.dynamics.propagation.PropagationTerminationDetailsFromHybridCondition
   :members:

.. autoclass:: tudatpy.dynamics.propagation.RotationalProperModeDampingResults
   :members:

.. autoclass:: tudatpy.dynamics.propagation.dependent_variable_dictionary.DependentVariableDictionary
   :members:

