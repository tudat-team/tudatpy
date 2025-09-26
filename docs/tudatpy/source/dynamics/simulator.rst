.. _simulator:

``simulator``
========================
Setup and execution of numerical simulations.

This module contains everything related to the setup of a numerical simulation, including submodules to define
environment, propagation, and integration settings. It also contains objects that perform the simulation using such
settings.




Functions
---------
.. currentmodule:: tudatpy.dynamics.simulator

.. autosummary::

   create_dynamics_simulator

   create_variational_equations_solver



.. autofunction:: tudatpy.dynamics.simulator.create_dynamics_simulator

.. autofunction:: tudatpy.dynamics.simulator.create_variational_equations_solver






Classes
-------
.. currentmodule:: tudatpy.dynamics.simulator

.. autosummary::

   DynamicsSimulator

   SingleArcSimulator

   MultiArcSimulator

   HybridArcSimulator

   VariationalSimulator
   
   SingleArcVariationalSimulator

   CombinedStateTransitionAndSensitivityMatrixInterface



.. autoclass:: tudatpy.dynamics.simulator.DynamicsSimulator
   :members:

.. autoclass:: tudatpy.dynamics.simulator.SingleArcSimulator
   :members:

.. autoclass:: tudatpy.dynamics.simulator.MultiArcSimulator
   :members:

.. autoclass:: tudatpy.dynamics.simulator.HybridArcSimulator
   :members:

.. autoclass:: tudatpy.dynamics.simulator.VariationalSimulator
   :members:

.. autoclass:: tudatpy.dynamics.simulator.SingleArcVariationalSimulator
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface
   :members:






