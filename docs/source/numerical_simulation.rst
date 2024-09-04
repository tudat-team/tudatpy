``numerical_simulation``
========================
Setup and execution of numerical simulations.

This module contains everything related to the setup of a numerical simulation, including submodules to define
environment, propagation, and integration settings. It also contains objects that perform the simulation using such
settings.










.. toctree::
   :maxdepth: 2
   :caption: Modules

   estimation_setup
   environment_setup
   propagation_setup
   estimation
   environment
   propagation


Functions
---------
.. currentmodule:: tudatpy.numerical_simulation

.. autosummary::

   create_dynamics_simulator



.. autofunction:: tudatpy.numerical_simulation.create_dynamics_simulator






Classes
-------
.. currentmodule:: tudatpy.numerical_simulation

.. autosummary::

   SingleArcSimulator

   SingleArcVariationalSimulator

   Estimator



.. autoclass:: tudatpy.numerical_simulation.SingleArcSimulator
   :members:
   :special-members: __init__



.. autoclass:: tudatpy.numerical_simulation.SingleArcVariationalSimulator
   :members:
   :special-members: __init__



.. autoclass:: tudatpy.numerical_simulation.Estimator
   :members:
   :special-members: __init__





