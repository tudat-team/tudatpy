.. _propagation_setup:

``propagation_setup``
=====================
Definition of propagation settings.

This module contains submodules to define propagation settings. It also contains a set of factory functions to create
dynamical models for translational state, rotational state, and mass rate.










.. toctree::
   :maxdepth: 2
   :caption: Modules

   /dynamics/propagation_setup/acceleration
   /dynamics/propagation_setup/dependent_variable
   /dynamics/propagation_setup/integrator
   /dynamics/propagation_setup/mass_rate
   /dynamics/propagation_setup/propagator
   /dynamics/propagation_setup/torque
   /dynamics/propagation_setup/thrust


Functions
---------
.. currentmodule:: tudatpy.dynamics.propagation_setup

.. autosummary::

   create_acceleration_models

   create_torque_models

   create_mass_rate_models



.. autofunction:: tudatpy.dynamics.propagation_setup.create_acceleration_models

.. autofunction:: tudatpy.dynamics.propagation_setup.create_torque_models

.. autofunction:: tudatpy.dynamics.propagation_setup.create_mass_rate_models







