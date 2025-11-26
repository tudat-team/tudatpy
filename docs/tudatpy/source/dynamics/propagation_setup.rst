.. _propagation_setup:

``propagation_setup``
=====================

This module contains submodules to define settings for the numerical propagation of states of natural and artificail bodies.
A detailed overview of the propagation framework in Tudat can be found in the `user guide <https://docs.tudat.space/en/latest/user-guide/state-propagation/propagation-setup.html>`_ .

As a short overview of how the submodules below are used in the context of a propagation:

* Defining state derivative model settings for propagation of translational, rotational and mass dynamics: :ref:`acceleration`, :ref:`torque` and :ref:`mass_rate`, with thrust-specific models in :ref:`thrust`
* Defining settings for the numerical integration: :ref:`integrator`
* Defining settings for saving additional (dependent) variables in addition to the state variables: :ref:`dependent_variable`
* Combining settings into a full propagation settings: :ref:`propagator`

More details on the procedure and options in creating environment models and bodies can be found in our `user guide <https://docs.tudat.space/en/latest/user-guide/state-propagation/environment-setup.html>`_. For the use of the environment models *during*
a numerical propagation (for instance for custom models) see `here <https://docs.tudat.space/en/latest/user-guide/state-propagation/environment-setup/custom-models/environment-during-propagation.html#environment-during-propagation>`_.

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







