.. _dynamics:

``dynamics``
============

This submodule contains the functionality for doing numerical state propagation in Tudat, which includes the models for physical environment of bodies (ephemerides, rotation models, gravity fiels, etc.).

* :ref:`environment_setup`/ :ref:`environment`: Functionality related to the physical environment (properties of natural and celestial bodies)
* :ref:`propagation_setup`/ :ref:`propagation`: Functionality related to numerical propagation of states (state types, acceleration models, output variables, *etc.*)
* :ref:`parameters_setup`/ :ref:`parameters`: Functionality related to the definition of parameters for which variational equations can be propagated (initial states, gravity field coefficients, etc.)
* :ref:`simulator` Functionality to combine the models and settings listed above to perform the actual numerical propagation.

The distinction between the ``environment`` and ``environment_setup`` libraries is the following (and similar for ``propagation`` and ``parameters``):

* The ``environment_setup`` submodule contains no actual functionality to perform any calculations. It contains a long list of functions to create *settings* that are used to create the models that do the actual calculations.
* The ``environment``  submodule contains the functionality to perform the actual calculations. Typically, the objects in this submodule are created from a list of objects created in the ``environment_setup`` library.

The reason for this split is that the environment models can have various interdependencies which are difficult to manually implement, but straightforward to conceptually define with a string, boolean, etc. in a settings object. A Tudat function (in this case the :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function) then takes care of mapping the user-definition of the models defined through ``environment_setup`` functionality to the actual objects performing calculations in the ``environment`` model. As one example: it is easy to state that a set of aerodynamic coefficients dependent on angle of attack (this is defined in the ``environment_setup`` submodule), while it is rather cumbersome to manually extract the angle of attack, and input it to the aerodynamic coefficient during every time step.

Each of the submpdules provide a description of how the link between the setup layer and the calculation layer is made, where settings can be defined, and where the resulting model can be extracted and used.

An overview of how to use the overall state propagation functionality can be found in the `user guide <https://docs.tudat.space/en/latest/user-guide/state-propagation.html>`_

.. toctree::
   :maxdepth: 2
   :caption: Modules

   /dynamics/environment_setup
   /dynamics/propagation_setup
   /dynamics/parameters_setup
   /dynamics/environment
   /dynamics/propagation
   /dynamics/parameters
   /dynamics/simulator
