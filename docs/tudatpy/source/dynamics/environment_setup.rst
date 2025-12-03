.. _environment_setup:

``environment_setup``
=====================

This module contains submodules to define environment settings. In Tudat, the definition of 'environment' is the collection
of all bodies (natural and artificial) and their physical properties. Each of the submodules below contains a set of
functions to create a specific type of environment settings.

This module contains the :class:`~tudatpy.dynamics.environment_setup.BodyListSettings`, which is the contained that holds
all the settings for the bodies that are to be created and used in a simulation/analysis. The
:class:`~tudatpy.dynamics.environment_setup.BodyListSettings` stores a list of :class:`~tudatpy.dynamics.environment_setup.BodySettings`
objects, each of which can be endowed with a list of environment (and system) settings, settings for which
are defined in the modules listed below . The creation of the actual bodies used in the simulation/analysis from the
settings defined in this module is done using the :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function.

More details on the procedure and options in creating environment models and bodies can be found in our `user guide <https://docs.tudat.space/en/latest/user-guide/state-propagation/environment-setup.html>`_. For the use of the environment models *during*
a numerical propagation (for instance for custom models) see `here <https://docs.tudat.space/en/latest/user-guide/state-propagation/environment-setup/custom-models/environment-during-propagation.html#environment-during-propagation>`_.


.. toctree::
   :maxdepth: 2
   :caption: Modules

   /dynamics/environment_setup/aerodynamic_coefficients
   /dynamics/environment_setup/atmosphere
   /dynamics/environment_setup/ephemeris
   /dynamics/environment_setup/gravity_field
   /dynamics/environment_setup/gravity_field_variation
   /dynamics/environment_setup/ground_station
   /dynamics/environment_setup/radiation_pressure
   /dynamics/environment_setup/rigid_body
   /dynamics/environment_setup/rotation_model
   /dynamics/environment_setup/shape
   /dynamics/environment_setup/shape_deformation
   /dynamics/environment_setup/vehicle_systems


Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup

.. autosummary::

   get_default_body_settings

   get_default_body_settings_time_limited

   get_default_single_body_settings

   get_default_single_body_settings_time_limited

   get_default_single_alternate_body_settings

   create_system_of_bodies

   create_simplified_system_of_bodies

   create_body_ephemeris

   add_aerodynamic_coefficient_interface

   add_radiation_pressure_target_model

   add_flight_conditions

   add_rotation_model

   add_mass_properties_model

   add_rigid_body_properties

   add_engine_model

   add_variable_direction_engine_model



.. autofunction:: tudatpy.dynamics.environment_setup.get_default_body_settings

.. autofunction:: tudatpy.dynamics.environment_setup.get_default_body_settings_time_limited

.. autofunction:: tudatpy.dynamics.environment_setup.get_default_single_body_settings

.. autofunction:: tudatpy.dynamics.environment_setup.get_default_single_body_settings_time_limited

.. autofunction:: tudatpy.dynamics.environment_setup.get_default_single_alternate_body_settings

.. autofunction:: tudatpy.dynamics.environment_setup.create_system_of_bodies

.. autofunction:: tudatpy.dynamics.environment_setup.create_simplified_system_of_bodies

.. autofunction:: tudatpy.dynamics.environment_setup.create_body_ephemeris

.. autofunction:: tudatpy.dynamics.environment_setup.add_aerodynamic_coefficient_interface

.. autofunction:: tudatpy.dynamics.environment_setup.add_radiation_pressure_target_model

.. autofunction:: tudatpy.dynamics.environment_setup.add_flight_conditions

.. autofunction:: tudatpy.dynamics.environment_setup.add_rotation_model

.. autofunction:: tudatpy.dynamics.environment_setup.add_mass_properties_model

.. autofunction:: tudatpy.dynamics.environment_setup.add_rigid_body_properties

.. autofunction:: tudatpy.dynamics.environment_setup.add_engine_model

.. autofunction:: tudatpy.dynamics.environment_setup.add_variable_direction_engine_model






Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup

.. autosummary::

   BodyListSettings

   BodySettings



.. autoclass:: tudatpy.dynamics.environment_setup.BodyListSettings
   :members:

   .. automethod:: __init__

.. autoclass:: tudatpy.dynamics.environment_setup.BodySettings
   :members:
   :exclude-members: radiation_pressure_settings



