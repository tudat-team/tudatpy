.. _environment_setup:

``environment_setup``
=====================
Definition of environment settings.

This module contains submodules to define environment settings. It also contains a set of factory functions to use
environment settings in a simulation by creating natural and artificial body objects.










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



