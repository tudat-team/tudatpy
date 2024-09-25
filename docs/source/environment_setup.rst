``environment_setup``
=====================
Definition of environment settings.

This module contains submodules to define environment settings. It also contains a set of factory functions to use
environment settings in a simulation by creating natural and artificial body objects.










.. toctree::
   :maxdepth: 2
   :caption: Modules

   aerodynamic_coefficients
   atmosphere
   ephemeris
   gravity_field
   gravity_field_variation
   ground_station
   radiation_pressure
   rigid_body
   rotation_model
   shape
   shape_deformation
   vehicle_systems


Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup

.. autosummary::

   get_default_body_settings

   get_default_body_settings_time_limited

   get_default_single_body_settings

   get_default_single_body_settings_time_limited

   get_default_single_alternate_body_settings

   get_default_single_body_settings_time_limited

   add_aerodynamic_coefficient_interface

   create_system_of_bodies

   create_simplified_system_of_bodies

   create_body_ephemeris

   add_radiation_pressure_interface

   add_flight_conditions

   add_rotation_model

   add_mass_properties_model

   add_engine_model

   add_variable_direction_engine_model



.. autofunction:: tudatpy.numerical_simulation.environment_setup.get_default_body_settings

.. autofunction:: tudatpy.numerical_simulation.environment_setup.get_default_body_settings_time_limited

.. autofunction:: tudatpy.numerical_simulation.environment_setup.get_default_single_body_settings

.. autofunction:: tudatpy.numerical_simulation.environment_setup.get_default_single_body_settings_time_limited

.. autofunction:: tudatpy.numerical_simulation.environment_setup.get_default_single_alternate_body_settings

.. autofunction:: tudatpy.numerical_simulation.environment_setup.get_default_single_body_settings_time_limited

.. autofunction:: tudatpy.numerical_simulation.environment_setup.add_aerodynamic_coefficient_interface

.. autofunction:: tudatpy.numerical_simulation.environment_setup.create_system_of_bodies

.. autofunction:: tudatpy.numerical_simulation.environment_setup.create_simplified_system_of_bodies

.. autofunction:: tudatpy.numerical_simulation.environment_setup.create_body_ephemeris

.. autofunction:: tudatpy.numerical_simulation.environment_setup.add_radiation_pressure_interface

.. autofunction:: tudatpy.numerical_simulation.environment_setup.add_flight_conditions

.. autofunction:: tudatpy.numerical_simulation.environment_setup.add_rotation_model

.. autofunction:: tudatpy.numerical_simulation.environment_setup.add_mass_properties_model

.. autofunction:: tudatpy.numerical_simulation.environment_setup.add_engine_model

.. autofunction:: tudatpy.numerical_simulation.environment_setup.add_variable_direction_engine_model






Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup

.. autosummary::

   BodyListSettings

   BodySettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.BodyListSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.BodySettings
   :members:



