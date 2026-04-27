.. _space_time:

``space_time``
==============

This module contains settings and helper functions for configuring
space-time metric models used in relativistic dynamics.

The main interface with Tudat is the :attr:`~tudatpy.dynamics.environment_setup.BodyListSettings.space_time_settings`
attribute (of type :class:`~tudatpy.dynamics.environment_setup.space_time.SpaceTimePropertiesSettings`) of the body-list settings.
**The functions in this submodule are used to create these settings objects.** When creating a system of bodies (typically using
the :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), an object of type
:class:`~tudatpy.dynamics.environment.SpaceTimeProperties` is created and attached to the resulting
:attr:`~tudatpy.dynamics.environment.SystemOfBodies.space_time_properties` attribute.

Unlike most submodules in :ref:`environment_setup`, the settings in this module are global to the full
:class:`~tudatpy.dynamics.environment.SystemOfBodies` (not properties of individual
:class:`~tudatpy.dynamics.environment.Body` objects). The resulting space-time properties are used by relativistic
models such as :func:`~tudatpy.dynamics.propagation_setup.acceleration.relativistic_correction`,
:func:`~tudatpy.dynamics.propagation_setup.acceleration.einstein_infeld_hofmann`, and relativistic time-converter setup
functions (:func:`~tudatpy.dynamics.environment_setup.direct_relativistic_time_converter_settings`,
:func:`~tudatpy.dynamics.environment_setup.set_relativistic_time_converters`).

The following code block gives an overview of the steps to define, create, and extract system-level space-time properties:

.. code-block:: python

  from tudatpy.dynamics import environment_setup

  # Create body settings
  body_settings = environment_setup.get_default_body_settings( ... )

  # Define global space-time settings (base class type SpaceTimePropertiesSettings)
  body_settings.space_time_settings = environment_setup.space_time.space_time_properties_settings(
      metric_settings=environment_setup.space_time.solar_system_metric_settings(
          first_order_bodies=["Sun", "Earth"],
          second_order_bodies=["Sun"]),
      ppn_parameter_set=environment_setup.space_time.ppn_parameter_set(
          parameter_gamma=1.0,
          parameter_beta=1.0 ) )

  # Create bodies
  bodies = environment_setup.create_system_of_bodies(body_settings)

  # Extract global space-time properties from the SystemOfBodies
  space_time_properties = bodies.space_time_properties

Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.space_time

.. autosummary::

   ppn_parameter_set

   space_time_properties_settings

   schwarzschild_metric_settings

   solar_system_metric_settings



.. autofunction:: tudatpy.dynamics.environment_setup.space_time.ppn_parameter_set

.. autofunction:: tudatpy.dynamics.environment_setup.space_time.space_time_properties_settings

.. autofunction:: tudatpy.dynamics.environment_setup.space_time.schwarzschild_metric_settings

.. autofunction:: tudatpy.dynamics.environment_setup.space_time.solar_system_metric_settings




Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.space_time

.. autosummary::

   PPNParameterSet

   SpaceTimeMetricType

   SpaceTimeMetricSettings

   SpaceTimePropertiesSettings

   SchwarzschildSpaceTimeMetricSettings

   SolarSystemSpaceTimeMetricSettings



.. autoclass:: tudatpy.dynamics.environment_setup.space_time.PPNParameterSet
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.space_time.SpaceTimeMetricType
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.space_time.SpaceTimeMetricSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.space_time.SpaceTimePropertiesSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.space_time.SchwarzschildSpaceTimeMetricSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.space_time.SolarSystemSpaceTimeMetricSettings
   :members:
