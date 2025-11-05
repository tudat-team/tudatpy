.. _radiation_pressure:

``radiation_pressure``
======================
This module contains a set of factory functions for setting up the
radiation pressure models of bodies in an environment, including relevant models
solar luminosity, solar system body albedo and emissivity, spacecraft surface reaction to radiation pressure.


The main interfaces with Tudat are the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.radiation_source_settings` and
:attr:`~tudatpy.dynamics.environment_setup.BodySettings.radiation_pressure_target_settings` attributes of the body settings, which define
settings for radiation pressure source and target settings of a body. **The functions in this submodule are used to create these settings objects.**. When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), objects of type
:class:`~tudatpy.dynamics.environment.RadiationPressureTargetModel` and :class:`~tudatpy.dynamics.environment.RadiationPressureSourceModel`  (or a derived class) are created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings objects, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.radiation_pressure_source_model` and
:attr:`~tudatpy.dynamics.environment.Body.radiation_pressure_target_models` attributes.

Note that, while a body can have only a
single source model (which may include multiple contributions such as albedo and planetary radiation pressure), a body can have
multiple target models, and the target model that is used can be specified when choosing the specific :func:`tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure` acceleration

For isotropic source models, the :func:`~tudatpy.dynamics.environment_setup.radiation_pressure.isotropic_radiation_source` is used, which requires an input of type :class:`tudatpy.dynamics.environment_setup.radiation_pressure.LuminosityModelSettings` (or
derived class) to define an irradiance model. For bodies with a variable emission over the surface (albedo, infrared), the :func:`~tudatpy.dynamics.environment_setup.radiation_pressure.panelled_extended_radiation_source` model is
used, which requires a list of :class:`tudatpy.dynamics.environment_setup.radiation_pressure.PanelRadiosityModelSettings` (or derived classes) to define an irradiance model.

More details on the link between different aspects of radiation pressure in Tudat are described on `a dedicated page <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/radiation_pressure_acceleration.html>`_. in the user guide.


Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.radiation_pressure

.. autosummary::

   constant_luminosity

   irradiance_based_constant_luminosity

   time_variable_luminosity

   irradiance_based_time_variable_luminosity

   constant_surface_property_distribution

   spherical_harmonic_surface_property_distribution

   predefined_spherical_harmonic_surface_property_distribution

   knocke_type_surface_property_distribution

   predefined_knocke_type_surface_property_distribution

   custom_surface_property_distribution

   constant_radiosity

   constant_albedo_surface_radiosity

   variable_albedo_surface_radiosity

   thermal_emission_angle_based_radiosity

   thermal_emission_blackbody_constant_emissivity

   thermal_emission_blackbody_variable_emissivity

   specular_diffuse_body_panel_reflection

   lambertian_body_panel_reflection

   isotropic_radiation_source

   panelled_extended_radiation_source

   cannonball_radiation_target

   panelled_radiation_target



.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.constant_luminosity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.irradiance_based_constant_luminosity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.time_variable_luminosity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.irradiance_based_time_variable_luminosity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.constant_surface_property_distribution

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.spherical_harmonic_surface_property_distribution

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.predefined_spherical_harmonic_surface_property_distribution

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.knocke_type_surface_property_distribution

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.predefined_knocke_type_surface_property_distribution

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.custom_surface_property_distribution

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.constant_radiosity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.constant_albedo_surface_radiosity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.variable_albedo_surface_radiosity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.thermal_emission_angle_based_radiosity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.thermal_emission_blackbody_constant_emissivity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.thermal_emission_blackbody_variable_emissivity

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.specular_diffuse_body_panel_reflection

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.lambertian_body_panel_reflection

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.isotropic_radiation_source

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.panelled_extended_radiation_source

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.cannonball_radiation_target

.. autofunction:: tudatpy.dynamics.environment_setup.radiation_pressure.panelled_radiation_target




Enumerations
------------
.. currentmodule:: tudatpy.dynamics.environment_setup.radiation_pressure

.. autosummary::

   KnockeTypeSurfacePropertyDistributionModel

   SphericalHarmonicsSurfacePropertyDistributionModel

   RadiationPressureType



.. autoclass:: tudatpy.dynamics.environment_setup.radiation_pressure.KnockeTypeSurfacePropertyDistributionModel
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.radiation_pressure.SphericalHarmonicsSurfacePropertyDistributionModel
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.radiation_pressure.RadiationPressureType
   :members:


Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.radiation_pressure

.. autosummary::

   LuminosityModelSettings

   SurfacePropertyDistributionSettings

   PanelRadiosityModelSettings

   BodyPanelReflectionLawSettings

   RadiationSourceModelSettings

   RadiationPressureTargetModelSettings



.. autoclass:: tudatpy.dynamics.environment_setup.radiation_pressure.LuminosityModelSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.radiation_pressure.SurfacePropertyDistributionSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.radiation_pressure.PanelRadiosityModelSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.radiation_pressure.BodyPanelReflectionLawSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.radiation_pressure.RadiationSourceModelSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.radiation_pressure.RadiationPressureTargetModelSettings
   :members:



