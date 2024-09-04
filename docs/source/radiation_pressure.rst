``radiation_pressure``
======================
This module contains a set of factory functions for setting up the
radiation pressure models of celestial bodies in an environment, including relevant models
solar luminosity, solar system body albedo and emissivity, spacecraft surface reaction to radiation pressure.












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.radiation_pressure

.. autosummary::

   constant_luminosity

   cannonball_radiation_target

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

   thermal_emission_blackbody_constant_emissivity

   thermal_emission_blackbody_variable_emissivity

   thermal_emission_angle_based_radiosity

   specular_diffuse_body_panel_reflection

   lambertian_body_panel_reflection

   isotropic_radiation_source

   panelled_radiation_target

   panelled_extended_radiation_source



.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.constant_luminosity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.cannonball_radiation_target

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.irradiance_based_constant_luminosity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.time_variable_luminosity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.irradiance_based_time_variable_luminosity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.constant_surface_property_distribution

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.spherical_harmonic_surface_property_distribution

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.predefined_spherical_harmonic_surface_property_distribution

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.knocke_type_surface_property_distribution

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.predefined_knocke_type_surface_property_distribution

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.custom_surface_property_distribution

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.constant_radiosity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.constant_albedo_surface_radiosity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.variable_albedo_surface_radiosity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.thermal_emission_blackbody_constant_emissivity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.thermal_emission_blackbody_variable_emissivity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.thermal_emission_angle_based_radiosity

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.specular_diffuse_body_panel_reflection

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.lambertian_body_panel_reflection

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.isotropic_radiation_source

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.panelled_radiation_target

.. autofunction:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.panelled_extended_radiation_source




Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.radiation_pressure

.. autosummary::

   KnockeTypeSurfacePropertyDistributionModel

   SphericalHarmonicSurfacePropertyDistribution



.. autoclass:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.KnockeTypeSurfacePropertyDistributionModel
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.SphericalHarmonicSurfacePropertyDistribution
   :members:




Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.radiation_pressure

.. autosummary::

   LuminosityModelSettings

   SurfacePropertyDistributionSettings

   PanelRadiosityModelSettings

   BodyPanelReflectionLawSettings

   RadiationSourceModelSettings

   RadiationPressureTargetModelSettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.LuminosityModelSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.SurfacePropertyDistributionSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.PanelRadiosityModelSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.BodyPanelReflectionLawSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.RadiationSourceModelSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.radiation_pressure.RadiationPressureTargetModelSettings
   :members:



