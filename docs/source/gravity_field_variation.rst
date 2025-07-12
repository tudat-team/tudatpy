.. _gravity_field_variation:

``gravity_field_variation``
===========================
This module contains a set of factory functions for setting up the
gravity field variations of a (spherical harmonic) gravity field.

The main interfaces with Tudat is the :attr:`~tudatpy.numerical_simulation.environment_setup.BodySettings.gravity_field_variation_settings`
list attribute of the body settings, which defines settings for time-variations of spherical harmonic coefficients of a body.
The functions in this submodule are used to create the settings objects that go into this list.

Once created, the gravity field variation settings defined through the settings in this submodule each compute a :math:`\Delta C_{j,lm}`
and :math:`\Delta \bar{S}_{j,lm}` variation to the cosine and sine coefficients at degree :math:`l` and order :math:`m`.
Each of the :math:`N` gravity field variation models defined for a given body is evaluated at each time step to produce:

.. math::
   \bar{C}_{lm}=\bar{C}_{0,lm}+\sum_{j=1}^{N}\Delta \bar{C}_{j,lm}(t)\\
   \bar{S}_{lm}=\bar{S}_{0,lm}+\sum_{j=1}^{N}\Delta \bar{S}_{j,lm}(t)

where :math:`\bar{C}_{0,lm}` and :math:`\bar{S}_{0,lm}` are the static coefficients (defined through :attr:`~tudatpy.numerical_simulation.environment_setup.BodySettings.gravity_field_settings` for a given body)

Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation

.. autosummary::

   solid_body_tide

   solid_body_tide_complex_k

   solid_body_tide_degree_variable_k

   solid_body_tide_degree_variable_complex_k

   solid_body_tide_degree_order_variable_k

   solid_body_tide_degree_order_variable_complex_k

   solid_multi_body_tide_degree_order_variable_k

   mode_coupled_solid_body_tide

   periodic

   single_period_periodic

   polynomial

   single_power_polynomial

   tabulated

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_complex_k

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_variable_k

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_variable_complex_k

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_complex_k

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_multi_body_tide_degree_order_variable_k

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.mode_coupled_solid_body_tide

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.periodic

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.single_period_periodic

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.polynomial

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.single_power_polynomial

.. autofunction:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.tabulated















Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation

.. autosummary::

   BodyDeformationTypes
   


.. autoclass:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BodyDeformationTypes
   :members:



Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation

.. autosummary::

   GravityFieldVariationSettings

   BasicSolidBodyGravityFieldVariationSettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings
   :members:



