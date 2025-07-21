.. _aerodynamic_coefficients:

``aerodynamic_coefficients``
============================
This module contains the factory functions for setting up the
aerodynamic coefficient models of (artificial and natural) bodies in an environment.

The main interfaces with Tudat is the :attr:`~tudatpy.numerical_simulation.environment_setup.BodySettings.aerodynamic_coefficient_settings`
attribute of the body settings, which defines settings for the aerodynamic coefficients of a body.
The functions in this submodule are used to create this settings objects.

The coefficient models create from the settings defined by these settings are used by the
:func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.aerodynamic` acceleration model and
:func:`~tudatpy.numerical_simulation.propagation_setup.torque.aerodynamic` torque model.


The functions in this module can be used to define force and/or moment coefficients. The frame in which
the coefficients are to be provided can be defined by the user through the appropriate :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientFrames` input. For instance, when wanting to provide :math:`C_{D},C_{S},C_{L}`,
``positive_body_fixed_frame_coefficients`` is used. When wanting to provide :math:`C_{X},C_{Y},C_{Z}`,
``negative_aerodynamic_frame_coefficients`` is used.


In the definition of moment coefficients, a reference point w.r.t. which aerodynamic moment coefficients are defined must be provided.
This ``moment_reference_point`1 variable is used to calculate the contribution of the aerodynamic force coefficients to the effective moment coefficients
See the ``add_force_contribution_to_moments`` attribute of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for more details. If the present input is set to NaN (as is the default), the reference point is left undefined, and the aerodynamic moments are computed without computing any force coefficient contribution to the moment coefficients.

Within the aerodynamic coefficient settings that are created, the user may also add settings for control surface aerodynamic
coefficient settings. These settings must use the same reference area, reference length and reference point as the vehicle
coefficient settings. See :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.custom_control_surface`
for an example of how to create and set control surface coefficients. If control surface coefficient are set, the
control surface deflection set through :attr:`~tudatpy.numerical_simulation.environment.VehicleSystems.set_control_surface_deflection`
can be used during every time step to control the surface deflection :math:`\delta`. The coefficients (:math:`C_{L}` shown for example)
used during the propagation then become

.. math::
   C_{L}(...)=C_{L,0}(...)+\sum_{i}^{N}\Delta C_{L,i}(..., \delta_{i})

where :math:`C_{L,0}` are the nominal coefficients (with control surfaces) defined by the coefficient settings, :math:`...` denote
the independent variables (except control surface deflection) as a function of which the coefficients are calculated,
:math:`\Delta C_{L,i}` denotes the control surface coefficient increment due to control surface :math:`i` (for :math:`N` control surfaces),
and :math:`\delta_{i}` is the current deflection of this surface.


Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients

.. autosummary::

   constant

   constant_force_and_moment

   custom_aerodynamic_force_coefficients

   custom_aerodynamic_force_and_moment_coefficients

   tabulated

   tabulated_force_only

   tabulated_from_files

   tabulated_force_only_from_files

   scaled_by_constant

   scaled_by_vector

   scaled_by_vector_function

   custom_control_surface

   tabulated_from_files_control_surface



.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.constant

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.constant_force_and_moment

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_coefficients

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_and_moment_coefficients

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated_force_only

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated_from_files

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated_force_only_from_files

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.scaled_by_constant

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.scaled_by_vector

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.scaled_by_vector_function

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.custom_control_surface

.. autofunction:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated_from_files_control_surface




Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients

.. autosummary::

   AerodynamicsReferenceFrameAngles
   
   AerodynamicsReferenceFrames

   AerodynamicCoefficientFrames

   AerodynamicCoefficientsIndependentVariables
   
   AtmosphericCompositionSpecies


.. autoclass:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicsReferenceFrameAngles
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicsReferenceFrames
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientFrames
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientsIndependentVariables
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AtmosphericCompositionSpecies
   :members:



Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients

.. autosummary::

   AerodynamicCoefficientSettings

   ConstantAerodynamicCoefficientSettings



.. autoclass:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ConstantAerodynamicCoefficientSettings
   :members:



