.. _parameters_setup:

``parameters_setup``
===========================



This module contains a set of factory functions for setting up the
observation models, for use in the tudat estimation framework

This module and its constituents are in many cases documented under the assumption that its functionalities are used in the context of an estimation problem.
However, since estimatable parameter settings are firstly used to set up variational equations of the dynamical / observation model w.r.t. the estimatable parameters,
the functionality of this module can be relevant in any context in which variational equations are required.




Functions
---------
.. currentmodule:: tudatpy.dynamics.parameters_setup

.. autosummary::

   create_parameter_set

   initial_states

   constant_drag_coefficient

   arcwise_constant_drag_coefficient

   drag_component_scaling

   side_component_scaling

   lift_component_scaling

   radiation_pressure_coefficient

   arcwise_radiation_pressure_coefficient

   radiation_pressure_target_direction_scaling

   radiation_pressure_target_perpendicular_direction_scaling

   empirical_accelerations

   arcwise_empirical_accelerations

   constant_empirical_acceleration_terms

   full_empirical_acceleration_terms

   arcwise_constant_empirical_acceleration_terms

   quasi_impulsive_shots

   gravitational_parameter

   spherical_harmonics_c_coefficients

   spherical_harmonics_s_coefficients

   spherical_harmonics_c_coefficients_block

   spherical_harmonics_s_coefficients_block

   yarkovsky_parameter

   constant_rotation_rate

   rotation_pole_position

   order_invariant_k_love_number

   order_varying_k_love_number

   mode_coupled_k_love_numbers

   polynomial_gravity_field_variation_amplitudes

   periodic_gravity_field_variation_amplitudes

   monomial_gravity_field_variation_amplitudes

   monomial_full_block_gravity_field_variation_amplitudes

   direct_tidal_dissipation_time_lag

   inverse_tidal_quality_factor

   mean_moment_of_inertia

   periodic_spin_variations

   polar_motion_amplitudes

   scaled_longitude_libration_amplitude

   core_factor

   free_core_nutation_rate

   absolute_observation_bias

   relative_observation_bias

   arcwise_absolute_observation_bias

   arcwise_relative_observation_bias

   ground_station_position

   reference_point_position

   ppn_parameter_gamma

   ppn_parameter_beta

   global_polynomial_clock_corrections

   arc_wise_polynomial_clock_corrections

   time_drift_observation_bias

   arcwise_time_drift_observation_bias

   constant_time_bias

   arcwise_time_bias

   custom_parameter

   custom_analytical_partial

   custom_numerical_partial


.. autofunction:: tudatpy.dynamics.parameters_setup.create_parameter_set

.. autofunction:: tudatpy.dynamics.parameters_setup.initial_states

.. autofunction:: tudatpy.dynamics.parameters_setup.constant_drag_coefficient

.. autofunction:: tudatpy.dynamics.parameters_setup.arcwise_constant_drag_coefficient

.. autofunction:: tudatpy.dynamics.parameters_setup.drag_component_scaling

.. autofunction:: tudatpy.dynamics.parameters_setup.side_component_scaling

.. autofunction:: tudatpy.dynamics.parameters_setup.lift_component_scaling

.. autofunction:: tudatpy.dynamics.parameters_setup.radiation_pressure_coefficient

.. autofunction:: tudatpy.dynamics.parameters_setup.arcwise_radiation_pressure_coefficient

.. autofunction:: tudatpy.dynamics.parameters_setup.radiation_pressure_target_direction_scaling

.. autofunction:: tudatpy.dynamics.parameters_setup.radiation_pressure_target_perpendicular_direction_scaling

.. autofunction:: tudatpy.dynamics.parameters_setup.empirical_accelerations

.. autofunction:: tudatpy.dynamics.parameters_setup.arcwise_empirical_accelerations

.. autofunction:: tudatpy.dynamics.parameters_setup.constant_empirical_acceleration_terms

.. autofunction:: tudatpy.dynamics.parameters_setup.full_empirical_acceleration_terms

.. autofunction:: tudatpy.dynamics.parameters_setup.arcwise_constant_empirical_acceleration_terms

.. autofunction:: tudatpy.dynamics.parameters_setup.quasi_impulsive_shots

.. autofunction:: tudatpy.dynamics.parameters_setup.gravitational_parameter

.. autofunction:: tudatpy.dynamics.parameters_setup.spherical_harmonics_c_coefficients

.. autofunction:: tudatpy.dynamics.parameters_setup.spherical_harmonics_s_coefficients

.. autofunction:: tudatpy.dynamics.parameters_setup.spherical_harmonics_c_coefficients_block

.. autofunction:: tudatpy.dynamics.parameters_setup.spherical_harmonics_s_coefficients_block

.. autofunction:: tudatpy.dynamics.parameters_setup.yarkovsky_parameter

.. autofunction:: tudatpy.dynamics.parameters_setup.constant_rotation_rate

.. autofunction:: tudatpy.dynamics.parameters_setup.rotation_pole_position

.. autofunction:: tudatpy.dynamics.parameters_setup.order_invariant_k_love_number

.. autofunction:: tudatpy.dynamics.parameters_setup.order_varying_k_love_number

.. autofunction:: tudatpy.dynamics.parameters_setup.mode_coupled_k_love_numbers

.. autofunction:: tudatpy.dynamics.parameters_setup.polynomial_gravity_field_variation_amplitudes

.. autofunction:: tudatpy.dynamics.parameters_setup.periodic_gravity_field_variation_amplitudes

.. autofunction:: tudatpy.dynamics.parameters_setup.monomial_gravity_field_variation_amplitudes

.. autofunction:: tudatpy.dynamics.parameters_setup.monomial_full_block_gravity_field_variation_amplitudes

.. autofunction:: tudatpy.dynamics.parameters_setup.direct_tidal_dissipation_time_lag

.. autofunction:: tudatpy.dynamics.parameters_setup.inverse_tidal_quality_factor

.. autofunction:: tudatpy.dynamics.parameters_setup.mean_moment_of_inertia

.. autofunction:: tudatpy.dynamics.parameters_setup.periodic_spin_variations

.. autofunction:: tudatpy.dynamics.parameters_setup.polar_motion_amplitudes

.. autofunction:: tudatpy.dynamics.parameters_setup.scaled_longitude_libration_amplitude

.. autofunction:: tudatpy.dynamics.parameters_setup.core_factor

.. autofunction:: tudatpy.dynamics.parameters_setup.free_core_nutation_rate

.. autofunction:: tudatpy.dynamics.parameters_setup.absolute_observation_bias

.. autofunction:: tudatpy.dynamics.parameters_setup.relative_observation_bias

.. autofunction:: tudatpy.dynamics.parameters_setup.arcwise_absolute_observation_bias

.. autofunction:: tudatpy.dynamics.parameters_setup.arcwise_relative_observation_bias

.. autofunction:: tudatpy.dynamics.parameters_setup.ground_station_position

.. autofunction:: tudatpy.dynamics.parameters_setup.reference_point_position

.. autofunction:: tudatpy.dynamics.parameters_setup.ppn_parameter_gamma

.. autofunction:: tudatpy.dynamics.parameters_setup.ppn_parameter_beta

.. autofunction:: tudatpy.dynamics.parameters_setup.global_polynomial_clock_corrections

.. autofunction:: tudatpy.dynamics.parameters_setup.arc_wise_polynomial_clock_corrections

.. autofunction:: tudatpy.dynamics.parameters_setup.time_drift_observation_bias

.. autofunction:: tudatpy.dynamics.parameters_setup.arcwise_time_drift_observation_bias

.. autofunction:: tudatpy.dynamics.parameters_setup.constant_time_bias

.. autofunction:: tudatpy.dynamics.parameters_setup.arcwise_time_bias

.. autofunction:: tudatpy.dynamics.parameters_setup.custom_parameter

.. autofunction:: tudatpy.dynamics.parameters_setup.custom_analytical_partial

.. autofunction:: tudatpy.dynamics.parameters_setup.custom_numerical_partial


Enumerations
------------
.. currentmodule:: tudatpy.dynamics.parameters_setup

.. autosummary::

   EstimatableParameterTypes

   EmpiricalAccelerationComponents

   EmpiricalAccelerationFunctionalShapes



.. autoclass:: tudatpy.dynamics.parameters_setup.EstimatableParameterTypes
   :members:

.. autoclass:: tudatpy.dynamics.parameters_setup.EmpiricalAccelerationComponents
   :members:

.. autoclass:: tudatpy.dynamics.parameters_setup.EmpiricalAccelerationFunctionalShapes
   :members:


Classes
-------
.. currentmodule:: tudatpy.dynamics.parameters_setup

.. autosummary::

   EstimatableParameterSettings

   CustomAccelerationPartialSettings

.. autoclass:: tudatpy.dynamics.parameters_setup.EstimatableParameterSettings
   :members:

.. autoclass:: tudatpy.dynamics.parameters_setup.CustomAccelerationPartialSettings
   :members:
