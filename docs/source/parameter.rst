``parameter``
=============
This module contains a set of factory functions for setting up the
observation models, for use in the tudat estimation framework

This module and its constituents are in many cases documented under the assumption that its functionalities are used in the context of an estimation problem.
However, since estimatable parameter settings are firstly used to set up variational equations of the dynamical / observation model w.r.t. the estimatable parameters,
the functionality of this module can be relevant in any context in which variational equations are required.












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.estimation_setup.parameter

.. autosummary::

   initial_states

   constant_drag_coefficient

   arcwise_constant_drag_coefficient

   radiation_pressure_coefficient

   arcwise_radiation_pressure_coefficient

   empirical_accelerations

   arcwise_empirical_accelerations

   constant_empirical_acceleration_terms

   arcwise_constant_empirical_acceleration_terms

   quasi_impulsive_shots

   gravitational_parameter

   spherical_harmonics_c_coefficients

   spherical_harmonics_s_coefficients

   spherical_harmonics_c_coefficients_block

   spherical_harmonics_s_coefficients_block

   constant_rotation_rate

   rotation_pole_position

   mean_moment_of_inertia

   periodic_spin_variations

   polar_motion_amplitudes

   core_factor

   free_core_nutation_rate

   absolute_observation_bias

   relative_observation_bias

   arcwise_absolute_observation_bias

   arcwise_relative_observation_bias

   ground_station_position

   ppn_parameter_gamma

   ppn_parameter_beta



.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.initial_states

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.constant_drag_coefficient

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.arcwise_constant_drag_coefficient

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.radiation_pressure_coefficient

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.arcwise_radiation_pressure_coefficient

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.empirical_accelerations

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.arcwise_empirical_accelerations

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.constant_empirical_acceleration_terms

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.arcwise_constant_empirical_acceleration_terms

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.quasi_impulsive_shots

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.gravitational_parameter

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_c_coefficients

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_s_coefficients

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_c_coefficients_block

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_s_coefficients_block

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.constant_rotation_rate

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.rotation_pole_position

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.mean_moment_of_inertia

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.periodic_spin_variations

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.polar_motion_amplitudes

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.core_factor

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.free_core_nutation_rate

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.absolute_observation_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.relative_observation_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.arcwise_absolute_observation_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.arcwise_relative_observation_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.ground_station_position

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.ppn_parameter_gamma

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.parameter.ppn_parameter_beta




Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.estimation_setup.parameter

.. autosummary::

   EstimatableParameterTypes



.. autoclass:: tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterTypes
   :members:




Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.estimation_setup.parameter

.. autosummary::

   EstimatableParameterSettings



.. autoclass:: tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings
   :members:



