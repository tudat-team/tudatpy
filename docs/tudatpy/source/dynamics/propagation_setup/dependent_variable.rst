.. _dependent_variable:

``dependent_variable``
======================
This module provides the functionality for creating dependent variable settings. Dependent variables are quantities
other than the propagated states that are saved when doing a propagation.

Each dependent variable is created via a dedicated function, the settings for each of which is defined by an object of type :class:`~tudatpy.dynamics.propagation_setup.dependent_variable.SingleDependentVariableSaveSettings`.  **The functions in this submodule are used to create these settings objects.** After propagation, the computed dependent variables can be extracted from the :attr:`~tudatpy.dynamics.propagation.SingleArcSimulationResults.dependent_variable_history` attribute of the :class:`~tudatpy.dynamics.propagation.SingleArcSimulationResults`. AN object of this type (which saves all propagation results) is extracted from the dynamics simulator using the  :attr:`~tudatpy.dynamics.simulator.SingleArcSimulator.propagation_results` attribute.

See the example below for the definition of four separate dependent variables:

.. code-block:: python

    dependent_variables_to_save = [
        propagation_setup.dependent_variable.total_acceleration( "Delfi-C3" ),
        propagation_setup.dependent_variable.keplerian_state( "Delfi-C3", "Earth" ),
        propagation_setup.dependent_variable.latitude( "Delfi-C3", "Earth" ),
        propagation_setup.dependent_variable.longitude( "Delfi-C3", "Earth" )
    ]


This list can then be added to the propagation settings by adding them as input argument to the propagator settings functions in the :ref:`propagator` module (such as :func:`~tudatpy.dynamics.propagation_setup.propagator.translational`, :func:`~tudatpy.dynamics.propagation_setup.propagator.rotational` and
:func:`~tudatpy.dynamics.propagation_setup.propagator.multitype`).

Note that *all* dependent variables output is in SI units (meters, radians, seconds). All epochs are provided in seconds since J2000.











Functions
---------
.. currentmodule:: tudatpy.dynamics.propagation_setup.dependent_variable

.. autosummary::

   get_dependent_variable_id

   get_dependent_variable_size

   get_dependent_variable_shape

   mach_number

   altitude

   airspeed

   body_fixed_airspeed_velocity

   body_fixed_groundspeed_velocity

   density

   temperature

   dynamic_pressure

   local_aerodynamic_g_load

   relative_position

   relative_distance

   relative_velocity

   relative_speed

   keplerian_state

   modified_equinoctial_state

   single_acceleration

   single_acceleration_norm

   total_acceleration_norm

   total_acceleration

   single_torque_norm

   single_torque

   total_torque_norm

   total_torque

   spherical_harmonic_terms_acceleration

   spherical_harmonic_terms_acceleration_norm

   aerodynamic_force_coefficients

   aerodynamic_moment_coefficients

   latitude

   geodetic_latitude

   longitude

   heading_angle

   flight_path_angle

   angle_of_attack

   sideslip_angle

   bank_angle

   radiation_pressure

   total_gravity_field_variation_acceleration

   single_gravity_field_variation_acceleration

   single_per_term_gravity_field_variation_acceleration

   inertial_to_body_fixed_rotation_frame

   tnw_to_inertial_rotation_matrix

   rsw_to_inertial_rotation_matrix

   inertial_to_body_fixed_313_euler_angles

   intermediate_aerodynamic_rotation_matrix_variable

   periapsis_altitude

   apoapsis_altitude

   central_body_fixed_spherical_position

   central_body_fixed_cartesian_position

   body_mass

   radiation_pressure_coefficient

   total_mass_rate

   gravity_field_potential

   gravity_field_laplacian_of_potential

   minimum_body_distance

   minimum_visible_station_body_distances

   custom_dependent_variable

   received_irradiance

   received_irradiance_shadow_function



.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.get_dependent_variable_id

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.get_dependent_variable_size

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.get_dependent_variable_shape

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.mach_number

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.altitude

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.airspeed

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.body_fixed_airspeed_velocity

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.body_fixed_groundspeed_velocity

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.density

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.temperature

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.dynamic_pressure

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.local_aerodynamic_g_load

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.relative_position

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.relative_distance

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.relative_velocity

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.relative_speed

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.keplerian_state

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.modified_equinoctial_state

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.single_acceleration

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.single_acceleration_norm

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.total_acceleration_norm

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.total_acceleration

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.single_torque_norm

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.single_torque

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.total_torque_norm

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.total_torque

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.spherical_harmonic_terms_acceleration

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.spherical_harmonic_terms_acceleration_norm

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.aerodynamic_force_coefficients

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.aerodynamic_moment_coefficients

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.latitude

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.geodetic_latitude

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.longitude

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.heading_angle

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.flight_path_angle

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.angle_of_attack

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.sideslip_angle

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.bank_angle

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.radiation_pressure

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.total_gravity_field_variation_acceleration

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.single_gravity_field_variation_acceleration

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.single_per_term_gravity_field_variation_acceleration

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.inertial_to_body_fixed_rotation_frame

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.tnw_to_inertial_rotation_matrix

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.rsw_to_inertial_rotation_matrix

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.inertial_to_body_fixed_313_euler_angles

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.intermediate_aerodynamic_rotation_matrix_variable

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.periapsis_altitude

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.apoapsis_altitude

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.central_body_fixed_spherical_position

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.central_body_fixed_cartesian_position

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.body_mass

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.radiation_pressure_coefficient

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.total_mass_rate

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.gravity_field_potential

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.gravity_field_laplacian_of_potential

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.minimum_body_distance

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.minimum_visible_station_body_distances

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.custom_dependent_variable

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.received_irradiance

.. autofunction:: tudatpy.dynamics.propagation_setup.dependent_variable.received_irradiance_shadow_function




Enumerations
------------
.. currentmodule:: tudatpy.dynamics.propagation_setup.dependent_variable

.. autosummary::

   PropagationDependentVariables



.. autoclass:: tudatpy.dynamics.propagation_setup.dependent_variable.PropagationDependentVariables
   :members:




Classes
-------
.. currentmodule:: tudatpy.dynamics.propagation_setup.dependent_variable

.. autosummary::

   VariableSettings

   SingleDependentVariableSaveSettings

   SingleAccelerationDependentVariableSaveSettings



.. autoclass:: tudatpy.dynamics.propagation_setup.dependent_variable.VariableSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.dependent_variable.SingleDependentVariableSaveSettings
   :members:

.. autoclass:: tudatpy.dynamics.propagation_setup.dependent_variable.SingleAccelerationDependentVariableSaveSettings
   :members:



