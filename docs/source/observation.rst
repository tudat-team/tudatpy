``observation``
===============
This module contains a set of factory functions for setting up the
observation models, for use in the tudat estimation framework












Functions
---------
.. currentmodule:: tudatpy.numerical_simulation.estimation_setup.observation

.. autosummary::

   link_definition

   body_origin_link_end_id

   body_reference_point_link_end_id

   one_way_downlink_link_ends

   one_way_uplink_link_ends

   light_time_convergence_settings

   first_order_relativistic_light_time_correction

   absolute_bias

   relative_bias

   arcwise_absolute_bias

   arcwise_absolute_bias

   arcwise_relative_bias

   arcwise_relative_bias

   time_drift_bias

   arc_wise_time_drift_bias

   arc_wise_time_drift_bias

   combined_bias

   one_way_range

   n_way_range

   n_way_range_from_one_way_links

   two_way_range

   two_way_range_from_one_way_links

   angular_position

   relative_angular_position

   one_way_doppler_instantaneous

   two_way_doppler_instantaneous

   two_way_doppler_instantaneous_from_one_way_links

   one_way_doppler_averaged

   n_way_doppler_averaged

   n_way_doppler_averaged_from_one_way_links

   two_way_doppler_averaged

   two_way_doppler_averaged_from_one_way_links

   cartesian_position

   cartesian_velocity

   313_euler_angles

   elevation_angle_viability

   elevation_angle_viability_list

   body_avoidance_viability

   body_avoidance_viability_list

   body_occultation_viability

   body_occultation_viability_list

   doppler_ancilliary_settings

   two_way_range_ancilliary_settings

   two_way_doppler_ancilliary_settings

   n_way_range_ancilliary_settings

   n_way_doppler_ancilliary_settings

   tabulated_simulation_settings

   tabulated_simulation_settings_list

   get_default_reference_link_end

   continuous_arc_simulation_settings

   continuous_arc_simulation_settings_list

   add_gaussian_noise_to_all

   add_gaussian_noise_to_observable

   add_gaussian_noise_to_observable_for_link_ends

   add_viability_check_to_all

   add_viability_check_to_observable

   add_viability_check_to_observable_for_link_ends

   add_dependent_variables_to_all

   add_dependent_variables_to_observable

   add_dependent_variables_to_obs_for_links_end

   add_noise_function_to_all

   add_noise_function_to_observable

   add_noise_function_to_observable_for_link_ends



.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.link_definition

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.body_origin_link_end_id

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.body_reference_point_link_end_id

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.one_way_downlink_link_ends

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.one_way_uplink_link_ends

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.light_time_convergence_settings

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.first_order_relativistic_light_time_correction

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.relative_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.arcwise_absolute_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.arcwise_absolute_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.arcwise_relative_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.arcwise_relative_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.time_drift_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.arc_wise_time_drift_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.arc_wise_time_drift_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.combined_bias

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.one_way_range

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.n_way_range

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_from_one_way_links

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.two_way_range

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.two_way_range_from_one_way_links

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.angular_position

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.relative_angular_position

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_instantaneous

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_instantaneous

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_instantaneous_from_one_way_links

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_averaged

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_averaged_from_one_way_links

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_averaged

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_averaged_from_one_way_links

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.cartesian_position

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.cartesian_velocity

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.313_euler_angles

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.elevation_angle_viability

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.elevation_angle_viability_list

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.body_avoidance_viability

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.body_avoidance_viability_list

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.body_occultation_viability

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.body_occultation_viability_list

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.doppler_ancilliary_settings

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.two_way_range_ancilliary_settings

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_ancilliary_settings

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_ancilliary_settings

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_ancilliary_settings

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings_list

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.get_default_reference_link_end

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.continuous_arc_simulation_settings

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.continuous_arc_simulation_settings_list

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_all

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_observable

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_observable_for_link_ends

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_viability_check_to_all

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_viability_check_to_observable

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_viability_check_to_observable_for_link_ends

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_dependent_variables_to_all

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_dependent_variables_to_observable

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_dependent_variables_to_obs_for_links_end

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_noise_function_to_all

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_noise_function_to_observable

.. autofunction:: tudatpy.numerical_simulation.estimation_setup.observation.add_noise_function_to_observable_for_link_ends




Enumerations
------------
.. currentmodule:: tudatpy.numerical_simulation.estimation_setup.observation

.. autosummary::

   LinkEndType

   ObservableType

   ObservationViabilityType

   LightTimeFailureHandling



.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.ObservableType
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilityType
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.LightTimeFailureHandling
   :members:




Classes
-------
.. currentmodule:: tudatpy.numerical_simulation.estimation_setup.observation

.. autosummary::

   LinkEndId

   LinkDefinition

   DopplerProperTimeRateSettings

   ObservationSettings

   OneWayDopplerObservationSettings

   LightTimeCorrectionSettings

   LightTimeConvergenceCriteria

   ObservationBiasSettings

   ObservationSimulationSettings

   TabulatedObservationSimulationSettings

   ObservationViabilitySettings

   ObservationDependentVariableSettings

   ObservationAncilliarySimulationSettings



.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.DopplerProperTimeRateSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.OneWayDopplerObservationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.LightTimeConvergenceCriteria
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings
   :members:

.. autoclass:: tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings
   :members:



