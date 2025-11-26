.. _observations_dependent_variables:

``observations_dependent_variables``
====================================

As is the case with the state propagation (see :ref:`dependent_variables`), you can define any number of so-called dependent variable to be computed along with the observations themselves. These include distances between link ends, angles between link ends, and a variety of other options. These quantities can be useful in understanding the behaviour of e.g. observation noise, tracking geometry or coverage.

The functions in this module are used create objects of type :class:`~tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings` that contain settings for these observation dependent variables. The main interface with Tudat is that these objects are can be added to :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects, which are in turn used as input to observation simulation (see :ref:`observations_simulation_settings`). The dependent variable settings can be set in these objects directly using the :attr:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings.dependent_variable_settings_list` attribute, or by using the :func:`~tudatpy.estimation.observations_setup.observations_dependent_variables.add_dependent_variables_to_all`, :func:`~tudatpy.estimation.observations_setup.observations_dependent_variables.add_dependent_variables_to_observable` and :func:`~tudatpy.estimation.observations_setup.observations_dependent_variables.add_dependent_variables_check_to_observable_for_link_ends` functions (which are functions of convenience to add the same ancillary settings to a set of observation simulation settings)

After observation simulation, the dependent variables can be extracted from the :class:`~tudatpy.estimation.observations.ObservationCollection` or :class:`~tudatpy.estimation.observations.SingleObservationSet` objects.

Functions
---------
.. currentmodule:: tudatpy.estimation.observations_setup.observations_dependent_variables

.. autosummary::

   elevation_angle_dependent_variable

   azimuth_angle_dependent_variable

   target_range_between_link_ends_dependent_variable

   avoidance_angle_dependent_variable

   body_center_distance_dependent_variable
   
   body_limb_distance_dependent_variable

   angle_wrt_orbital_plane_dependent_variable

   integration_time_dependent_variable

   retransmission_delays_dependent_variable

   add_dependent_variables_to_all

   add_dependent_variables_to_observable

   add_dependent_variables_to_observable_for_link_ends

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.elevation_angle_dependent_variable

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.azimuth_angle_dependent_variable

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.target_range_between_link_ends_dependent_variable

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.avoidance_angle_dependent_variable

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.body_center_distance_dependent_variable
   
.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.body_limb_distance_dependent_variable

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.angle_wrt_orbital_plane_dependent_variable

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.integration_time_dependent_variable

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.retransmission_delays_dependent_variable

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.add_dependent_variables_to_all

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.add_dependent_variables_to_observable

.. autofunction:: tudatpy.estimation.observations_setup.observations_dependent_variables.add_dependent_variables_to_observable_for_link_ends

Enumerations
------------
.. currentmodule:: tudatpy.estimation.observations_setup.observations_dependent_variables

.. autosummary::

   IntegratedObservationPropertyHandling


.. autoclass:: tudatpy.estimation.observations_setup.observations_dependent_variables.IntegratedObservationPropertyHandling
   :members:




Classes
-------
.. currentmodule:: tudatpy.estimation.observations_setup.observations_dependent_variables

.. autosummary::

   ObservationDependentVariableSettings

.. autoclass:: tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
   :members: