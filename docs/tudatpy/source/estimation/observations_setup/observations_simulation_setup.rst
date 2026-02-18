.. _observations_simulation_settings:

``observations_simulation_settings``
====================================

This module contains a set of functions for creating settings for simulating observations in Tudat. The observation simulation settings are stored in objects of type :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings`. Note that the actual models that compute the observations are *not* defined in this module, but in the :ref:`model_settings` module. The functionality in this module defines how to use the observation models to simulate observations (epochs at which to simulate, noise levels, observability constraints, *etc.*).

The :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects can be created using some of the factory functions on this module. Most notably the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings` and :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.continuous_arc_simulation_settings` functions. Alternatively, if a set of actual observations are already available (loaded from real data, or from a simulation), the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.observation_settings_from_collection` function below can be used to generated observation simulation settings with settings (e.g. obsevation epochs, link ends, etc.) identical to those in the real data.

The primary link with the rest of Tudat is that the :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects are used as input to the :func`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` to generate simulated observations. More details and examples on the procedure to simulate observations are given on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/creating-observations/simulating-observations.html#simulating-observations>`_.


Functions
---------
.. currentmodule:: tudatpy.estimation.observations_setup.observations_simulation_settings

.. autosummary::

   tabulated_simulation_settings

   tabulated_simulation_settings_list

   continuous_arc_simulation_settings

   continuous_arc_simulation_settings_list

   observation_settings_from_collection

   change_simulation_settings_observable_types

   create_observation_simulators

.. autofunction:: tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings

.. autofunction:: tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings_list

.. autofunction:: tudatpy.estimation.observations_setup.observations_simulation_settings.continuous_arc_simulation_settings

.. autofunction:: tudatpy.estimation.observations_setup.observations_simulation_settings.continuous_arc_simulation_settings_list

.. autofunction:: tudatpy.estimation.observations_setup.observations_simulation_settings.observation_settings_from_collection

.. autofunction:: tudatpy.estimation.observations_setup.observations_simulation_settings.change_simulation_settings_observable_types

.. autofunction:: tudatpy.estimation.observations_setup.observations_simulation_settings.create_observation_simulators

Classes
-------
.. currentmodule:: tudatpy.estimation.observations_setup.observations_simulation_settings

.. autosummary::

   ObservationSimulationSettings

   TabulatedObservationSimulationSettings

.. autoclass:: tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings
   :members:

.. autoclass:: tudatpy.estimation.observations_setup.observations_simulation_settings.TabulatedObservationSimulationSettings
   :members: