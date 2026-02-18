.. _observable_models_setup:

``observable_models_setup``
===========================

This module and its submodules contain settings for creating observation models in Tudat. These observation models are for instance used during a state estimation/orbit determination from real or simulated tracking data.

The :ref:`model_settings` module provides a list of obsevation model settings (range, Doppler, etc.) to use, while the :ref:`biases`, :ref:`light_time_corrections` and :ref:`links` provide interfaces for creating the various input to the observation model settings. A more detailed overview of the different setting types and their interaction is given on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-model-setup.html>`_ .

The main interface with the rest of Tudat for this module is the :func:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` class, which are created by the various functions in the :ref:`model_settings` module.  These settings are used to create observation models either directly using the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.create_observation_simulators` function, or when creating an :class:`~tudatpy.estimation.estimation_analysis.Estimator` object (which creates the observation models internally, and can be extracted using the :attr:`~tudatpy.estimation.estimation_analysis.Estimator.observation_simulators` attribute). Objects of the :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationModel` class are used to compute the observations for a single type and set of link ends, with each type of observation implemented in its dedicated derived class. Models for all link ends specified of a single observable type are stored in an :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` object.

The use of the observation models to simulate observations is described on a dedicated `user guide page <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/creating-observations/simulating-observations.html#simulating-observations>`_

.. toctree::
   :maxdepth: 2
   :caption: Modules

   /estimation/observable_models_setup/biases
   /estimation/observable_models_setup/light_time_corrections
   /estimation/observable_models_setup/links
   /estimation/observable_models_setup/model_settings