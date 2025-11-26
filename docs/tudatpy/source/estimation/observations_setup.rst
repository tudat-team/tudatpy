.. _observations_setup:

``observations_setup``
===========================

This module and its submodules contain settings for simulating/loading observation models in Tudat. These observations are for instance used as input for a state estimation/orbit determination.

The :ref:`observations_simulation_setup` module provides the top-level functionality for defining obsevarvation simulation settings, while the :ref:`observations_dependent_variables`, :ref:`ancillary_settings;`, :ref:`random_noise` and :ref:`viability` modules provide interfaces for creating the various input to customize these settings. The :ref:`observations_wrapper` module provides various interfaces to load data from various sources. A more detailed overview of generating observations in Tudat is given on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/creating-observations.html>`_ .

The main interfaces with the rest of Tudat for this module are the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` class, which are created by the various functions in the :ref:`observations_simulation_setup` module and define settings for simulating observations, and :class:`~tudatpy.estimation.observations.ObservationCollection` objects, which are produced by various functions in the :ref:`observations_wrapper` module and contain the actual observations used downstream in Tudat.

``observations_setup``
======================

.. toctree::
   :maxdepth: 2
   :caption: Modules

   /estimation/observations_setup/ancillary_settings
   /estimation/observations_setup/observations_dependent_variables
   /estimation/observations_setup/observations_simulation_setup
   /estimation/observations_setup/observations_wrapper
   /estimation/observations_setup/random_noise
   /estimation/observations_setup/viability