.. _estimation:

``estimation``
==============

This submodule contains the functionality for doing numerical state and parameter estimation in Tudat, including models for observations and settings for simulating observations

* :ref:`observable_models_setup`/ :ref:`observable_models`: Functionality related to observation models (e.g. a model for one-way range, two way Doppler, etc.)
* :ref:`observations_setup`/ :ref:`observations`: Functionality related to loading observations (from various source) and simulating observations (e.g. defining settings how to use the observation models using modules from the previous point)
* :ref:`estimation_analysis` Functionality to combine the models and settings listed above (combined with the functionality in :ref:`dynamics`) to perform the actual state and parameter estimation

The distinction between the ``observable_models`` and ``observable_models_setup`` libraries is the following (and similar for ``observations``/``observations_setup``)

* The ``observable_models_setup`` submodule contains no actual functionality to perform any calculations. It contains a long list of functions to create *settings* that are used to create the observation models that do the actual calculations.
* The ``observable_models``  submodule contains the functionality to compute the actual observations. Typically, the objects in this submodule are created from a list of settings objects created in the ``observable_models_setup`` library.

The reason for this split is that the observation models can have many complex dependencies which are difficult to manually implement, but straightforward to conceptually define with a string, boolean, etc. in a settings object. A Tudat function (in this case the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.create_observation_simulators` function) then takes care of mapping the user-definition of the models defined through ``observable_models_setup`` functionality to the actual objects performing calculations in the ``environment`` model.

Each of the submpdules provide a description of how the link between the setup layer and the calculation layer is made, where settings can be defined, and where the resulting model can be extracted and used.

An overview of how to use the overall state propagation functionality can be found in the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation.html>`_

``estimation``
==============

.. toctree::
   :maxdepth: 2
   :caption: Modules

   /estimation/estimation_analysis
   /estimation/observable_models_setup
   /estimation/observations_setup
   /estimation/observable_models
   /estimation/observations
