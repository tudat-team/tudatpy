.. _observations:

``observations``
=====================

This module contains the objects in Tudat that are used to store and process observations. The :class:`~tudatpy.estimation.observations.ObservationCollection` is the object for which a complete set of observations (e.g. all observations to be used in a single estimation) are stored. This object is composed of a list of :class:`~tudatpy.estimation.observations.SingleObservationSet` objects, each of which stores observations of a common type, link ends, etc. Dealing with observations in Tudat is discussed extensively on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation.html>`_. In addition to the observations, these objects also store observation residuals and dependent variables. These are normally computed during an estimation (see :ref:`estimation_analysis`), but can also be computed outside an estimation loop using the function :func:`~tudatpy.estimation.observations.compute_residuals_and_dependent_variables`

In addition to the top-level classes, this module contains functionality to manipulate objects of these classes (:func:`~tudatpy.estimation.observations.merge_observation_collections`, :func:`~tudatpy.estimation.observations.split_observation_collection`, :func:`~tudatpy.estimation.observations.create_new_observation_collection`). For a more detailed discussion on the use of these functions, as well as the filtering of observations using :func:`~tudatpy.estimation.observations.create_filtered_observation_collection`, :func:`~tudatpy.estimation.observations.create_filtered_observation_set` (e.g. removing observations outside a certaint time range, with a residual higher than a given value, *etc.*, see our `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/observation-collection-manipulation/processing-observations.html#filtering-observations>`




.. toctree::
   :maxdepth: 2
   :caption: Modules

   /estimation/observations/observations_geometry
   /estimation/observations/observations_processing


Functions
---------
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   compute_residuals_and_dependent_variables

   single_observation_set

   create_single_observation_set

   filter_observations

   create_filtered_observation_set

   split_observation_set

   create_filtered_observation_collection

   merge_observation_collections

   split_observation_collection

   create_new_observation_collection


.. autofunction:: tudatpy.estimation.observations.compute_residuals_and_dependent_variables

.. autofunction:: tudatpy.estimation.observations.single_observation_set

.. autofunction:: tudatpy.estimation.observations.create_single_observation_set

.. autofunction:: tudatpy.estimation.observations.filter_observations

.. autofunction:: tudatpy.estimation.observations.create_filtered_observation_set

.. autofunction:: tudatpy.estimation.observations.split_observation_set

.. autofunction:: tudatpy.estimation.observations.merge_observation_collections

.. autofunction:: tudatpy.estimation.observations.create_filtered_observation_collection

.. autofunction:: tudatpy.estimation.observations.split_observation_collection

.. autofunction:: tudatpy.estimation.observations.create_new_observation_collection


Classes
-------
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   SingleObservationSet

   ObservationCollection

.. autoclass:: tudatpy.estimation.observations.SingleObservationSet
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationCollection
   :members: