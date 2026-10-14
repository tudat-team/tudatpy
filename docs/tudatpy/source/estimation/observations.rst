.. _observations:

``observations``
=====================

This module contains the objects in Tudat that are used to store and process observations. The :class:`~tudatpy.estimation.observations.ObservationCollection` is the object for which a complete set of observations (e.g. all observations to be used in a single estimation) are stored. This object is composed of a list of :class:`~tudatpy.estimation.observations.SingleObservationSet` objects, each of which stores observations of a common type, link ends, etc. Dealing with observations in Tudat is discussed extensively on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation.html>`_. In addition to the observations, these objects also store observation residuals and dependent variables. These are normally computed during an estimation (see :ref:`estimation_analysis`), but can also be computed outside an estimation loop using the function :func:`~tudatpy.estimation.observations.compute_residuals_and_dependent_variables`

In addition to the top-level classes, this module contains functionality to manipulate objects of these classes (:func:`~tudatpy.estimation.observations.merge_observation_collections`, :func:`~tudatpy.estimation.observations.split_observation_collection`, :func:`~tudatpy.estimation.observations.create_new_observation_collection`). For a more detailed discussion on the use of these functions, as well as the filtering of observations using :func:`~tudatpy.estimation.observations.create_filtered_observation_collection`, :func:`~tudatpy.estimation.observations.create_filtered_observation_set` (e.g. removing observations outside a certain time range, with a residual higher than a given value, etc.), see our `corresponding user guide section <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/observation-collection-manipulation/processing-observations.html#filtering-observations>`_.




Selecting observations of one type
----------------------------------

Use an observation parser with the collection methods to retrieve, weight or
filter a subset of a mixed observation collection. For example, given a
collection named ``collection`` whose residuals have already been computed:

.. code-block:: python

   from tudatpy.estimation.observable_models_setup import model_settings
   from tudatpy.estimation.observations import observations_processing

   pixels = observations_processing.observation_parser(
       model_settings.pixel_coordinates_type
   )
   pixel_residuals = collection.get_concatenated_residuals(pixels)

   sigma_pixels = 2.0
   collection.set_constant_weight(1.0 / sigma_pixels**2, pixels)

   residual_filter = observations_processing.observation_filter(
       observations_processing.ObservationFilterType.residual_filtering, 10.0
   )
   collection.filter_observations(residual_filter, pixels)

The same parser can be passed to ``get_concatenated_observations`` and
``get_concatenated_weights``. Component-wise uncertainties can be supplied by
passing a vector of inverse variances to ``set_constant_weight``.


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

   create_single_observation_set

   create_filtered_observation_set

   split_observation_set

   create_filtered_observation_collection

   merge_observation_collections

   split_observation_collection

   create_new_observation_collection

   create_observation_collection_from_tracking_data

   create_observation_collection_from_arrays

   create_single_type_observation_collection_from_arrays

   simulate_observations

   simulate_pseudo_observations

   set_tracking_supplementary_data_in_bodies


.. autofunction:: tudatpy.estimation.observations.compute_residuals_and_dependent_variables

.. autofunction:: tudatpy.estimation.observations.create_single_observation_set

.. autofunction:: tudatpy.estimation.observations.create_filtered_observation_set

.. autofunction:: tudatpy.estimation.observations.split_observation_set

.. autofunction:: tudatpy.estimation.observations.merge_observation_collections

.. autofunction:: tudatpy.estimation.observations.create_filtered_observation_collection

.. autofunction:: tudatpy.estimation.observations.split_observation_collection

.. autofunction:: tudatpy.estimation.observations.create_new_observation_collection

.. autofunction:: tudatpy.estimation.observations.create_observation_collection_from_tracking_data

.. autofunction:: tudatpy.estimation.observations.create_observation_collection_from_arrays

.. autofunction:: tudatpy.estimation.observations.create_single_type_observation_collection_from_arrays

.. autofunction:: tudatpy.estimation.observations.simulate_observations

.. autofunction:: tudatpy.estimation.observations.simulate_pseudo_observations

.. autofunction:: tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies


Classes
-------
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   SingleObservationSet

   ObservationCollection

.. autoclass:: tudatpy.estimation.observations.SingleObservationSet
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.estimation.observations.ObservationCollection
   :members:
   :special-members: __init__
