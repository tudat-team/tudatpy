.. _tracking_data_gaia:

``gaia``
========

This module loads Gaia solar-system astrometry, asteroid states, and covariance
data from the online archive or local parquet archives. :class:`GaiaAstrometry`
prepares TDB epochs, angular observations, and transit-correlated weights.
Convert it to generic tracking data with :meth:`GaiaAstrometry.to_tracking_data`,
or directly to an :class:`~tudatpy.estimation.observations.ObservationDataset`
with :meth:`GaiaAstrometry.to_observation_dataset`.

.. currentmodule:: tudatpy.data_input.tracking_data.gaia

.. autosummary::

   GaiaAstrometry
   generate_astrometry_parquet
   generate_asteroid_parquet
   gaia_object_catalog
   get_state_from_gaia_archive
   get_state_covariance_from_gaia_archive
   get_kepler_covariance_from_gaia_archive

.. autoclass:: GaiaAstrometry
   :members:

.. autofunction:: generate_astrometry_parquet
.. autofunction:: generate_asteroid_parquet
.. autofunction:: gaia_object_catalog
.. autofunction:: get_state_from_gaia_archive
.. autofunction:: get_state_covariance_from_gaia_archive
.. autofunction:: get_kepler_covariance_from_gaia_archive
