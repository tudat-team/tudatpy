.. _tracking_data_gaia:

``gaia``
========

This module loads Gaia solar-system astrometry, asteroid states, and covariance
data from the online archive or local parquet archives. :class:`GaiaAstrometry`
prepares TDB epochs, angular observations, and transit-correlated weights.
Convert it to generic tracking data with :meth:`GaiaAstrometry.to_tracking_data`.
Each Gaia transit becomes one tracking-data object and one observation set.
``read_gaia_astrometry_data`` reads an existing local parquet file or queries the
AIP FPR mirror using a source-ID lookup, then saves the raw response if an
archive path is supplied. It returns the filtered observations directly as
tracking data and supplementary data.
The returned supplementary data contains Gaia's geocentric state history and
is installed with
:func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`.
Requested photocenter and light-deflection corrections are stored as plain
settings and evaluated by
:func:`~tudatpy.estimation.observations.create_observation_dataset_from_tracking_data`
when ``apply_corrections=True``; source observations remain unchanged.

.. currentmodule:: tudatpy.data_input.tracking_data.gaia

.. autosummary::

   GaiaAstrometry
   read_gaia_astrometry_data
   download_gaia_archive
   gaia_object_catalog
   get_state_from_gaia_archive
   get_state_covariance_from_gaia_archive
   get_kepler_covariance_from_gaia_archive

.. autoclass:: GaiaAstrometry
   :members:

.. autofunction:: read_gaia_astrometry_data
.. autofunction:: download_gaia_archive
.. autofunction:: gaia_object_catalog
.. autofunction:: get_state_from_gaia_archive
.. autofunction:: get_state_covariance_from_gaia_archive
.. autofunction:: get_kepler_covariance_from_gaia_archive
