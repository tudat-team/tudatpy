.. _tracking_data_gaia:

``gaia``
========

This module contains functionality to load Gaia astrometric observations of asteroids, as well as
asteroid states and covariances, from the Gaia archive or from local parquet files. The astrometric
observations of an asteroid are loaded with
:func:`~tudatpy.data_input.tracking_data.gaia.load_gaia_astrometry`, which reads a local parquet
file if it exists, and otherwise retrieves the observations from the online archive (saving them to
the parquet file if a path is provided). The observations are stored in a
:class:`~tudatpy.data_input.tracking_data.gaia.GaiaAstrometry` object, with observation times in
TDB, and with weights that include the correlations between the observations of each transit. They
are converted to tracking data with
:meth:`~tudatpy.data_input.tracking_data.gaia.GaiaAstrometry.to_tracking_data`, which creates one
:class:`~tudatpy.data_input.tracking_data.TrackingData` object (and consequently one observation
set) per transit, together with supplementary data containing the geocentric state history of Gaia.
This supplementary data is set in the bodies with
:func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`. The photocenter
and light deflection corrections requested in
:meth:`~tudatpy.data_input.tracking_data.gaia.GaiaAstrometry.to_tracking_data` are computed when
creating the observation dataset with
:func:`~tudatpy.estimation.observations.create_observation_dataset_from_tracking_data` (with
``apply_corrections=True``); the observations stored in the tracking data are not modified.

.. currentmodule:: tudatpy.data_input.tracking_data.gaia

.. autosummary::

   GaiaAstrometry
   load_gaia_astrometry
   generate_astrometry_parquet
   generate_asteroid_parquet
   gaia_object_catalog
   get_state_from_gaia_archive
   get_state_covariance_from_gaia_archive
   get_kepler_covariance_from_gaia_archive

.. autoclass:: GaiaAstrometry
   :members:

.. autofunction:: load_gaia_astrometry
.. autofunction:: generate_astrometry_parquet
.. autofunction:: generate_asteroid_parquet
.. autofunction:: gaia_object_catalog
.. autofunction:: get_state_from_gaia_archive
.. autofunction:: get_state_covariance_from_gaia_archive
.. autofunction:: get_kepler_covariance_from_gaia_archive
