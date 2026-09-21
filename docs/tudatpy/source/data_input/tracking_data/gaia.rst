.. _tracking_data_gaia:

``gaia``
========

This module loads Gaia solar-system astrometry from the online archive or a
local parquet archive. :class:`GaiaAstrometry` prepares TDB epochs, angular
observations, source metadata, and transit-correlated weights. Convert it to
generic tracking data with :meth:`GaiaAstrometry.to_tracking_data`, or directly
to an :class:`~tudatpy.estimation.observations.ObservationDataset` with
:meth:`GaiaAstrometry.to_observation_dataset`.

.. currentmodule:: tudatpy.data_input.tracking_data.gaia

.. autosummary::

   GaiaAstrometry
   generate_astrometry_parquet

.. autoclass:: GaiaAstrometry
   :members:

.. autofunction:: generate_astrometry_parquet
