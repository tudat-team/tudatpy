.. _tracking_data_jpl_radar:

``jpl_radar``
=============

.. automodule:: tudatpy.data_input.tracking_data.jpl_radar

This submodule retrieves small-body radar astrometry from the
`JPL Small-Body Radar Astrometry API <https://ssd-api.jpl.nasa.gov/doc/sb_radar.html>`_.
Delay and Doppler measurements are converted to Tudat tracking-data objects,
and known JPL radar station identifiers are mapped to their MPC observatory
codes. The :func:`read_jpl_radar_data` function provides the direct conversion
interface.

.. currentmodule:: tudatpy.data_input.tracking_data.jpl_radar

.. autofunction:: read_jpl_radar_data

Supporting API
--------------

The query class and catalogue function below expose the canonical radar table,
raw API records and station coordinates before conversion to Tudat tracking
data.

.. autofunction:: get_available_radar_targets

.. autoclass:: JPLRadarQuery
   :members:
