.. _tracking_data_obs_80_cols:

``obs_80_cols``
===============

.. automodule:: tudatpy.data_input.tracking_data.obs_80_cols

This submodule loads optical astrometry, space-based astrometry and radar data
from MPC 80-column files. The supported fixed-width record layouts are
described by the Minor Planet Center for
`optical and space-based observations <https://minorplanetcenter.net/iau/info/OpticalObs.html>`_
and `radar observations <https://minorplanetcenter.net/iau/info/RadarObs.html>`_. The
:func:`read_80_column_data` function is the main interface for loading the data
and converting it to objects that Tudat can process further; see also
:ref:`tracking_data`.

.. currentmodule:: tudatpy.data_input.tracking_data.obs_80_cols

.. autofunction:: read_80_column_data

Supporting API
--------------

The function below is used internally by :func:`read_80_column_data` to parse
the fixed-width MPC records into an astropy table. It does not create Tudat
tracking-data objects, but can be used to inspect parsed optical and
space-based astrometry. Radar observations are available in the returned
table's metadata as a canonical radar table.

.. autofunction:: parse_80cols_file
