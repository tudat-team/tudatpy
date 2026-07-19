.. _tracking_data_obs_80_cols:

``obs_80_cols``
===============

.. automodule:: tudatpy.data_input.tracking_data.obs_80_cols

This submodule contains functionality to load tracking data from MPC 80-column
optical astrometry files. These files use the fixed-width astrometric
observation format described by the Minor Planet Center at
https://minorplanetcenter.net/iau/info/OpticalObs.html. The
:func:`read_80_column_data` function is the main interface for loading the data
and converting it to objects that Tudat can process further; see also
:ref:`tracking_data`. All other functionality in this module is
reserved for better understanding what data is being loaded, and in some cases
manipulating it, before it is processed into Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.obs_80_cols

.. autofunction:: read_80_column_data

Supporting API
--------------

The function below is used internally by :func:`read_80_column_data` to parse
the fixed-width MPC records into an astropy table. It does not create Tudat
tracking-data objects, but can be used to inspect the parsed optical astrometry
before conversion.

.. autofunction:: parse_80cols_file
