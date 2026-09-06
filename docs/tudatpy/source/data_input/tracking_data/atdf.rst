.. _tracking_data_atdf:

``atdf``
========

.. automodule:: tudatpy.data_input.tracking_data.atdf

This submodule contains functionality to load tracking data from ATDF/TRK-2-25
files. These NASA Deep Space Network (DSN) closed-loop archival tracking data
files are the earliest closed-loop DSN radio-science products, encoding
Doppler, range, and ramp records. The ``atdf2ascii`` Python package is used to
decode the binary ATDF files into intermediate ASCII tables. The
:func:`read_atdf_data` function is the main interface for loading the decoded
data and converting it to objects that Tudat can process further; see also
:ref:`tracking_data`.

.. currentmodule:: tudatpy.data_input.tracking_data.atdf

.. autofunction:: read_atdf_data

Supporting API
--------------

The class below exposes the lower-level ATDF processing controls used by
:func:`read_atdf_data`. It can be used to inspect or customize how the tables
decoded by ``atdf2ascii`` are converted, before the tracking-data and
supplementary-data objects are returned.

.. autoclass:: AtdfTrackingDataProcessor
   :members:
