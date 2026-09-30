.. _tracking_data_atdf:

``atdf``
========

.. automodule:: tudatpy.data_input.tracking_data.atdf

This submodule contains functionality to load tracking data from ATDF/TRK-2-25
files, based on the ``atdf2ascii`` :cite:p:`verma2022PythonbasedToolConstructing` Python
package, which is used to decode the binary ATDF files into intermediate ASCII tables.
The :func:`read_atdf_data` function decodes the binary ATDF files and converts them to
Tudat-compatible tracking data objects; see also :ref:`tracking_data`.

Alternatively, the lower-level :class:`AtdfTrackingDataProcessor` class can be used to
decode the binary ATDF files into ASCII tables, using the :meth:`AtdfTrackingDataProcessor.convert_atdf_to_ascii` method,
and then convert the ASCII tables to tracking-data and supplementary-data objects, using the :meth:`AtdfTrackingDataProcessor.process_ascii_tables` method.

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
