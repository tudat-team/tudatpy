.. _tracking_data_odf:

``odf``
=======

.. automodule:: tudatpy.data_input.tracking_data.odf

This submodule contains functionality to load tracking data from ODF files.
ODF files are binary files produced by the Deep Space Network (DSN) and encode
radiometric spacecraft tracking data used for orbit and parameter estimation.
The format is described in ``820-013, TRK-2-18 Tracking System Interfaces Orbit
Data File Interface, Revision E, 2008, JPL/DSN``. The :func:`read_odf_data`
function is the main interface for loading the data and converting it to objects
that Tudat can process further; see also :ref:`tracking_data`. All
other functionality in this module is reserved for better understanding what
data is being loaded, and in some cases manipulating it, before it is processed
into Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.odf

.. autofunction:: read_odf_data

Supporting API
--------------

The functions and classes below expose the raw contents of ODF binary files in
dedicated containers. They are used internally by :func:`read_odf_data` and do
not perform the full conversion to Tudat tracking-data objects, but can be used
to inspect the parsed ODF contents if needed.

.. autofunction:: read_raw_odf_file_contents

.. autoclass:: RawOdfFileContents
   :members:

.. autoclass:: OdfCommonDataBlock
   :members:

.. autoclass:: OdfDataBlock
   :members:

.. autoclass:: OdfDataSpecificBlock
   :members:

.. autoclass:: OdfDopplerDataBlock
   :members:

.. autoclass:: OdfRampBlock
   :members:

.. autoclass:: OdfDataType
