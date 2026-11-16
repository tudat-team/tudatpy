.. _tracking_data_generic_text_file:

``generic_text_file``
=====================

.. automodule:: tudatpy.data_input.tracking_data.generic_text_file

This submodule contains functionality to load tracking data from generic text
files with user-defined tracking-data columns. The
:func:`read_generic_text_data` function is the main interface for loading the
data and converting it to parsed tracking text-file contents that Tudat can
process further; see also :ref:`tracking_data`. All other
functionality in this module is reserved for better understanding what data is
being loaded, and in some cases manipulating it, before it is processed into
Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.generic_text_file

.. autofunction:: read_generic_text_data

Supporting API
--------------

The classes below define the parsed generic text-file contents and the column
or filter identifiers used while reading these files. They can be used to
inspect the loaded columns and metadata before any later conversion step.

.. autoclass:: TrackingTxtFileContents
   :members:

.. autoclass:: TrackingDataType

.. autoclass:: TrackingTxtFileReadFilterType
