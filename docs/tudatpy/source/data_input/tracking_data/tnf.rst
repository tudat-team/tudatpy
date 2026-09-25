.. _tracking_data_tnf:

``tnf``
=======

.. automodule:: tudatpy.data_input.tracking_data.tnf

This submodule contains functionality to load tracking data from TNF/TRK-2-34
files. These NASA Deep Space Network (DSN) Tracking and Navigation Files encode
radiometric tracking data, including Doppler, range, and ramp records.
The format is described in `820-013, TRK-2-34 DSN Tracking System Data Archival
Format <https://pds-geosciences.wustl.edu/radiosciencedocs/urn-nasa-pds-radiosci_documentation/dsn_trk-2-34/dsn_trk-2-34.2021-06-03.pdf>`_.
The `pytrk234 <https://github.com/NASA-PDS/PyTrk234/>`_ Python package is used heavily for parsing of the binary files and extraction of the raw data records.

The :func:`read_tnf_data` function is the main interface for loading
the decoded data and converting it to objects that Tudat can process further;
see also :ref:`tracking_data`. All other functionality in this module is
reserved for better understanding what data is being loaded, and in some cases
manipulating it, before it is processed into Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.tnf

.. autofunction:: read_tnf_data

Supporting API
--------------

The classes below expose the lower-level TNF processing controls used by
:func:`read_tnf_data`. They can be used to inspect or customize how records
decoded by ``pytrk234`` are converted, in particular how open ramp intervals are
handled before the tracking-data and supplementary-data objects are returned.

.. autoclass:: OpenRampHandling
   :members:

.. autoclass:: TnfTrackingDataProcessor
   :members:
