.. _tracking_data_udl:

udl
===

.. automodule:: tudatpy.data_input.tracking_data.udl

This submodule reads time- and frequency-difference-of-arrival (TDOA/FDOA)
JSON records exported from the Unified Data Library (UDL). The
:func:`read_utas_data` function converts one or more single-target files to the
common :ref:`tracking_data` containers. Use :class:`BatchUTAS` directly to
inspect grouped measurements or create receiving stations from the coordinates
in the source records.

.. currentmodule:: tudatpy.data_input.tracking_data.udl

.. autofunction:: read_utas_data

Supporting API
--------------

.. autoclass:: BatchUTAS
   :members:

.. autoclass:: StationPairObservations
   :members:

.. autoclass:: UTASMetadata
   :members:
