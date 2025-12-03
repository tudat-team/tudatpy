.. _processTrk234:

``processTrk234``
=================
Processing of NASA Deep Space Network (DSN) tracking data files (TRK-2-34 format).

This module provides tools to read, parse, and convert DSN radiometric tracking data
(Doppler, range, and ramp data) from TNF/TRK files into Tudat observation formats.


Classes
-------
.. currentmodule:: tudatpy.data.processTrk234

.. autosummary::

   Trk234Processor


Converter Classes
~~~~~~~~~~~~~~~~~

.. autosummary::

   converters.Converter
   converters.RadioBase
   converters.DerivedDopplerConverter
   converters.DerivedSraRangeConverter
   converters.RampConverter


Main Processor
~~~~~~~~~~~~~~

.. autoclass:: tudatpy.data.processTrk234.Trk234Processor
   :members:
   :special-members: __init__


Converters
~~~~~~~~~~

.. autoclass:: tudatpy.data.processTrk234.converters.Converter
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.data.processTrk234.converters.RadioBase
   :members:
   :special-members: __init__
   :show-inheritance:

.. autoclass:: tudatpy.data.processTrk234.converters.DerivedDopplerConverter
   :members:
   :special-members: __init__
   :show-inheritance:

.. autoclass:: tudatpy.data.processTrk234.converters.DerivedSraRangeConverter
   :members:
   :special-members: __init__
   :show-inheritance:

.. autoclass:: tudatpy.data.processTrk234.converters.RampConverter
   :members:
   :special-members: __init__
   :show-inheritance: