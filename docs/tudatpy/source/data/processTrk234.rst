.. _processTrk234:

``processTrk234``
=================
Processing of NASA Deep Space Network (DSN) tracking data files (TRK-2-34 format).

This module provides tools to read, parse, and convert DSN radiometric tracking data
(Doppler, range, and ramp data) from TNF/TRK files into Tudat observation formats.


Enums
-----
.. currentmodule:: tudatpy.data.processTrk234

.. autosummary::

   converters.OpenRampHandling

.. autoclass:: tudatpy.data.processTrk234.converters.OpenRampHandling
   :members:

Classes
-------
.. currentmodule:: tudatpy.data.processTrk234

.. autosummary::

   Trk234Processor
   converters.Converter
   converters.RadioBase
   converters.DerivedDopplerConverter
   converters.DerivedSraRangeConverter
   converters.RampConverter


Main Processor
~~~~~~~~~~~~~~

.. autoclass:: tudatpy.data.processTrk234.Trk234Processor
   :members:
   :exclude-members: process_observation_collection
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
