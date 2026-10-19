.. _environment_data:

``environment_data``
====================

This module contains interfaces for loading and querying external data that can
be used to define environment models in Tudat. It includes interfaces for SPICE
kernels, JPL Horizons, JPL's Small-Body Database, Space-Track, DISCOS,
space-weather data, comet-coma data, and mission-specific environment products.

The functions and classes in this module retrieve or load data. The conversion
of these data into Tudat environment models is handled by the corresponding
functions in :ref:`environment_setup`.

.. toctree::
   :maxdepth: 2
   :caption: Modules

   environment_data/coma
   environment_data/discos
   environment_data/horizons
   environment_data/missions
   environment_data/sbdb
   environment_data/space_weather
   environment_data/spacetrack
   environment_data/spice

Functions
---------
.. currentmodule:: tudatpy.data_input.environment_data

.. autosummary::

   read_dsn_amc_file

.. autofunction:: tudatpy.data_input.environment_data.read_dsn_amc_file

Classes
-------
.. currentmodule:: tudatpy.data_input.environment_data

.. autosummary::

   DsnAmcRecord

.. autoclass:: tudatpy.data_input.environment_data.DsnAmcRecord
   :members:
