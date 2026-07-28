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
   environment_data/sp3
   environment_data/spice
