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

   environment_data/atmospheric_loading
   environment_data/station_loading
   environment_data/coma
   environment_data/discos
   environment_data/earth_orientation
   environment_data/gravity_field
   environment_data/horizons
   environment_data/missions
   environment_data/sbdb
   environment_data/ocean_loading
   environment_data/space_weather
   environment_data/station_history
   environment_data/spacetrack
   environment_data/spice
   environment_data/terrestrial_frames

SP3 ephemeris products
---------------------

.. currentmodule:: tudatpy.data_input.environment_data

The SP3 reader loads uncompressed SP3c/d files as raw ephemeris records. Positions
and velocities are converted once to metres and metres per second. Clock
corrections are seconds and clock rates are seconds per second; missing clock
values are NaN. Missing positions/velocities and predicted-orbit flags remain
explicit. The reader neither interpolates an ephemeris nor converts its frame
or time system.

Some ILRSA combined products omit the literal ``EOF`` record. A physical end of
file is accepted when the declared epoch count and all final satellite records
are complete; missing samples and incomplete position/velocity pairs still fail.

Each record's :class:`~tudatpy.astro.time_representation.Time` represents the
coordinate printed in the file's declared time system. In particular, an ILRS
SP3 declaring UTC must not be treated as a GPS-time product. Frame labels and
provider comments must be interpreted before constructing an environment model;
a frame label alone does not specify an available reference-frame conversion.

.. autosummary::

   read_sp3_file
   Sp3Data
   Sp3Record

.. autofunction:: read_sp3_file

.. autoclass:: Sp3Data
   :members:

.. autoclass:: Sp3Record
   :members:
