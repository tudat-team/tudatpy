.. _earth_orientation_data:

``earth_orientation``
=====================

Load Earth-orientation parameters (EOP) as source data before constructing an
Earth rotation model. This module preserves the supplied epochs and metadata;
interpolation and decisions about missing or out-of-range values belong to the
consuming model.

.. code-block:: python

   from tudatpy.data_input.environment_data import earth_orientation

   data = earth_orientation.read_earth_orientation_data(
       "eopc04.txt", earth_orientation.EarthOrientationFormat.c04)
   first_mjd = min(data.records)
   record = data.records[first_mjd]
   print(record.polar_motion)       # xp, yp in radians
   print(record.ut1_minus_utc)      # seconds
   print(record.celestial_pole_offsets)  # dX, dY in radians, or None

The supported formats are legacy 16-column C04, modern 21-column C04 with an
hour field and polar-motion rates, and fixed-width ``finals2000A`` Bulletin A.
The format argument selects C04 or finals2000A; C04 column counts are detected
within that choice. The modern C04 MJD already includes its hour field.
Bulletin B columns in finals2000A are not substituted for Bulletin A values.
Older dPsi/dEps celestial-pole conventions are rejected rather than interpreted
as dX/dY.

Record keys and ``utc_modified_julian_date`` are UTC modified Julian dates.
Angles and their uncertainties use radians, polar-motion rates use radians per
second, and UT1-minus-UTC and excess length of day use seconds. Missing optional
quantities remain ``None``. A supplied zero is retained as a measured value.
The separate ``polar_motion_status``, ``ut1_status`` and
``celestial_pole_status`` retain Bulletin A's ``I``/``P`` flags for interpolated
or predicted data. These flags are empty for C04, which does not supply them.

The data object retains the source name, original headers, format description
and frame/convention metadata. A terrestrial-frame declaration is preserved
from the C04 header when available; a more specific realization is not inferred
from the filename. Future date-only placeholders in finals2000A do not extend
the data coverage. The reader adds no endpoint padding.

Malformed records raise an exception with source and line information. Checks
include calendar/MJD consistency, duplicate or unsorted epochs, incomplete
coordinate pairs, truncated fixed-width records, nonfinite values and negative
uncertainties. Both file and text readers apply the same checks.

The C++ equivalents are ``tudat::input_output::EarthOrientationData`` and
``readEarthOrientationData`` in ``tudat/io/readEarthOrientationData.h``. The
separately linkable ``tudat_environment_data`` target uses C++17 and Eigen
and has no dynamics, estimation, interpolation or SOFA dependency.

.. currentmodule:: tudatpy.data_input.environment_data.earth_orientation

.. autoclass:: EarthOrientationFormat
   :members:

.. autoclass:: EarthOrientationData
   :members:

.. autoclass:: EarthOrientationRecord
   :members:

.. autofunction:: read_earth_orientation_data
