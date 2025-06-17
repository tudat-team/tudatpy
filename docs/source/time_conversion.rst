.. _time_conversion:

``time_conversion``
===================
Conversions and computation on date and time.


This module provide a variety of functions to convert times and dates from Julian Dates, calendar dates and epochs.
Conversion between different time scales (UTC, TDB, TT, etc.) are also possible. Conversions between time scales make extensive use of the SOFA library.
Different helper functions are also included to ease the most common operation on dates and times.






Notes
-----
* Unless specified otherwise, the time used in Tudatpy is in seconds since J2000, noon on the 1st of January 2000 (e.g., this epoch defines :math:`t=0`)
* Tudat uses two different classes for date and time. One is the regular ``datetime`` from the Python ``datetime`` library. The other is a Tudat-native :class:`~tudatpy.astro.time_conversion.DateTime` class, which allows finer resolution for time definitions, and allows for easier conversion to time representations (seconds since epoch, Julian day, modified Julian day). You can convert between the two using the :meth:`~tudatpy.astro.time_conversion.DateTime.from_python_datetime` and :meth:`~tudatpy.astro.time_conversion.DateTime.to_python_datetime` methods.
* A number of conversion functions take the current Julian day or modified Julian day (as a ``float``) as input or output. This represents the number of days (86400 seconds) since noon January 1st 4713 BC, or midnight November 17th 1858 AD, respectively.
* A number of conversion functions take the seconds/days/... "since epoch" as input or output.




References
----------
*  Chapter 2 of: Kaplan, G. `United States Naval Observatory Circular No. 179 <https://www.usno.navy.mil/USNO/astronomical-applications/publications/Circular_179.pdf/view>`_, The IAU Resolutions on Astronomical Reference Systems, Time Scales, and Earth Rotation Models.
*  SOFA Documentation - `SOFA Time Scale and Calendar Tools <https://www.iausofa.org/sofa_ts_c.pdf>`_, disseminated by the International Astronomical Union
*  Secions 5.5.3 and Chapter 10 of the `IERS 2010 Conventions <https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36.pdf>`_, disseminated by the International Earth Rotation Service






Functions
---------
.. currentmodule:: tudatpy.astro.time_conversion

.. autosummary::

   date_time_components_to_epoch

   iso_string_to_epoch

   julian_day_to_python_datetime

   julian_day_to_seconds_since_epoch

   julian_day_to_modified_julian_day

   modified_julian_day_to_julian_day

   python_datetime_to_julian_day

   python_datetime_to_days_since_epoch

   seconds_since_epoch_to_julian_day

   seconds_since_epoch_to_julian_years_since_epoch

   seconds_since_epoch_to_julian_centuries_since_epoch

   calculate_seconds_in_current_julian_day

   is_leap_year

   get_days_in_month

   TCB_to_TDB

   TDB_to_TCB

   TCG_to_TT

   TT_to_TCG

   TAI_to_TT

   TT_to_TAI

   TT_to_TDB_approximate
   
   TT_to_TDB
   
   TDB_to_TT
   
   default_time_scale_converter


.. autofunction:: tudatpy.astro.time_conversion.date_time_components_to_epoch

.. autofunction:: tudatpy.astro.time_conversion.iso_string_to_epoch

.. autofunction:: tudatpy.astro.time_conversion.julian_day_to_python_datetime

.. autofunction:: tudatpy.astro.time_conversion.julian_day_to_seconds_since_epoch

.. autofunction:: tudatpy.astro.time_conversion.julian_day_to_modified_julian_day

.. autofunction:: tudatpy.astro.time_conversion.modified_julian_day_to_julian_day

.. autofunction:: tudatpy.astro.time_conversion.python_datetime_to_julian_day

.. autofunction:: tudatpy.astro.time_conversion.python_datetime_to_days_since_epoch

.. autofunction:: tudatpy.astro.time_conversion.seconds_since_epoch_to_julian_day

.. autofunction:: tudatpy.astro.time_conversion.seconds_since_epoch_to_julian_years_since_epoch

.. autofunction:: tudatpy.astro.time_conversion.seconds_since_epoch_to_julian_centuries_since_epoch

.. autofunction:: tudatpy.astro.time_conversion.calculate_seconds_in_current_julian_day

.. autofunction:: tudatpy.astro.time_conversion.is_leap_year

.. autofunction:: tudatpy.astro.time_conversion.get_days_in_month

.. autofunction:: tudatpy.astro.time_conversion.TCB_to_TDB

.. autofunction:: tudatpy.astro.time_conversion.TDB_to_TCB

.. autofunction:: tudatpy.astro.time_conversion.TCG_to_TT

.. autofunction:: tudatpy.astro.time_conversion.TT_to_TCG

.. autofunction:: tudatpy.astro.time_conversion.TAI_to_TT

.. autofunction:: tudatpy.astro.time_conversion.TT_to_TAI

.. autofunction:: tudatpy.astro.time_conversion.TT_to_TDB_approximate

.. autofunction:: tudatpy.astro.time_conversion.TT_to_TDB

.. autofunction:: tudatpy.astro.time_conversion.TDB_to_TT

.. autofunction:: tudatpy.astro.time_conversion.default_time_scale_converter


Enumerations
------------
.. currentmodule:: tudatpy.astro.time_conversion

.. autosummary::

   TimeScales



.. autoclass:: tudatpy.astro.time_conversion.TimeScales
   :members:





Classes
-------
.. currentmodule:: tudatpy.astro.time_conversion

.. autosummary::

   DateTime
   
   TimeScaleConverter

.. autoclass:: tudatpy.astro.time_conversion.DateTime
   :members:
   :special-members: __init__
   :exclude-members: iso_string, day_of_year, epoch, julian_day, modified_julian_day

.. autoclass:: tudatpy.astro.time_conversion.TimeScaleConverter
   :members:




