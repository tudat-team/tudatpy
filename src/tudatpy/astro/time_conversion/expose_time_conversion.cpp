/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_time_conversion.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/basic_astro/timeConversions.h>
#include <tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h>
#include <tudat/math/basic/mathematicalConstants.h>
#include <tudat/basics/deprecationWarnings.h>
#include <pybind11/operators.h>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <chrono>
#include <ctime>

#include "scalarTypes.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/basic_astro/physicalConstants.h"

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tsi = tudat::sofa_interface;
namespace pc = tudat::physical_constants;
namespace teo = tudat::earth_orientation;
namespace tutil = tudat::utilities;

namespace tudat
{

namespace earth_orientation
{
std::shared_ptr< TerrestrialTimeScaleConverter > createDefaultTimeConverterPy( )
{
    return createDefaultTimeConverter( );
}

}  // namespace earth_orientation

namespace basic_astrodynamics
{

}  // namespace basic_astrodynamics

}  // namespace tudat

// Convert from Gregorian date to time_point (Python datetime). Only
// year/month/day, no time.
std::chrono::system_clock::time_point dateTimeToTimePoint( const tba::DateTime& dateTime )
{
    tutil::printDeprecationWarning( "datetime_to_python", "DateTime.to_python_datetime" );
    return dateTime.timePoint( );
}

// Convert Julian day to calendar date. This code ensures that the value
// returned is a time_point (Python datetime).
std::chrono::system_clock::time_point convertJulianDayToCalendarDatePy( const double julianDay )
{
    tutil::printDeprecationWarning( "julian_day_to_calendar_date", "DateTime.from_julian_day(...).to_python_datetime()" );

    tba::DateTime dateTime = tba::DateTime::fromTime< double >( tudat::timeFromJulianDay< double >( julianDay ) );

    return dateTime.timePoint( );
}

// Convert calendar date to Julian day since a given epoch. This code allows for
// the calendar date to be a time_point (Python datetime).
template< typename TimeScalarType = double >
TimeScalarType convertCalendarDateToJulianDayPy( const std::chrono::system_clock::time_point calendarDate )
{
    tutil::printDeprecationWarning( "calendar_date_to_julian_day", "DateTime.from_python_datetime(...).to_julian_day()" );

    tba::DateTime dateTime = tba::DateTime::fromTimePoint( calendarDate );
    return dateTime.julianDay< TimeScalarType >( );
}

template< typename TimeScalarType = double >
TimeScalarType convertCalendarDateToJulianDaySinceEpochPy(
        const std::chrono::system_clock::time_point calendarDate,
        const TimeScalarType epochSinceJulianDayZero = tba::getJulianDayOnJ2000< TimeScalarType >( ) )
{
    tutil::printDeprecationWarning( "calendar_date_to_days_since_epoch",
                                    "DateTime.from_python_datetime(...).to_days_since_reference_julian_day()" );
    tba::DateTime dateTime = tba::DateTime::fromTimePoint( calendarDate );
    return dateTime.julianDay< TimeScalarType >( ) - epochSinceJulianDayZero;
}

namespace tudatpy
{

namespace astro
{
namespace time_conversion
{

void expose_time_conversion( py::module& m )
{

    py::class_< tudat::Time >( m, "Time", R"doc(No documentation found.)doc" )
            .def( py::init< const int, const long double >( ),
                  py::arg( "full_periods" ),
                  py::arg( "seconds_into_full_period" ) )
            .def( py::init< const double >( ), py::arg( "seconds_since_j2000" ) )
            .def( "to_float",
                  &tudat::Time::getSeconds< double >,
                  R"doc(No documentation found.)doc" )
            .def( py::self + py::self )
            .def( py::self + double( ) )
            .def( double( ) + py::self )
            .def( py::self += py::self )
            .def( py::self += double( ) )
            .def( py::self - py::self )
            .def( py::self - double( ) )
            .def( py::self -= py::self )
            .def( py::self -= double( ) )
            .def( double( ) - py::self )
            .def( py::self * double( ) )
            .def( double( ) * py::self )
            .def( py::self *= double( ) )
            .def( py::self / double( ) )
            .def( py::self /= double( ) )
            .def( py::self == py::self )
            .def( double( ) == py::self )
            .def( py::self == double( ) )
            .def( py::self != py::self )
            .def( py::self != double( ) )
            .def( double( ) != py::self )
            .def( py::self < py::self )
            .def( py::self < double( ) )
            .def( double( ) < py::self )
            .def( py::self > py::self )
            .def( py::self > double( ) )
            .def( double( ) > py::self )
            .def( py::self <= py::self )
            .def( py::self <= double( ) )
            .def( double( ) <= py::self )
            .def( py::self >= py::self )
            .def( double( ) >= py::self )
            .def( py::self >= double( ) );

    //    m.attr("default_time_converter") =
    //    tudat::earth_orientation::defaultTimeConverter;

    py::enum_< tba::TimeScales >( m,
                                  "TimeScales",
                                  R"doc(

 Enumeration of available time scales between which the :class:`~TimeScaleConverter` can automaticaly convert.

 )doc" )
            .value( "tai_scale", tba::tai_scale, R"doc(
 )doc" )
            .value( "tt_scale", tba::tt_scale, R"doc(
 )doc" )
            .value( "tdb_scale", tba::tdb_scale, R"doc(
 )doc" )
            .value( "utc_scale", tba::utc_scale, R"doc(
 )doc" )
            .value( "ut1_scale", tba::ut1_scale, R"doc(
 )doc" )
            .export_values( );

    py::class_< teo::TerrestrialTimeScaleConverter, std::shared_ptr< teo::TerrestrialTimeScaleConverter > >( m,
                                                                                                             "TimeScaleConverter",
                                                                                                             R"doc(

 Class to convert between different time scales (TAI, TT, TDB, UTC, UT1)

 Class to convert between different time scales (TAI, TT, TDB, UTC, UT1), as per algorithms described in (for instance) IERS 2010 Conventions and USNO circular no. 179.
 The algorithms used for the conversion are (where there is any choice in models):

 * Conversion between TDB and TT uses the SOFA function ``iauDtdb`` function (equivalently to :func:`TT_to_TDB` and :func:`TDB_to_TT`)
 * Conversion between UTC and UT1 applies (semi-)diurnal and daily measured variations depending on settings during object creation (typically as per :func:`~default_time_scale_converter`)
 * Leap seconds as per the latest SOFA package

 )doc" )
            .def( "convert_time",
                  &teo::TerrestrialTimeScaleConverter::getCurrentTime< TIME_TYPE >,
                  py::arg( "input_scale" ),
                  py::arg( "output_scale" ),
                  py::arg( "input_value" ),
                  py::arg( "earth_fixed_position" ) = Eigen::Vector3d::Zero( ) )
            .def( "get_time_difference",
                  &teo::TerrestrialTimeScaleConverter::getCurrentTimeDifference< double >,
                  py::arg( "input_scale" ),
                  py::arg( "output_scale" ),
                  py::arg( "input_value" ),
                  py::arg( "earth_fixed_position" ) = Eigen::Vector3d::Zero( ) );

    py::class_< tba::DateTime >( m, "DateTime", R"doc(

 Class to store a calendar date and time of day, with high resolution.

 Class to store a calendar date and time of day, with high resolution compared to Python datetime.datetime. This class
 stores the seconds as a ``long double`` variable in the C++ implementation, corresponding to about
 16 or 19 digits of precision (depending on the compiler used). In either case, this will be sufficient for sub-femtosecond
 resolution. In addition, this class allows easy conversion to typical time representations in astrodynamics (seconds since J2000,
 Julian day, and modified Julian day).





 )doc" )
            .def( py::init< const int, const int, const int, const int, const int, const long double >( ),
                  py::arg( "year" ),
                  py::arg( "month" ),
                  py::arg( "day" ),
                  py::arg( "hour" ) = 12,
                  py::arg( "minute" ) = 0,
                  py::arg( "seconds" ) = 0.0L )
            .def( "__str__", []( tba::DateTime& datetime ) { return datetime.isoString( ); } )
            .def( "__repr__",
                  []( const tba::DateTime& datetime ) {
                      return "DateTime(" + std::to_string( datetime.getYear( ) ) + ", " + std::to_string( datetime.getMonth( ) ) + ", " +
                              std::to_string( datetime.getDay( ) ) + ", " + std::to_string( datetime.getHour( ) ) + ", " +
                              std::to_string( datetime.getMinute( ) ) + ", " + std::to_string( datetime.getSeconds( ) ) + ")";
                  } )
            .def_property( "year", &tba::DateTime::getYear, &tba::DateTime::setYear, R"doc(

 Calendar year


 :type: int
 )doc" )
            .def_property( "month", &tba::DateTime::getMonth, &tba::DateTime::setMonth, R"doc(

 Calendar month (value must be 1-12)


 :type: int
 )doc" )
            .def_property( "day", &tba::DateTime::getDay, &tba::DateTime::setDay, R"doc(

 Calendar day in current month, value must be larger than 0, and smaller or equal to the number of days in the month


 :type: int
 )doc" )
            .def_property( "hour", &tba::DateTime::getHour, &tba::DateTime::setHour, R"doc(

 Full hours into the current day (value must be 0-23)


 :type: int
 )doc" )
            .def_property( "minute", &tba::DateTime::getMinute, &tba::DateTime::setMinute, R"doc(

 Full minutes into the current hour (value must be 0-59)


 :type: int
 )doc" )
            .def_property( "seconds", &tba::DateTime::getSeconds, &tba::DateTime::setSeconds, R"doc(

 Number of seconds into the current minute. Note that this value is stored as ``long double`` in Tudat, which may be 64-bit or 80-bit (16 or 19 digits) depending on the compiler used.


 :type: float
 )doc" )
            .def( "iso_string",
                  &tba::DateTime::isoString,
                  py::arg( "add_T" ) = false,
                  py::arg( "number_of_digits_seconds" ) = 15,
                  R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_iso_string` instead.

 Function to get the ISO-compatible string.


 Function to get the current date and time as an ISO-compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..") where the seconds may be provided with any number of digits. The 'T' entry separating the date from the time may be omitted by setting the ``add_T`` parameter to false


 Parameters
 ----------
 add_T : bool
    Boolean denoting whether to use a 'T' or a blank space to separate the date from the time
 number_of_digits_seconds : int, default = 15
    Number of digits to use after the decimal separator (trailing zeros will be truncated)

 Returns
 -------
 str
     ISO-compatible string representing the date and time





 )doc" )
            .def( "to_iso_string",
                  &tba::DateTime::isoString,
                  py::arg( "add_T" ) = false,
                  py::arg( "number_of_digits_seconds" ) = 15,
                  R"doc(

 Function to get the ISO-compatible string.


 Function to get the current date and time as an ISO-compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..") where the seconds may be provided with any number of digits. The 'T' entry separating the date from the time may be omitted by setting the ``add_T`` parameter to false


 Parameters
 ----------
 add_T : bool
 Boolean denoting whether to use a 'T' or a blank space to separate the date from the time

 number_of_digits_seconds : int, default = 15
 Number of digits to use after the decimal separator (trailing zeros will be truncated)

 Returns
 -------
 str
     ISO-compatible string representing the date and time





 )doc" )
            .def( "day_of_year",
                  &tba::DateTime::dayOfYear,
                  R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_day_of_year` instead.

 Function to get the day number in the current year


 Returns
 -------
 int
     Day number in the current year





 )doc" )
            .def( "to_day_of_year",
                  &tba::DateTime::dayOfYear,
                  R"doc(

 Function to get the day number in the current year


 Returns
 -------
 int
     Day number in the current year





 )doc" )
            .def( "epoch",
                  &tba::DateTime::epoch< TIME_TYPE >,
                  R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_epoch` instead.

 Function to get the epoch in seconds since J2000 for the current date and time


 Returns
 -------
 float
     Current epoch in seconds since J2000





 )doc" )
            .def( "to_epoch",
                  &tba::DateTime::epoch< TIME_TYPE >,
                  R"doc(

 Function to get the epoch in seconds since J2000 for the current date and time


 Returns
 -------
 float
     Current epoch in seconds since J2000





 )doc" )
            .def( "julian_day",
                  &tba::DateTime::julianDay< double >,
                  R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_julian_day` instead.

 Function to get the epoch as Julian day for the current date and time


 Returns
 -------
 float
     Current Julian day





 )doc" )
            .def( "to_julian_day",
                  &tba::DateTime::julianDay< double >,
                  R"doc(

 Function to get the epoch as Julian day for the current date and time


 Returns
 -------
 float
     Current Julian day





 )doc" )
            .def( "modified_julian_day",
                  &tba::DateTime::modifiedJulianDay< double >,
                  R"doc(
 
 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_modified_julian_day` instead.

 Function to get the epoch as modified Julian day for the current date and time


 Returns
 -------
 float
     Current modified Julian day





 )doc" )
            .def( "to_modified_julian_day",
                  &tba::DateTime::modifiedJulianDay< double >,
                  R"doc(

 Function to get the epoch as modified Julian day for the current date and time


 Returns
 -------
 float
     Current modified Julian day





 )doc" )
            .def_static( "from_python_datetime", &tba::DateTime::fromTimePoint, py::arg( "datetime" ), R"doc(
            
Function to convert a Python `datetime.datetime` object to a Tudat :class:`DateTime` object. The Tudat-native alternative has the advantage of providing sub-femtosecond resolution, as opposed to the microsecond resolution of the Python version.

.. warning::

    This function uses the C++ `std::chrono` library, which is limited in the time range it can represent. If the range is exceeded, the conversion will overflow and **NOT** throw an exception.

    The exact range is platform-dependent. On Windows, dates between 1970-01-01 and 3000-12-31 are allowed, on MacOS dates after 1900-01-01 are allowed and on Linux dates between 1678-01-01 and 2261-12-31 are allowed.

Parameters
----------
datetime : datetime.datetime
    Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
Returns
-------
DateTime
    DateTime object defined in Tudat

Examples
--------
In this example, the Tudat DateTime object is constructed from python native datetime object.

.. code-block:: python

    from datetime import datetime
    from tudatpy.astro.time_conversion import DateTime

    python_datetime = datetime(2025, 1, 1, 0, 0, 0)

    dt = DateTime.from_python_datetime(python_datetime)
    print(dt) # prints 2025-01-01 00:00:00.000000000000000

                        )doc" )
            .def_static( "from_year_and_day_of_year",
                         &tba::DateTime::fromYearAndDaysInYear,
                         py::arg( "year" ),
                         py::arg( "day_of_year" ),
                         R"doc(
                         
Create the Tudat :class:`DateTime` from the year and the number of days in the year.

Parameters
----------
year : int
    Calendar year.
day_of_year : int
    Number of days that have passed in the year.
Returns
-------
DateTime
    Corresponding calendar date as a :class:`DateTime` object. Note: the hours, minutes and seconds in the object are set to 0 when calling this function.

Examples
--------
In this example, the calendar date corresponding to when 122 days have passed in 2020 is computed.

.. code-block:: python

    # Compute the calendar date when 122 days have passed in 2020
    currentDate = time_conversion.DateTime.from_year_and_day_of_year(2020, 122)
    # Print the converted output
    print(currentDate)  # prints (2020, 5, 2, 0, 0)
                         
                         )doc" )
            .def_static( "from_iso_string", &tba::DateTime::fromIsoString, py::arg( "iso_time" ), R"doc(
            
 Creates a Tudat-native :class:`DateTime` object from an ISO datetime string.

 Parameters
 ----------
 iso_datetime : str
     Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)

 Returns
 -------
 DateTime
     Tudat ``DateTime`` object.

 Examples
 --------
 In this example, the datetime is constructed from the iso string.
 
 .. code-block:: python
 
     from tudatpy.astro.time_conversion import DateTime

     dt = DateTime.from_iso_string("2025-01-01T00:00:00.000")
     print(dt) # prints 2025-01-01 00:00:00.000000000000000
                         
                         )doc" )
            .def_static( "from_epoch",
                         &tba::DateTime::fromTime< TIME_TYPE >,
                         py::arg( "epoch" ),
                         R"doc(

 Creates a Tudat-native :class:`DateTime` object from the seconds since J2000.

 Parameters
 ----------
 epoch : float
     Seconds since J2000

 Returns
 -------
 DateTime
     Tudat ``DateTime`` object.

 Examples
 --------
 In this example, the datetime is constructed from an epoch in seconds since J2000.
 
 .. code-block:: python
 
     from tudatpy.astro.time_conversion import DateTime

     epoch_et = 788961600.0

     dt = DateTime.from_epoch(epoch_et)
     print(dt) # prints 2025-01-01 00:00:00.000000000000000
                         
                         )doc" )
            .def_static( "from_julian_day",
                         &tba::DateTime::fromJulianDay,
                         py::arg( "julian_day" ),
                         R"doc(

 Creates a Tudat-native :class:`DateTime` object from a Julian day.

 Parameters
 ----------
 julian_day : float
     Date in Julian days since January 1st 4713 BC.

 Returns
 -------
 DateTime
     Tudat ``DateTime`` object.

 Examples
 --------
 In this example, the DateTime is constructed from a Julian day since January 1st 4713 BC.
 
 .. code-block:: python
 
     from tudatpy.astro.time_conversion import DateTime

     julian_day = 2451545.0 

     dt = DateTime.from_julian_day(julian_day)
     print(dt) # prints 2000-01-01 12:00:00.000000000000000
                         
                         )doc" )
            .def_static( "from_modified_julian_day",
                         &tba::DateTime::fromModifiedJulianDay,
                         py::arg( "modified_julian_day" ),
                         R"doc(

 Creates a Tudat-native :class:`DateTime` object from a modified Julian day.

 Parameters
 ----------
 modified_julian_day : float
     Date in modified Julian days (number of days since November 17th 1858).

 Returns
 -------
 DateTime
     Tudat ``DateTime`` object.

 Examples
 --------
 In this example, the DateTime is constructed from a modified Julian day since November 17th 1858.
 
 .. code-block:: python
 
     from tudatpy.astro.time_conversion import DateTime

     modified_julian_day = 51544.5

     dt = DateTime.from_modified_julian_day(modified_julian_day)
     print(dt) # prints 2000-01-01 12:00:00.000000000000000
                         
                         )doc" )
            .def( "to_python_datetime", &tba::DateTime::timePoint, R"doc(
                
Method to convert retrieve a Python datetime.datetime object from the Tudat :class:`DateTime` object. This is the inverse of the :meth:`~tudatpy.astro.time_conversion.DateTime.from_python_datetime` method.

.. note::

    The conversion uses the C++ `std::chrono` library, which is limited the time range it can represent. If the range is exceeded, the conversion will fail and throw an exception.
    
    The exact range is platform-dependent. On Windows, dates between 1970-01-01 and 3000-12-31 are allowed, on MacOS dates after 1900-01-01 are allowed and on Linux dates between 1678-01-01 and 2261-12-31 are allowed.

Returns
-------
datetime.datetime
    Datetime object, using the Python datetime library

    )doc" )
            .def( "to_days_since_reference_julian_day",
                  &tba::DateTime::daySinceReferenceJulianDay< double >,
                  py::arg( "reference_julian_day" ) = tba::JULIAN_DAY_ON_J2000,
                  R"doc(

 Convert a DateTime to Julian days since a given reference Julian date.


 Parameters
 ----------
 reference_julian_day : float, default = constants.JULIAN_DAY_ON_J2000
     Reference epoch (in days) since when the Julian days have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0) corresponding to the 1st of January 2000.
 Returns
 -------
 float
     Date in Julian days since the given reference epoch.

 Examples
 --------
 In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days since J2000 (the 1st of January 2000).

 .. code-block:: python

   from tudatpy.astro.time_conversion import DateTime

   # Define the calendar date using datetime
   dt = DateTime(2022, 5, 21, 13, 52, 41)
   # Convert the calendar date to Julian days since J2000
   julian_date = dt.to_days_since_reference_julian_day()
   # Print the converted output
   print(julian_date)  # prints 8176.07825231459


     )doc" )
            .def( "add_seconds", &tba::DateTime::addSecondsToDateTime< TIME_TYPE >, py::arg( "seconds_to_add" ), R"doc(
            
 Function to create a new Tudat :class:`DateTime` object by adding a number of seconds to an existing Tudat :class:`DateTime` object.

 .. note::
    
    This method does not modify the original :class:`DateTime` object, but returns a new one with the added seconds.

 Parameters
 ----------
 seconds_to_add : float
     Number of seconds to add
 Returns
 -------
 DateTime
     Tudat-native Datetime object created by adding the given number of seconds to the original DateTime
 
 Examples
 --------
 In this example, 86400 seconds are added to a DateTime object to construct a new DateTime.
 
 .. code-block:: python
 
     from tudatpy.astro.time_conversion import DateTime

     dt = DateTime(2025, 1, 1, 0, 0, 0.0)
     dt_seconds_added = dt.add_seconds(86400.0)
     print(f"Original dt: {dt}")
     print(f"dt with seconds added: {dt_seconds_added}")
     # prints:
     # Original dt: 2025-01-01 00:00:00.000000000000000
     # dt with seconds added: 2025-01-02 00:00:00.000000000000000   

            )doc" )
            .def( "add_days", &tba::DateTime::addDaysToDateTime< TIME_TYPE >, py::arg( "days_to_add" ), R"doc(
            
 Function to create a new Tudat :class:`DateTime` object by adding a number of days (86400 seconds) to an existing Tudat :class:`DateTime` object

 .. note::
    
    This method does not modify the original :class:`DateTime` object, but returns a new one with the added days.

 Parameters
 ----------
 days_to_add : float
     Number of days to add
 Returns
 -------
 DateTime
     Tudat-native Datetime object created by adding the given number of days to the original DateTime
 Examples
 --------
 In this example, 1 day is added to a DateTime object to construct a new DateTime.
 
 .. code-block:: python
 
     from tudatpy.astro.time_conversion import DateTime

     dt = DateTime(2025, 1, 1, 0, 0, 0.0)
     dt_days_added = dt.add_days(1.0)
     print(f"Original dt: {dt}")
     print(f"dt with days added: {dt_days_added}")
     # prints:
     # Original dt: 2025-01-01 00:00:00.000000000000000
     # dt with days added: 2025-01-02 00:00:00.000000000000000   

            )doc" );

    m.def( "julian_day_to_seconds_since_epoch",
           &tba::convertJulianDayToSecondsSinceEpoch< double >,
           py::arg( "julian_day" ),
           py::arg( "days_since_julian_day_zero" ) = tba::JULIAN_DAY_ON_J2000,
           R"doc(

 Convert Julian days to seconds since a given epoch.


 Parameters
 ----------
 julian_day : float
     Date in Julian days since January 1st 4713 BC.
 days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
     Reference epoch (in days since January 1st 4713 BC) since when the number of seconds have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0), corresponding to the 1st of January 2000.
 Returns
 -------
 float
     Seconds since the Julian date and the given epoch.





 Examples
 --------
 In this example, the Julian date `2459721.0783` (in days since January 1st 4713 BC), is converted to seconds since J2000 (January 1st 2000).

 .. code-block:: python

   # Define the Julian date in days since January 1st 4713 BC
   julian_date = 2459721.0783
   # Convert the Julian date to the number of seconds since J2000
   seconds_since_J2000 = time_conversion.julian_day_to_seconds_since_epoch(julian_date)
   # Print the converted output
   print(seconds_since_J2000)  # prints 706413165.1200145



     )doc" );

    m.def( "seconds_since_epoch_to_julian_day",
           &tba::convertSecondsSinceEpochToJulianDay< double >,
           py::arg( "seconds_since_epoch" ),
           py::arg( "days_since_julian_day_zero" ) = tba::JULIAN_DAY_ON_J2000,
           R"doc(

 Convert seconds since a given reference epoch to a Julian day.

 Inverse function of :func:`julian_day_to_seconds_since_epoch`.

 Parameters
 ----------
 seconds_since_epoch : float
     Seconds since ``days_since_julian_day_zero`` which are to be converted to date in Julian days.
 days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
     Reference epoch (in days since January 1st 4713 BC) since when the number of seconds have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0), corresponding to the 1st of January 2000.
 Returns
 -------
 float
     Date in Julian days since January 1st 4713 BC, as computed from the input parameters





 Examples
 --------
 In this example, an amount of seconds since J2000 (January 1st 2000) is converted to the Julian date (in days since January 1st 4713 BC).

 .. code-block:: python
 
   # Define the amount of seconds since January 1st 2000
   seconds_since_J2000 = 706413165.1200145
   # Convert the amount of seconds since J2000 to the Julian date
   julian_date = time_conversion.seconds_since_epoch_to_julian_day(seconds_since_J2000)
   # Print the converted output
   print(julian_date)  # prints 2459721.0783


     )doc" );

    m.def( "seconds_since_epoch_to_julian_years_since_epoch",
           &tba::convertSecondsSinceEpochToJulianYearsSinceEpoch< double >,
           py::arg( "seconds_since_epoch" ),
           R"doc(

 Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch.

 Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch. This is equivalent to converting a time interval in seconds to Julian years

 Parameters
 ----------
 seconds_since_epoch : float
     Seconds elapsed since a given (unspecified) epoch.
 Returns
 -------
 float
     Julian years since the specified epoch.

     Since this is a float, not a integer, meaning that the fraction of the year is also included.






 Examples
 --------
 In this example, `706413165.12` seconds since a given epoch are converted to Julian years since the same epoch.

 .. code-block:: python

   # Define the number of seconds elapsed
   seconds_since_epoch = 706413165.12
   # Convert the number of seconds to Julian years
   julian_years = time_conversion.seconds_since_epoch_to_julian_years_since_epoch(seconds_since_epoch)
   # Print the converted output
   print(julian_years)  # prints 22.38488240930869


     )doc" );

    m.def( "seconds_since_epoch_to_julian_centuries_since_epoch",
           &tba::convertSecondsSinceEpochToJulianCenturiesSinceEpoch< double >,
           py::arg( "seconds_since_epoch" ),
           R"doc(

 Convert the number of seconds since a given (unspecified) epoch to Julian centuries since the same epoch.

 Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch. This is equivalent to converting a time interval in seconds to Julian centuries

 Parameters
 ----------
 seconds_since_epoch : float
     Seconds elapsed since a given (unspecified) epoch.
 Returns
 -------
 float
     Julian centuries since the specified epoch.

     Since this is a float, not a integer, meaning that the fraction of the century is also included.






 Examples
 --------
 In this example, `706413165.12` seconds since a given epoch are converted to Julian centuries since the same epoch.

 .. code-block:: python

   # Define the number of seconds elapsed
   seconds_since_epoch = 706413165.12
   # Convert the number of seconds to Julian centuries
   julian_centuries = time_conversion.seconds_since_epoch_to_julian_centuries_since_epoch(seconds_since_epoch)
   # Print the converted output
   print(julian_centuries)  # prints 0.2238488240930869


     )doc" );

    m.def( "julian_day_to_modified_julian_day",
           &tba::convertJulianDayToModifiedJulianDay< double >,
           py::arg( "julian_day" ),
           R"doc(

 Convert a Julian day to a Modified Julian day.


 Parameters
 ----------
 julian_day : float
     Date in Julian days (number of days since January 1st 4713 BC).
 Returns
 -------
 float
     Date in modified Julian days (number of days since November 17th 1858).





 Examples
 --------
 In this example, the Julian date `2451545.0` (J2000) is converted to a modified Julian date.

 .. code-block:: python

   # Convert from Julian Days to Modified Julian Days
   MJD = time_conversion.julian_day_to_modified_julian_day(constants.JULIAN_DAY_ON_J2000)
   # Print the converted output
   print(MJD)  # prints 51544.5


     )doc" );

    m.def( "modified_julian_day_to_julian_day",
           &tba::convertModifiedJulianDayToJulianDay< double >,
           py::arg( "modified_julian_day" ),
           R"doc(

 Convert a Modified Julian day to a Julian day.

 Inverse function of :func:`julian_day_to_modified_julian_day`.

 Parameters
 ----------
 modified_julian_day : float
     Date in modified Julian days (number of days since November 17th 1858).
 Returns
 -------
 float
     Date in Julian days (number of days since January 1st 4713 BC).





 Examples
 --------
 In this example, the Modified Julian date `51544.5` ( corresponding to J2000) is converted to a modified Julian date.

 .. code-block:: python

   # Define J2000 in Modified Julian Days
   J2000_MJD = 51544.5
   # Convert from Modified Julian Days to Julian Days
   J2000 = time_conversion.modified_julian_day_to_julian_day(J2000_MJD)
   # Print the converted output
   print(J2000)  # prints 2451545.0


     )doc" );

    m.def( "is_leap_year",
           &tba::isLeapYear,
           py::arg( "year" ),
           R"doc(

 Assess wether a year is a leap year or not.


 Parameters
 ----------
 year : int
     Calendar year.
 Returns
 -------
 bool
     A value of True means that the year is a leap year.





 Examples
 --------
 In this example, the first list should contains only `True`, and the second `False`, since the first list uses leap years and the second does not.

 .. code-block:: python

   # Check known leap years
   leap_years = [time_conversion.is_leap_year(year) for year in [2020, 2016, 2000, 2400]]
   # Print the converted output
   print(leap_years)  # prints [True, True, True, True]
   # Check known non-leap years
   non_leap_years = [time_conversion.is_leap_year(year) for year in [2021, 2022, 2100, 2001]]
   # Print the converted output
   print(non_leap_years)  # prints [False, False, False, False]


     )doc" );

    m.def( "get_days_in_month",
           &tba::getDaysInMonth,
           py::arg( "month" ),
           py::arg( "year" ),
           R"doc(

 Get the number of days in the month of a given year.


 Parameters
 ----------
 month : int
     Calendar month.
 year : int
     Calendar year.
 Returns
 -------
 int
     Number of days in the month of the given year.





 Examples
 --------
 In this example, the number of days in February for both 2021 and 2020 are computed.

 .. code-block:: python

   # Check the number of days in February 2021
   days_feb_2021 = time_conversion.get_days_in_month(2, 2021)
   # Print the converted output
   print(days_feb_2021)  # prints 28
   # Check the number of days in February 2022
   days_feb_2020 = time_conversion.get_days_in_month(2, 2020)
   # Print the converted output
   print(days_feb_2020)  # prints 29


     )doc" );

    m.def( "calculate_seconds_in_current_julian_day",
           &tba::calculateSecondsInCurrentJulianDay,
           py::arg( "julian_day" ),
           R"doc(

 Determine the number of seconds that have elapsed in the given Julian day.


 Parameters
 ----------
 julian_day : float
     Date in Julian days (number of days since January 1st 4713 BC).
 Returns
 -------
 float
     Number of seconds that have passed in the given Julian day.





 Examples
 --------
 In this example, the number of seconds that have elapsed at the Julian day `2451545.2` is computed.

 .. code-block:: python

   # Compute the number of seconds that have passed in the given Julian day
   seconds_passed = time_conversion.calculate_seconds_in_current_julian_day(constants.JULIAN_DAY_ON_J2000)
   # Print the converted output
   print(seconds_passed)  # prints 43200.0


     )doc" );

    // Time scales conversion (inputs and outputs are always time in
    // seconds since J2000)
    m.def( "TCB_to_TDB",
           &tba::convertTcbToTdb< double >,
           py::arg( "TCB_time" ),
           R"doc(

 Convert time from the TCB scale to the TDB scale.

 The TCB scale is the Barycentric Coordinate Time, and the TDB scale is the Barycentric Dynamical Time.

 Parameters
 ----------
 TCB_time : float
     Time in seconds since J2000, in the TCB time scale.
 Returns
 -------
 float
     Time in seconds since J2000, in the TDB time scale.





 Examples
 --------
 In this example, the calendar date of the 17th of February 2022, at 15:41 and 2 seconds is first converted to Julian seconds since J2000.
 Then, this date and time is converted from the TCB scale to the TDB scale.

 .. code-block:: python

   # Define the date and time
   date = datetime.datetime(2022, 2, 17, 15, 41, 2)
   # Convert it in Julian days since J2000
   date_J2000 = time_conversion.python_datetime_to_julian_day(date)
   # Convert it in Julian seconds since J2000
   date_J2000_sec = time_conversion.julian_day_to_seconds_since_epoch(date_J2000)
   # Check the date from the TCB scale to the TDB scale
   date_TDB_scale = time_conversion.TCB_to_TDB(date_J2000_sec)
   # Print the converted output
   print(date_TDB_scale)  # prints 698384439.9176273


     )doc" );

    m.def( "TDB_to_TCB",
           &tba::convertTdbToTcb< double >,
           py::arg( "TDB_time" ),
           R"doc(

 Convert time from the TBD scale to the TCB scale.

 The TDB scale is the Barycentric Dynamical Time, and the TCB scale is the Barycentric Coordinate Time.

 Inverse function of :func:`TCB_to_TDB`.


 Parameters
 ----------
 TDB_time : float
     Time in seconds since J2000, in the TDB time scale.
 Returns
 -------
 float
     Time in seconds since J2000, in the TCB time scale.






     )doc" );

    m.def( "TCG_to_TT",
           &tba::convertTcgToTt< double >,
           py::arg( "TCG_time" ),
           R"doc(

 Convert time from the TCG scale to the TT scale.

 The TCG scale is the Geocentric Coordinate Time, and the TT scale is the Terrestrial Time.

 Parameters
 ----------
 TCG_time : float
     Time in seconds since J2000, in the TCG time scale.
 Returns
 -------
 float
     Time in seconds since J2000, in the TT time scale.






     )doc" );

    m.def( "TT_to_TCG",
           &tba::convertTtToTcg< double >,
           py::arg( "TT_time" ),
           R"doc(

 Convert time from the TT scale to the TCG scale.

 The TT scale is the Terrestrial Time, and the TCG scale is the Geocentric Coordinate Time.

 Inverse function of :func:`TCG_to_TT`.


 Parameters
 ----------
 TT_time : float
     Time in seconds since J2000, in the TT time scale.
 Returns
 -------
 float
     Time in seconds since J2000, in the TCG time scale.






     )doc" );

    m.def( "TAI_to_TT",
           &tba::convertTAItoTT< double >,
           py::arg( "TAI_time" ),
           R"doc(

 Convert time from the TAI scale to the TT scale.

 The TAI scale is the International Atomic Time, and the TT scale is the Terrestrial Time.

 Parameters
 ----------
 TAI_time : float
     Time in seconds since J2000, in the TAI time scale.
 Returns
 -------
 float
     Time in seconds since J2000, in the TT time scale.






     )doc" );

    m.def( "TT_to_TAI",
           &tba::convertTTtoTAI< double >,
           py::arg( "TT_time" ),
           R"doc(

 Convert time from the TT scale to the TAI scale.

 The TT scale is the Terrestrial Time, and the TAI scale is the International Atomic Time.

 Inverse function of :func:`TAI_to_TT`.


 Parameters
 ----------
 TT_time : float
     Time in seconds since J2000, in the TT time scale.
 Returns
 -------
 float
     Time in seconds since J2000, in the TAI time scale.






     )doc" );

    m.def( "TT_to_TDB_approximate",
           &tba::approximateConvertTTtoTDB,
           py::arg( "TT_time" ),
           R"doc(

 Approximately convert time from the TT scale to the TDB scale, using the following first-order approximation:

 .. math::
     t_{\text{TDB}} = t_{\text{TT}} + 0.001657  * \sin( 628.3076 T_{\text{TT}} + 6.2401 );

 with :math:`t` the epoch in seconds since J2000, :math:`T` the epoch in centuries since J2000.

 The TT scale is the Terrestrial Time, and the TDB scale is the Barycentric Dynamical Time.

 Parameters
 ----------
 TT_time : float
     Time in seconds since J2000, in the TT time scale.
 Returns
 -------
 float
     Time in seconds since J2000, in the TDB time scale.






     )doc" );

    m.def( "TT_to_TDB",
           &tsi::convertTTtoTDB< TIME_TYPE >,
           py::arg( "TT_time" ),
           py::arg( "earth_fixed_position" ),
           R"doc(

 Convert time from the TT scale to the TDB scale.

 Convert time from the TT scale to the TDB scale, using the `iauDtdb <https://www2.mpia-hd.mpg.de/~mathar/progs/sofa_api/group__SR.html#gaeab39417bb16e66c570232102a055f2d>`_ function in Sofa, which has an accuracy of 3 nanoseconds or better in the time span 1950-2050
 To call the Sofa function, we assume UTC=UT1 for the UT1 input. For this function, UT1 is only used in this sofa function to compute the local solar time, which in turn is used to compute the time-modulation of the topocentric term
 (dependent on the ground station position). Since these terms have an amplitude at the microsecond level, the error induced by this assumption is negligible given the inherent quality of the Sofa model.

 Parameters
 ----------
 TT_time : float
     Time in seconds since J2000, in the TT time scale.
 earth_fixed_position : numpy.ndarray, default=numpy.array([0, 0, 0])
     Earth-fixed position (e.g. in ITRF) that is used for detailed conversion between TDB and TT (induces a signature at the microsecond level)
 Returns
 -------
 float
     Time in seconds since J2000, in the TDB time scale.






     )doc" );

    m.def( "TDB_to_TT",
           &tsi::convertTDBtoTT< TIME_TYPE >,
           py::arg( "TDB_time" ),
           py::arg( "earth_fixed_position" ),
           R"doc(

 Convert time from the TDB scale to the TT scale.

 Convert time from the TT scale to the TDB scale, using the `iauDtdb <https://www2.mpia-hd.mpg.de/~mathar/progs/sofa_api/group__SR.html#gaeab39417bb16e66c570232102a055f2d>`_ function in Sofa, which has an accuracy of 3 nanoseconds or better in the time span 1950-2050
 To call the Sofa function, we assume UTC=UT1 and TT=TDB to compute UT1 from the input TDB. For this function, UT1 is only used in this sofa function to compute the local solar time, which in turn is used to compute the time-modulation of the topocentric term
 (dependent on the ground station position). Since these terms have an amplitude at the microsecond level, the error induced by this assumption is negligible given the inherent quality of the Sofa model.

 Parameters
 ----------
 TDB_time : float
     Time in seconds since J2000, in the TDB time scale.
 earth_fixed_position : numpy.ndarray, default=numpy.array([0, 0, 0])
     Earth-fixed position (e.g. in ITRF) that is used for detailed conversion between TDB and TT (induces a signature at the microsecond level)
 Returns
 -------
 float
     Time in seconds since J2000, in the TDB time scale.






     )doc" );

    m.def( "default_time_scale_converter",
           &teo::createDefaultTimeConverterPy,
           R"doc(

 Function to create a time-scale converter object with default settings. In particular, it uses default settings for conversion between UT1 and UTC:

 * Corrections for semi-diurnal variations due to libration for a non-rigid Earth as per Table 5.1b of IERS Conventions 2010
 * Corrections diurnal and semidiurnal variations due to ocean tides as per Tables 8.2a and 8.2b of the IERS Conventions 2010
 * For epoch 01-01-1962 and later: linear interpolation (correcting for discontinuities during days with leap seconds) of daily corrections for UTC-UT1 from the eopc04_14_IAU2000.62-now.txt file in the tudat-resources directory
 * For epochs before 01-01-1962, where UTC-UT1 is not available from these files, we use the values if for :math:\Delta T = UT1-TT` from `here <https://webspace.science.uu.nl/~gent0113/deltat/deltat.htm>`_ (back to year 1620). In this period, we set the approxomation UTC=UT1

 See :class:`~TimeScaleConverter` for specific functionality and options for time-scale conversions.

 Returns
 -------
 TimeScaleConverter
     Object to convert between different terrestrial time scales.





     )doc" );

    m.def( "date_time_components_to_epoch",
           &tba::timeFromDecomposedDateTime< TIME_TYPE >,
           py::arg( "year" ),
           py::arg( "month" ),
           py::arg( "day" ),
           py::arg( "hour" ),
           py::arg( "minute" ),
           py::arg( "seconds" ),
           R"doc(

 Computes the epoch as seconds since J2000 from the entries of the current date and time.

 Computes the epoch as seconds since J2000. This function is added for convenience, and is equivalent to ``DateTime(...).to_epoch()``.

 Parameters
 ----------
 year : int
     Calendar year

 month : int
     Calendar month (value must be 1-12)

 day : int
     Calendar day in current month, value must be larger than 0, and smaller or equal to the number of days in the month

 hour : int
     Full hours into the current day (value must be 0-23)

 minute : int
     Full minutes into the current hour (value must be 0-59)

 seconds : float
     Number of seconds into the current minute. Note that this value is stored as ``long double`` in Tudat, which may be 64-bit or 80-bit (16 or 19 digits) depending on the compiler used.

 Returns
 -------
 float
     Time in seconds since J2000.






     )doc" );

    m.def( "iso_string_to_epoch",
           &tba::timeFromIsoString< TIME_TYPE >,
           py::arg( "iso_datetime" ),
           R"doc(

 Computes the epoch as seconds since J2000 from an ISO datetime string.

 Computes the epoch as seconds since J2000. This function is added for convenience and is equivalent to ``DateTime.from_iso_string(...).to_epoch()``.

 Parameters
 ----------
 iso_datetime : str
     Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)

 Returns
 -------
 float
     Time in seconds since J2000.






     )doc" );

    //    m.def("epoch_from_julian_day",
    //          &tudat::timeFromJulianDay< TIME_TYPE >,
    //          py::arg("julian_day"),
    //          get_docstring("epoch_from_julian_day").c_str() );

    //    m.def("epoch_from_modified_julian_day",
    //          &tudat::timeFromModifiedJulianDay< TIME_TYPE >,
    //          py::arg("modified_julian_day"),
    //          get_docstring("epoch_from_modified_julian_day").c_str()
    //          );

    /////////////// DEPRECATED

    m.def( "date_time_from_epoch",
           &tba::DateTime::fromTime< TIME_TYPE >,
           py::arg( "epoch" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.from_epoch` instead.


 Creates a Tudat-native :class:`DateTime` object from the seconds since J2000.


 Parameters
 ----------
 epoch : float
     Seconds since J2000

 Returns
 -------
 DateTime
     Tudat ``DateTime`` object.






     )doc" );

    m.def( "date_time_from_iso_string",
           &tba::DateTime::fromIsoString,
           py::arg( "iso_string" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.from_iso_string` instead.

 Creates a Tudat-native :class:`DateTime` object from an ISO datetime string.


 Parameters
 ----------
 iso_datetime : str
     Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)

 Returns
 -------
 DateTime
     Tudat ``DateTime`` object.






     )doc" );

    m.def( "calendar_date_to_julian_day_since_epoch",
           &convertCalendarDateToJulianDaySinceEpochPy< double >,
           py::arg( "calendar_date" ),
           py::arg( "days_since_julian_day_zero" ) = tba::JULIAN_DAY_ON_J2000 );

    m.def( "calendar_date_to_days_since_epoch",
           &convertCalendarDateToJulianDaySinceEpochPy< double >,
           py::arg( "calendar_date" ),
           py::arg( "days_since_julian_day_zero" ) = tba::JULIAN_DAY_ON_J2000,
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_python_datetime(...).to_days_since_reference_julian_day()` instead.

 Convert a calendar date to Julian days since a given epoch.


 Parameters
 ----------
 calendar_date : datetime.datetime
     Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified. Milliseconds are ignored.
 days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
     Reference epoch (in days) since when the Julian days have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0) corresponding to the 1st of January 2000.
 Returns
 -------
 float
     Date in Julian days since the given epoch.





 Examples
 --------
 In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days since J2000 (the 1st of January 2000).

 .. code-block:: python

   # Define the calendar date using datetime
   calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
   # Convert the calendar date to Julian days since J2000
   julian_date = time_conversion.calendar_date_to_days_since_epoch(calendar_date)
   # Print the converted output
   print(julian_date)  # prints 8176.07825231459


     )doc" );

    m.def( "python_datetime_to_days_since_epoch",
           &convertCalendarDateToJulianDaySinceEpochPy< double >,
           py::arg( "datetime" ),
           py::arg( "days_since_julian_day_zero" ) = tba::JULIAN_DAY_ON_J2000,
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_python_datetime(...).to_days_since_reference_julian_day()` instead.


 Convert a calendar date to Julian days since a given epoch.


 Parameters
 ----------
 datetime : datetime.datetime
     Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified. Milliseconds are ignored.
 days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
     Reference epoch (in days) since when the Julian days have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0) corresponding to the 1st of January 2000.
 Returns
 -------
 float
     Date in Julian days since the given epoch.





 Examples
 --------
 In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days since J2000 (the 1st of January 2000).

 .. code-block:: python

   # Define the calendar date using datetime
   calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
   # Convert the calendar date to Julian days since J2000
   julian_date = time_conversion.python_datetime_to_days_since_epoch(calendar_date)
   # Print the converted output
   print(julian_date)  # prints 8176.07825231459


     )doc" );

    m.def( "julian_day_to_calendar_date",
           &convertJulianDayToCalendarDatePy,
           py::arg( "julian_day" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_julian_day(...).to_python_datetime()` instead.

 Convert Julian days to a calendar date.

 Inverse function of :func:`calendar_date_to_julian_day`.

 Parameters
 ----------
 julian_day : float
     Date in Julian days since January 1st 4713 BC.
 Returns
 -------
 datetime.datetime
     Datetime object, using the Python datetime library, containing the date and time corresponding to the Julian date input.





 Examples
 --------
 In this example, the Julian date `2459721.0783` (in days since January 1st 4713 BC), is converted to a calendar date.

 .. code-block:: python

   # Define the Julian date in days since January 1st 4713 BC
   julian_date = 2459721.0783
   # Convert the Julian date to a calendar date
   calendar_date = time_conversion.julian_day_to_calendar_date(julian_date)
   # Print the converted output
   print(calendar_date)  # prints datetime.datetime(2022, 5, 21, 13, 52, 45)


     )doc" );

    m.def( "julian_day_to_python_datetime",
           &convertJulianDayToCalendarDatePy,
           py::arg( "julian_day" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_julian_day(...).to_python_datetime()` instead.


 Convert Julian days to a calendar date.

 Inverse function of :func:`python_datetime_to_julian_day`.

 Parameters
 ----------
 julian_day : float
     Date in Julian days since January 1st 4713 BC.
 Returns
 -------
 datetime.datetime
     Datetime object, using the Python datetime library, containing the date and time corresponding to the Julian date input.





 Examples
 --------
 In this example, the Julian date `2459721.0783` (in days since January 1st 4713 BC), is converted to a calendar date.

 .. code-block:: python

   # Define the Julian date in days since January 1st 4713 BC
   julian_date = 2459721.0783
   # Convert the Julian date to a calendar date
   calendar_date = time_conversion.julian_day_to_python_datetime(julian_date)
   # Print the converted output
   print(calendar_date)  # prints datetime.datetime(2022, 5, 21, 13, 52, 45)


     )doc" );

    m.def( "calendar_date_to_julian_day",
           &convertCalendarDateToJulianDayPy< double >,
           py::arg( "calendar_date" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_python_datetime(...).to_julian_day()` instead.

 Convert a calendar date to Julian days.


 Parameters
 ----------
 calendar_date : datetime.datetime
     Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
 Returns
 -------
 float
     Julian day number (days since noon January 1st 4713 BC.)





 Examples
 --------
 In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days.

 .. code-block:: python

   # Define the calendar date using datetime
   calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
   # Convert the calendar date to Julian days since January 1st 4713 BC
   julian_date = time_conversion.calendar_date_to_julian_day(calendar_date)
   # Print the converted output
   print(julian_date)  # prints 2459721.0782523146


     )doc" );

    m.def( "python_datetime_to_julian_day",
           &convertCalendarDateToJulianDayPy< double >,
           py::arg( "datetime" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_python_datetime(...).to_julian_day()` instead.

 Convert a calendar date to Julian days.


 Parameters
 ----------
 datetime : datetime.datetime
     Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
 Returns
 -------
 float
     Julian day number (days since noon January 1st 4713 BC.)





 Examples
 --------
 In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days.

 .. code-block:: python

   # Define the calendar date using datetime
   calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
   # Convert the calendar date to Julian days since January 1st 4713 BC
   julian_date = time_conversion.python_datetime_to_julian_day(calendar_date)
   # Print the converted output
   print(julian_date)  # prints 2459721.0782523146


     )doc" );
    m.def( "datetime_to_tudat",
           &tba::DateTime::fromTimePoint,
           py::arg( "datetime" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.from_python_datetime` instead.

 Function to convert a Python datetime.datetime object to a Tudat :class:`DateTime` object. The Tudat-native alternative has the advantage of providing sub-femtosecond resolution, as opposed to the microsecond resolution of the Python version


 Parameters
 ----------
 datetime : datetime.datetime
     Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
 Returns
 -------
 DateTime
     DateTime object defined in Tudat

    )doc" );

    m.def( "year_and_days_in_year_to_calendar_date",
           &tba::DateTime::fromYearAndDaysInYear,
           py::arg( "year" ),
           py::arg( "days_in_year" ),
           R"doc(
        
.. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.from_year_and_day_of_year` instead.

Create the Tudat :class:`DateTime` from the year and the number of days in the year.

Parameters
----------
year : int
    Calendar year.
days_in_year : int
    Number of days that have passed in the year.
Returns
-------
DateTime
    Corresponding calendar date as a :class:`DateTime` object. Note: the hours, minutes and seconds in the object are set to 0 when calling this function.

Examples
--------
In this example, the calendar date corresponding to when 122 days have passed in 2020 is computed.

.. code-block:: python

    # Compute the calendar date when 122 days have passed in 2020
    currentDate = time_conversion.year_and_days_in_year_to_calendar_date(2020, 122)
    # Print the converted output
    print(currentDate)  # prints (2020, 5, 2, 0, 0)

     )doc" );

    m.def( "datetime_to_python",
           &dateTimeToTimePoint,
           py::arg( "datetime" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_python_datetime` instead.


 Function to convert a Tudat :class:`DateTime` object to a Python datetime.datetime object. This is the inverse of the :func:`datetime_to_tudat` function

 Parameters
 ----------
 datetime : DateTime
     Tudat-native Datetime object. Both the date and the time (hour, minutes, and seconds), can be specified, up to sub-femtosecond resolution.
 Returns
 -------
 datetime.datetime
     Datetime object, using the Python datetime library
     )doc" );

    m.def( "add_seconds_to_datetime",
           &tba::addSecondsToDateTime< TIME_TYPE >,
           py::arg( "datetime" ),
           py::arg( "seconds_to_add" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.add_seconds()` instead.

 Function to create a new Tudat :class:`DateTime` object by adding a number of seconds to an existing Tudat :class:`DateTime` object


 Parameters
 ----------
 datetime : DateTime
     Tudat-native Datetime object to which a number of seconds are to be added
 seconds_to_add : float
     Number of seconds to add
 Returns
 -------
 DateTime
     Tudat-native Datetime object created by adding the given number of seconds to the original DateTime






     )doc" );

    m.def( "add_days_to_datetime",
           &tba::addDaysToDateTime< TIME_TYPE >,
           py::arg( "datetime" ),
           py::arg( "days_to_add" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.add_days()` instead.

 Function to create a new Tudat :class:`DateTime` object by adding a number of days (86400 seconds) to an existing Tudat :class:`DateTime` object


 Parameters
 ----------
 datetime : DateTime
     Tudat-native Datetime object to which a number of days are to be added
 days_to_add : float
     Number of days to add
 Returns
 -------
 DateTime
     Tudat-native Datetime object created by adding the given number of days to the original DateTime






     )doc" );

    m.def( "epoch_from_date_time_components",
           &tba::timeFromDecomposedDateTime< TIME_TYPE >,
           py::arg( "year" ),
           py::arg( "month" ),
           py::arg( "day" ),
           py::arg( "hour" ),
           py::arg( "minute" ),
           py::arg( "seconds" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`date_time_components_to_epoch` instead.

 Computes the epoch as seconds since J2000 from the entries of the current date and time.

 Computes the epoch as seconds since J2000. This function is added for convenience, and creates a :class:`DateTime` object, and subsequently calls its ``epoch`` function

 Parameters
 ----------
 year : int
     Calendar year

 month : int
     Calendar month (value must be 1-12)

 day : int
     Calendar day in current month, value must be larger than 0, and smaller or equal to the number of days in the month

 hour : int
     Full hours into the current day (value must be 0-23)

 minute : int
     Full minutes into the current hour (value must be 0-59)

 seconds : float
     Number of seconds into the current minute. Note that this value is stored as ``long double`` in Tudat, which may be 64-bit or 80-bit (16 or 19 digits) depending on the compiler used.

 Returns
 -------
 float
     Time in seconds since J2000.






     )doc" );

    m.def( "epoch_from_date_time_iso_string",
           &tba::timeFromIsoString< TIME_TYPE >,
           py::arg( "iso_datetime" ),
           R"doc(

 .. warning::

    This function is deprecated and will be removed in a future version of Tudat. Use :func:`iso_string_to_epoch` instead.

 Computes the epoch as seconds since J2000 from an ISO datetime string.

 Computes the epoch as seconds since J2000. This function is added for convenience, and creates a :class:`DateTime` object, and subsequently calls its ``epoch`` function

 Parameters
 ----------
 iso_datetime : str
     Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)

 Returns
 -------
 float
     Time in seconds since J2000.

     )doc" );
}
}  // namespace time_conversion
}  // namespace astro
}  // namespace tudatpy
