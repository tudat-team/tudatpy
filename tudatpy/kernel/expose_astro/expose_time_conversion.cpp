/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudatpy/docstrings.h"

#include "expose_time_conversion.h"

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/basic_astro/timeConversions.h>

#include <pybind11/chrono.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tba = tudat::basic_astrodynamics;
namespace pc  = tudat::physical_constants;

// Convert a time_point (which is automatically converted from a Python datetime) to gregorian::date.
// Mind that the Gregorian date only contain the year, month, and day. Hour, minutes, and seconds, are lost in the conversion.
boost::gregorian::date timePointToGregorianDate(const std::chrono::system_clock::time_point datetime) {
  std::time_t tt = std::chrono::system_clock::to_time_t(datetime);
  std::tm local_tm = *localtime(&tt);  
  return boost::gregorian::date(local_tm.tm_year + 1900, local_tm.tm_mon + 1, local_tm.tm_mday);
}

// Convert from Gregorian date to time_point (Python datetime). Only year/month/day, no time.
std::chrono::system_clock::time_point GregorianDateToTimePoint(const boost::gregorian::date gregorianDate) {
  std::tm local_tm = boost::gregorian::to_tm(gregorianDate);
  std::time_t tt = std::mktime(&local_tm);
  return std::chrono::system_clock::from_time_t(tt);
}

// Convert Julian day to calendar date. This code ensures that the value returned is a time_point (Python datetime).
std::chrono::system_clock::time_point convertJulianDayToCalendarDatePy(const double julianDay) {
  boost::gregorian::date gregorianDate = tba::convertJulianDayToCalendarDate(julianDay);
  return GregorianDateToTimePoint(gregorianDate);
}

// Convert calendar date to Julian day. This code allows for the input to be a time_point (Python datetime).
template< typename TimeScalarType = double >
TimeScalarType convertCalendarDateToJulianDayPy(const std::chrono::system_clock::time_point calendarDate) {
  std::time_t tt = std::chrono::system_clock::to_time_t(calendarDate);
  std::tm local_tm = *localtime(&tt);
  return tba::convertCalendarDateToJulianDay(local_tm.tm_year + 1900, local_tm.tm_mon + 1, local_tm.tm_mday,
    local_tm.tm_hour, local_tm.tm_min, local_tm.tm_sec);
}

// Convert calendar date to Julian day since a given epoch. This code allows for the calendar date to be a time_point (Python datetime).
template< typename TimeScalarType = double >
TimeScalarType convertCalendarDateToJulianDaySinceEpochPy( const std::chrono::system_clock::time_point calendarDate,
                                                           const TimeScalarType epochSinceJulianDayZero = tba::JULIAN_DAY_ON_J2000) {
  std::time_t tt = std::chrono::system_clock::to_time_t(calendarDate);
  std::tm local_tm = *localtime(&tt);
  const boost::gregorian::date gregorianDate = boost::gregorian::date(local_tm.tm_year + 1900, local_tm.tm_mon + 1, local_tm.tm_mday);
  const TimeScalarType fractionOfDay = (local_tm.tm_hour*3600+local_tm.tm_min*60+local_tm.tm_sec)/pc::JULIAN_DAY;
  return tba::calculateJulianDaySinceEpoch(gregorianDate, fractionOfDay, epochSinceJulianDayZero);
}

// Compute which day of the year a given calendar date is. This block allows for the calendar date to be a time_point (Python datetime).
double convertDayMonthYearToDayOfYearPy(std::chrono::system_clock::time_point calendarDate) {
  boost::gregorian::date gregorianDate = timePointToGregorianDate(calendarDate);
  return tba::convertDayMonthYearToDayOfYear(gregorianDate);
}

// Convert a year and number of days in a year to a calendar date. This code convert the output to a time_point (Python datetime).
std::chrono::system_clock::time_point convertYearAndDaysInYearToDatePy(const int year, const int daysInYear) {
  boost::gregorian::date gregorianDate = tba::convertYearAndDaysInYearToDate(year, daysInYear);
  return GregorianDateToTimePoint(gregorianDate);
}

namespace tudatpy {

namespace astro {
namespace time_conversion {

void expose_time_conversion(py::module &m) {

    m.def("calendar_date_to_julian_day",
          &convertCalendarDateToJulianDayPy< double >,
          py::arg("calendar_date"),
          get_docstring("calendar_date_to_julian_day").c_str()
      );

    m.def("calendar_date_to_julian_day_since_epoch",
          &convertCalendarDateToJulianDaySinceEpochPy< double >,
          py::arg("calendar_date"),
          py::arg("epoch_since_julian_day_zero") = tba::JULIAN_DAY_ON_J2000,
          get_docstring("calendar_date_to_julian_day_since_epoch").c_str()
      );

    m.def("julian_day_to_calendar_date",
          &convertJulianDayToCalendarDatePy,
          py::arg("julian_day"),
          get_docstring("julian_day_to_calendar_date").c_str()
      );

    m.def("julian_day_to_seconds_since_epoch",
          &tba::convertJulianDayToSecondsSinceEpoch< double >,
          py::arg("julian_day"),
          py::arg("epoch_since_julian_day_zero") = tba::JULIAN_DAY_ON_J2000,
          get_docstring("julian_day_to_seconds_since_epoch").c_str()
      );

    m.def("seconds_since_epoch_to_julian_day",
          &tba::convertSecondsSinceEpochToJulianDay< double >,
          py::arg("seconds_since_epoch"),
          py::arg("epoch_since_julian_day_zero") = tba::JULIAN_DAY_ON_J2000,
          get_docstring("seconds_since_epoch_to_julian_day").c_str()
      );

    m.def("seconds_since_epoch_to_julian_years_since_epoch",
          &tba::convertSecondsSinceEpochToJulianYearsSinceEpoch< double >,
          py::arg("seconds_since_epoch"),
          get_docstring("seconds_since_epoch_to_julian_years_since_epoch").c_str()
      );

    m.def("seconds_since_epoch_to_julian_centuries_since_epoch",
          &tba::convertSecondsSinceEpochToJulianCenturiesSinceEpoch< double >,
          py::arg("seconds_since_epoch"),
          get_docstring("seconds_since_epoch_to_julian_centuries_since_epoch").c_str()
      );

    // m.def("calendar_date_to_julian_day_since_epoch",
    //       &tba::convertCalendarDateToJulianDaysSinceEpoch< double >,
    //       py::arg("calendar_year"),
    //       py::arg("calendar_month"),
    //       py::arg("calendar_day"),
    //       py::arg("calendar_hour"),
    //       py::arg("calendar_minutes"),
    //       py::arg("calendar_seconds"),
    //       py::arg("reference_julian_day")
    //   );

    // m.def("calendar_date_to_julian_day",
    //       &tba::convertCalendarDateToJulianDay< double >,
    //       py::arg("calendar_year"),
    //       py::arg("calendar_month"),
    //       py::arg("calendar_day"),
    //       py::arg("calendar_hour"),
    //       py::arg("calendar_minutes"),
    //       py::arg("calendar_seconds")
    //   );

    m.def("julian_day_to_modified_julian_day",
          &tba::convertJulianDayToModifiedJulianDay< double >,
          py::arg("julian_day"),
          get_docstring("julian_day_to_modified_julian_day").c_str()
      );

    m.def("modified_julian_day_to_julian_day",
          &tba::convertModifiedJulianDayToJulianDay< double >,
          py::arg("modified_julian_day"),
          get_docstring("modified_julian_day_to_julian_day").c_str()
      );

    m.def("is_leap_year",
          &tba::isLeapYear,
          py::arg("year"),
          get_docstring("is_leap_year").c_str()
      );

    m.def("get_days_in_month",
          &tba::getDaysInMonth,
          py::arg("month"),
          py::arg("year"),
          get_docstring("get_days_in_month").c_str()
      );

    m.def("calculate_seconds_in_current_julian_day",
          &tba::calculateSecondsInCurrentJulianDay,
          py::arg("julian_day"),
          get_docstring("calculate_seconds_in_current_julian_day").c_str()
      );

    m.def("calendar_date_to_day_of_year",
          &convertDayMonthYearToDayOfYearPy,
          py::arg("calendar_date"),
          get_docstring("calendar_date_to_day_of_year").c_str()
      );

    m.def("year_and_days_in_year_to_calendar_date",
          &convertYearAndDaysInYearToDatePy,
          py::arg("year"),
          py::arg("days_in_year"),
          get_docstring("year_and_days_in_year_to_calendar_date").c_str()
      );

    // Time scales conversion (inputs and outputs are always time in seconds since J2000)
    m.def("TCB_to_TDB",
          &tba::convertTcbToTdb< double >,
          py::arg("TCB_time"),
          get_docstring("TCB_to_TDB").c_str()
      );

    m.def("TBD_to_TCB",
          &tba::convertTdbToTcb< double >,
          py::arg("TBD_time"),
          get_docstring("TBD_to_TCB").c_str()
      );

    m.def("TCG_to_TT",
          &tba::convertTcgToTt< double >,
          py::arg("TCG_time"),
          get_docstring("TCG_to_TT").c_str()
      );

    m.def("TT_to_TCG",
          &tba::convertTtToTcg< double >,
          py::arg("TT_time"),
          get_docstring("TT_to_TCG").c_str()
      );

    m.def("TAI_to_TT",
          &tba::convertTAItoTT< double >,
          py::arg("TAI_time"),
          get_docstring("TAI_to_TT").c_str()
      );

    m.def("TT_to_TAI",
          &tba::convertTTtoTAI< double >,
          py::arg("TT_time"),
          get_docstring("TT_to_TAI").c_str()
      );

    m.def("TT_to_TDB_approximate",
          &tba::approximateConvertTTtoTDB,
          py::arg("TT_time"),
          get_docstring("TT_to_TDB_approximate").c_str()
      );

}
} // namespace time_conversion
} // namespace astro
} // namespace tudatpy
