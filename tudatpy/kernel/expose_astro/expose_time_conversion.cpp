/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <chrono>
#include <ctime>

#include "tudatpy/docstrings.h"

#include "expose_time_conversion.h"
#include "tudatpy/scalarTypes.h"

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/basic_astro/timeConversions.h>
#include <tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h>
#include <tudat/math/basic/mathematicalConstants.h>

#include <pybind11/chrono.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace pc  = tudat::physical_constants;
namespace teo = tudat::earth_orientation;

tba::DateTime timePointToDateTime(const std::chrono::system_clock::time_point datetime)
{
    std::time_t tt = std::chrono::system_clock::to_time_t(datetime);
    std::tm local_tm = *localtime(&tt);

    using namespace std::chrono;
    microseconds timeInMicroSeconds = duration_cast<microseconds>(datetime.time_since_epoch());
    long long fractional_seconds = timeInMicroSeconds.count() % 1000000LL;

    return tba::DateTime( local_tm.tm_year + 1900, local_tm.tm_mon + 1, local_tm.tm_mday,
        local_tm.tm_hour, local_tm.tm_min, static_cast< long double >( local_tm.tm_sec ) +
        static_cast< long double >( fractional_seconds ) /
        tudat::mathematical_constants::getFloatingInteger< long double >( 1000000 ) );
}

// Convert from Gregorian date to time_point (Python datetime). Only year/month/day, no time.
std::chrono::system_clock::time_point dateTimeToTimePoint(const tba::DateTime& dateTime )
{
    std::tm tm = {
        static_cast< int >( dateTime.getSeconds( ) ),
        dateTime.getMinute( ),
        dateTime.getHour( ),
        dateTime.getDay( ),
        dateTime.getMonth( ) - 1,
        dateTime.getYear( ) - 1900

    };
    tm.tm_isdst = -1;
    std::chrono::system_clock::time_point timePoint = std::chrono::system_clock::from_time_t(std::mktime(&tm));
    return timePoint + std::chrono::microseconds ( static_cast< int >( std::round(
                                                     ( dateTime.getSeconds( ) - static_cast< long double >( tm.tm_sec ) ) *
                                                     tudat::mathematical_constants::getFloatingInteger< long double >( 1E6 ) ) ) );
}

// Convert Julian day to calendar date. This code ensures that the value returned is a time_point (Python datetime).
std::chrono::system_clock::time_point convertJulianDayToCalendarDatePy( const double julianDay )
{
    tba::DateTime dateTime =
        tba::getCalendarDateFromTime< double >( tudat::timeFromJulianDay< double >( julianDay ) );

    return dateTimeToTimePoint(dateTime);
}

// Convert calendar date to Julian day since a given epoch. This code allows for the calendar date to be a time_point (Python datetime).
template< typename TimeScalarType = double >
TimeScalarType convertCalendarDateToJulianDayPy(
    const std::chrono::system_clock::time_point calendarDate )
{
    tba::DateTime dateTime = timePointToDateTime( calendarDate );
    return dateTime.julianDay< TimeScalarType >( );
}

namespace tudat
{

namespace earth_orientation
{
std::shared_ptr<TerrestrialTimeScaleConverter> createDefaultTimeConverterPy( )
{
    return createDefaultTimeConverter( );
}

}

}

namespace tudatpy {

namespace astro {
namespace time_conversion {

void expose_time_conversion(py::module &m) {

//    m.attr("default_time_converter") = tudat::earth_orientation::defaultTimeConverter;


    m.def("datetime_to_tudat",
          &timePointToDateTime,
          py::arg("datetime"),
          get_docstring("datetime_to_tudat").c_str()
    );

    m.def("datetime_to_python",
          &dateTimeToTimePoint,
          py::arg("datetime"),
          get_docstring("datetime_to_python").c_str()
    );

    m.def("add_seconds_to_datetime",
          &tba::addSecondsToDateTime< TIME_TYPE >,
          py::arg("datetime"),
          py::arg("seconds_to_add"),
          get_docstring("add_seconds_to_datetime").c_str()
    );

//    m.def("add_days_to_datetime",
//          &tba::addDaysToDateTime< TIME_TYPE >,
//          py::arg("datetime"),
//          py::arg("days_to_add"),
//          get_docstring("add_days_to_datetime").c_str()
//    );






    m.def("calendar_date_to_julian_day",
          &convertCalendarDateToJulianDayPy< double >,
          py::arg("calendar_date"),
          get_docstring("calendar_date_to_julian_day").c_str()
          );

    m.def("calendar_date_to_julian_day_since_epoch",
          &convertCalendarDateToJulianDayPy< double >,
          py::arg("calendar_date"),
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

    // Time scales conversion (inputs and outputs are always time in seconds since J2000)
    m.def("TCB_to_TDB",
          &tba::convertTcbToTdb< double >,
          py::arg("TCB_time"),
          get_docstring("TCB_to_TDB").c_str()
          );

    m.def("TDB_to_TCB",
          &tba::convertTdbToTcb< double >,
          py::arg("TDB_time"),
          get_docstring("TDB_to_TCB").c_str()
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


    py::enum_<tba::TimeScales>(
                m, "TimeScales" )
            .value("tai_scale", tba::tai_scale)
            .value("tt_scale", tba::tt_scale)
            .value("tdb_scale", tba::tdb_scale)
            .value("utc_scale", tba::utc_scale)
            .value("ut1_scale", tba::ut1_scale)
            .export_values();

    m.def("default_time_scale_converter",
          &teo::createDefaultTimeConverterPy,
          get_docstring("default_time_scale_converter").c_str() );

    py::class_<teo::TerrestrialTimeScaleConverter,
            std::shared_ptr<teo::TerrestrialTimeScaleConverter>>(
                m, "TimeScaleConverter",
                get_docstring("TimeScaleConverter").c_str( ) )
            .def( "convert_time", &teo::TerrestrialTimeScaleConverter::getCurrentTime< double >,
                  py::arg( "input_scale" ),
                  py::arg( "output_scale" ),
                  py::arg( "input_value" ),
                  py::arg( "earth_fixed_position" ) = Eigen::Vector3d::Zero( ) );

    py::class_< tba::DateTime >( m,"DateTime", get_docstring("DateTime").c_str())
        .def(py::init<
                 const int,
                 const int,
                 const int,
                 const int,
                 const int,
                 const long double>(),
             py::arg("year"), py::arg("month"), py::arg("day"),
             py::arg("hour") = 12, py::arg("minute") = 0, py::arg("seconds") = 0.0L )
        .def_property("year", &tba::DateTime::getYear, &tba::DateTime::setYear, get_docstring("DateTime.year").c_str() )
        .def_property("month", &tba::DateTime::getMonth, &tba::DateTime::setMonth, get_docstring("DateTime.month").c_str() )
        .def_property("day", &tba::DateTime::getDay, &tba::DateTime::setDay, get_docstring("DateTime.day").c_str() )
        .def_property("hour", &tba::DateTime::getHour, &tba::DateTime::setHour, get_docstring("DateTime.hour").c_str() )
        .def_property("minute", &tba::DateTime::getMinute, &tba::DateTime::setMinute, get_docstring("DateTime.minute").c_str() )
        .def_property("seconds", &tba::DateTime::getSeconds, &tba::DateTime::setSeconds, get_docstring("DateTime.seconds").c_str() )
        .def("iso_string",
             &tba::DateTime::isoString,
             py::arg("add_T") = false,
             py::arg("number_of_digits_seconds") = 15,
             get_docstring("DateTime.iso_string").c_str( ) )
        .def("day_of_year",
             &tba::DateTime::dayOfYear,
             get_docstring("DateTime.day_of_year").c_str( ) )
        .def("epoch",
             &tba::DateTime::epoch< TIME_TYPE >,
                 get_docstring("DateTime.epoch").c_str( ) )
        .def("julian_day",
             &tba::DateTime::julianDay< double >,
             get_docstring("DateTime.julian_day").c_str( ) )
        .def("modified_julian_day",
             &tba::DateTime::modifiedJulianDay< double >,
             get_docstring("DateTime.modified_julian_day").c_str( ) );


    m.def("epoch_from_date_time_components",
          &tba::timeFromDecomposedDateTime< TIME_TYPE >,
          py::arg("year"), py::arg("month"), py::arg("day"),
          py::arg("hour"), py::arg("minute"), py::arg("seconds"),
          get_docstring("epoch_from_date_time_components").c_str() );

    m.def("epoch_from_date_time_iso_string",
          &tba::timeFromIsoString< TIME_TYPE >,
          py::arg("iso_datetime"),
          get_docstring("epoch_from_date_time_iso_string").c_str() );


//    m.def("epoch_from_julian_day",
//          &tudat::timeFromJulianDay< TIME_TYPE >,
//          py::arg("julian_day"),
//          get_docstring("epoch_from_julian_day").c_str() );

//    m.def("epoch_from_modified_julian_day",
//          &tudat::timeFromModifiedJulianDay< TIME_TYPE >,
//          py::arg("modified_julian_day"),
//          get_docstring("epoch_from_modified_julian_day").c_str() );

    m.def("date_time_from_epoch",
          &tba::getCalendarDateFromTime< TIME_TYPE >,
          py::arg("epoch"),
          get_docstring("date_time_from_epoch").c_str() );

    m.def("date_time_from_iso_string",
          &tba::getCalendarDateFromTime< TIME_TYPE >,
          py::arg("iso_datetime"),
          get_docstring("date_time_from_iso_string").c_str() );



}
} // namespace time_conversion
} // namespace astro
} // namespace tudatpy
