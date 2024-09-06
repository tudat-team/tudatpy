/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <chrono>
#include <ctime>

#include "tudat/astro/basic_astro.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudatpy/scalarTypes.h"

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace pc = tudat::physical_constants;
namespace teo = tudat::earth_orientation;

tba::DateTime timePointToDateTime(
    const std::chrono::system_clock::time_point datetime) {
    std::time_t tt = std::chrono::system_clock::to_time_t(datetime);
    std::tm local_tm = *localtime(&tt);

    using namespace std::chrono;
    microseconds timeInMicroSeconds =
        duration_cast<microseconds>(datetime.time_since_epoch());
    long long fractional_seconds = timeInMicroSeconds.count() % 1000000LL;

    return tba::DateTime(
        local_tm.tm_year + 1900, local_tm.tm_mon + 1, local_tm.tm_mday,
        local_tm.tm_hour, local_tm.tm_min,
        static_cast<long double>(local_tm.tm_sec) +
            static_cast<long double>(fractional_seconds) /
                tudat::mathematical_constants::getFloatingInteger<long double>(
                    1000000LL));
}

// Convert from Gregorian date to time_point (Python datetime). Only
// year/month/day, no time.
std::chrono::system_clock::time_point dateTimeToTimePoint(
    const tba::DateTime& dateTime) {
    std::tm tm = {static_cast<int>(dateTime.getSeconds()),
                  dateTime.getMinute(),
                  dateTime.getHour(),
                  dateTime.getDay(),
                  dateTime.getMonth() - 1,
                  dateTime.getYear() - 1900

    };
    tm.tm_isdst = -1;
    std::chrono::system_clock::time_point timePoint =
        std::chrono::system_clock::from_time_t(std::mktime(&tm));
    return timePoint +
           std::chrono::microseconds(static_cast<int>(std::round(
               (dateTime.getSeconds() - static_cast<long double>(tm.tm_sec)) *
               tudat::mathematical_constants::getFloatingInteger<long double>(
                   1E6))));
}

// Convert Julian day to calendar date. This code ensures that the value
// returned is a time_point (Python datetime).
std::chrono::system_clock::time_point convertJulianDayToCalendarDatePy(
    const double julianDay) {
    tba::DateTime dateTime = tba::getCalendarDateFromTime<double>(
        tudat::timeFromJulianDay<double>(julianDay));

    return dateTimeToTimePoint(dateTime);
}

// Convert calendar date to Julian day since a given epoch. This code allows for
// the calendar date to be a time_point (Python datetime).
template <typename TimeScalarType = double>
TimeScalarType convertCalendarDateToJulianDayPy(
    const std::chrono::system_clock::time_point calendarDate) {
    tba::DateTime dateTime = timePointToDateTime(calendarDate);
    return dateTime.julianDay<TimeScalarType>();
}

template <typename TimeScalarType = double>
TimeScalarType convertCalendarDateToJulianDaySinceEpochPy(
    const std::chrono::system_clock::time_point calendarDate,
    const TimeScalarType epochSinceJulianDayZero =
        tba::getJulianDayOnJ2000<TimeScalarType>()) {
    tba::DateTime dateTime = timePointToDateTime(calendarDate);
    return dateTime.julianDay<TimeScalarType>() - epochSinceJulianDayZero;
}


namespace tudat {

    namespace earth_orientation {
        std::shared_ptr<TerrestrialTimeScaleConverter>
        createDefaultTimeConverterPy() {
            return createDefaultTimeConverter();
        }

    }  // namespace earth_orientation

}  // namespace tudat

namespace tudatpy {

    namespace astro {
        namespace time_conversion {

            PYBIND11_MODULE(expose_time_conversion, m) {
                //    m.attr("default_time_converter") =
                //    tudat::earth_orientation::defaultTimeConverter;


                m.def("datetime_to_tudat", &timePointToDateTime,
                      py::arg("datetime"),
"");

                m.def("datetime_to_python", &dateTimeToTimePoint,
                      py::arg("datetime"),
"");

                m.def("add_seconds_to_datetime",
                      &tba::addSecondsToDateTime<TIME_TYPE>,
                      py::arg("datetime"), py::arg("seconds_to_add"),
"");

                m.def("add_days_to_datetime",
                      &tba::addDaysToDateTime<TIME_TYPE>, py::arg("datetime"),
                      py::arg("days_to_add"),
"");


                m.def("calendar_date_to_julian_day",
                      &convertCalendarDateToJulianDayPy<double>,
                      py::arg("calendar_date"),
R"doc(Convert a calendar date to Julian days.



	:param calendar_date:
		Datetime, using the default Python library. Both the date and the time (hour, minutes, and seconds), can be specified. Milliseconds are ignored.
	:return:
		Date in Julian days since January 1st 4713 BC.
)doc");


                m.def("julian_day_to_calendar_date",
                      &convertJulianDayToCalendarDatePy, py::arg("julian_day"),
R"doc(Convert Julian days to a calendar date.

	Inverse function of :func:`calendar_date_to_julian_day`.

	:param julian_day:
		Date in Julian days since January 1st 4713 BC.
	:return:
		Calendar date using the Datetime Python library, containing the date and time corresponding to the Julian date input.
)doc");

                m.def(
                    "calendar_date_to_days_since_epoch",
                    &convertCalendarDateToJulianDaySinceEpochPy<double>,
                    py::arg("calendar_date"),
                    py::arg("days_since_julian_day_zero") =
                        tba::JULIAN_DAY_ON_J2000,
"");

                m.def(
                    "julian_day_to_seconds_since_epoch",
                    &tba::convertJulianDayToSecondsSinceEpoch<double>,
                    py::arg("julian_day"),
                    py::arg("days_since_julian_day_zero") =
                        tba::JULIAN_DAY_ON_J2000,
R"doc(Convert Julian days to seconds since a given epoch.



	:param julian_day:
		Date in Julian days since January 1st 4713 BC.
	:param epoch_since_julian_day_zero:
		Epoch since when the Julian days have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0), corresponding to the 1st of January 2000.
	:return:
		Seconds since the Julian date and the given epoch.
)doc");

                m.def(
                    "seconds_since_epoch_to_julian_day",
                    &tba::convertSecondsSinceEpochToJulianDay<double>,
                    py::arg("seconds_since_epoch"),
                    py::arg("days_since_julian_day_zero") =
                        tba::JULIAN_DAY_ON_J2000,
"");

                m.def("seconds_since_epoch_to_julian_years_since_epoch",
                      &tba::convertSecondsSinceEpochToJulianYearsSinceEpoch<
                          double>,
                      py::arg("seconds_since_epoch"),
R"doc(Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch.



	:param seconds_since_epoch:
		Seconds elapsed since a given (unspecified) epoch.
	:return:
		Julian years since the specified epoch.

		Since this is a float, not a integer, meaning that the fraction of the year is also included.

)doc");

                m.def("seconds_since_epoch_to_julian_centuries_since_epoch",
                      &tba::convertSecondsSinceEpochToJulianCenturiesSinceEpoch<
                          double>,
                      py::arg("seconds_since_epoch"),
R"doc(Convert the number of seconds since a given (unspecified) epoch to Julian centuries since the same epoch.



	:param seconds_since_epoch:
		Seconds elapsed since a given (unspecified) epoch.
	:return:
		Julian centuries since the specified epoch.

		Since this is a float, not a integer, meaning that the fraction of the century is also included.

)doc");

                m.def(
                    "julian_day_to_modified_julian_day",
                    &tba::convertJulianDayToModifiedJulianDay<double>,
                    py::arg("julian_day"),
R"doc(Convert a Julian day to a Modified Julian day.



	:param julian_day:
		Date in Julian days (number of days since January 1st 4713 BC).
	:return:
		Date in modified Julian days (number of days since November 17th 1858).
)doc");

                m.def(
                    "modified_julian_day_to_julian_day",
                    &tba::convertModifiedJulianDayToJulianDay<double>,
                    py::arg("modified_julian_day"),
R"doc(Convert a Modified Julian day to a Julian day.

	Inverse function of :func:`julian_day_to_modified_julian_day`.

	:param modified_julian_day:
		Date in modified Julian days (number of days since November 17th 1858).
	:return:
		Date in Julian days (number of days since January 1st 4713 BC).
)doc");

                m.def("is_leap_year", &tba::isLeapYear, py::arg("year"),
R"doc(Assess wether a year is a leap year or not.



	:param year:
		Calendar year.
	:return:
		A value of True means that the year is a leap year.
)doc");

                m.def("get_days_in_month", &tba::getDaysInMonth,
                      py::arg("month"), py::arg("year"),
R"doc(Get the number of days in the month of a given year.



	:param month:
		Calendar month.
	:param year:
		Calendar year.
	:return:
		Number of days in the month of the given year.
)doc");

                m.def("calculate_seconds_in_current_julian_day",
                      &tba::calculateSecondsInCurrentJulianDay,
                      py::arg("julian_day"),
R"doc(Determine the number of seconds that have elapsed in the given Julian day.



	:param julian_day:
		Date in Julian days (number of days since January 1st 4713 BC).
	:return:
		Number of seconds that have passed in the given Julian day.
)doc");

                // Time scales conversion (inputs and outputs are always time in
                // seconds since J2000)
                m.def("TCB_to_TDB", &tba::convertTcbToTdb<double>,
R"doc(Convert time from the TCB scale to the TDB scale.

	The TCB scale is the Barycentric Coordinate Time, and the TDB scale is the Barycentric Dynamical Time.

	:param TCB_time:
		Time in seconds since J2000, in the TCB time scale.
	:return:
		Time in seconds since J2000, in the TDB time scale.
)doc");

                m.def("TDB_to_TCB", &tba::convertTdbToTcb<double>,
R"doc(Convert time from the TBD scale to the TCB scale.

	The TDB scale is the Barycentric Dynamical Time, and the TCB scale is the Barycentric Coordinate Time.

	Inverse function of :func:`TCB_to_TDB`.


	:param TDB_time:
		Time in seconds since J2000, in the TDB time scale.
	:return:
		Time in seconds since J2000, in the TCB time scale.
)doc");

                m.def("TCG_to_TT", &tba::convertTcgToTt<double>,
R"doc(Convert time from the TCG scale to the TT scale.

	The TCG scale is the Geocentric Coordinate Time, and the TT scale is the Terrestrial Time.

	:param TCG_time:
		Time in seconds since J2000, in the TCG time scale.
	:return:
		Time in seconds since J2000, in the TT time scale.
)doc");

                m.def("TT_to_TCG", &tba::convertTtToTcg<double>,
R"doc(Convert time from the TT scale to the TCG scale.

	The TT scale is the Terrestrial Time, and the TCG scale is the Geocentric Coordinate Time.

	Inverse function of :func:`TCG_to_TT`.


	:param TT_time:
		Time in seconds since J2000, in the TT time scale.
	:return:
		Time in seconds since J2000, in the TCG time scale.
)doc");

                m.def("TAI_to_TT", &tba::convertTAItoTT<double>,
R"doc(Convert time from the TAI scale to the TT scale.

	The TAI scale is the International Atomic Time, and the TT scale is the Terrestrial Time.

	:param TAI_time:
		Time in seconds since J2000, in the TAI time scale.
	:return:
		Time in seconds since J2000, in the TT time scale.
)doc");

                m.def("TT_to_TAI", &tba::convertTTtoTAI<double>,
R"doc(Convert time from the TT scale to the TAI scale.

	The TT scale is the Terrestrial Time, and the TAI scale is the International Atomic Time.

	Inverse function of :func:`TAI_to_TT`.


	:param TT_time:
		Time in seconds since J2000, in the TT time scale.
	:return:
		Time in seconds since J2000, in the TAI time scale.
)doc");

                m.def("TT_to_TDB_approximate", &tba::approximateConvertTTtoTDB,
                      py::arg("TT_time"),
R"doc(Approximately convert time from the TT scale to the TDB scale.

	The TT scale is the Terrestrial Time, and the TDB scale is the Barycentric Dynamical Time.

	:param TT_time:
		Time in seconds since J2000, in the TT time scale.
	:return:
		Time in seconds since J2000, in the TDB time scale.
)doc");


                py::enum_<tba::TimeScales>(m, "TimeScales")
                    .value("tai_scale", tba::tai_scale)
                    .value("tt_scale", tba::tt_scale)
                    .value("tdb_scale", tba::tdb_scale)
                    .value("utc_scale", tba::utc_scale)
                    .value("ut1_scale", tba::ut1_scale)
                    .export_values();

                m.def("default_time_scale_converter",
                      &teo::createDefaultTimeConverterPy,
"");

                py::class_<teo::TerrestrialTimeScaleConverter,
                           std::shared_ptr<teo::TerrestrialTimeScaleConverter>>(
                    m, "TimeScaleConverter",
"")
                    .def("convert_time",
                         &teo::TerrestrialTimeScaleConverter::getCurrentTime<
                             double>,
                         py::arg("input_scale"), py::arg("output_scale"),
                         py::arg("input_value"),
                         py::arg("earth_fixed_position") =
                             Eigen::Vector3d::Zero())
                    .def("get_time_difference",
                         &teo::TerrestrialTimeScaleConverter::
                             getCurrentTimeDifference<double>,
                         py::arg("input_scale"), py::arg("output_scale"),
                         py::arg("input_value"),
                         py::arg("earth_fixed_position") =
                             Eigen::Vector3d::Zero());


                py::class_<tba::DateTime>(m, "DateTime",
"")
                    .def(py::init<const int, const int, const int, const int,
                                  const int, const long double>(),
                         py::arg("year"), py::arg("month"), py::arg("day"),
                         py::arg("hour") = 12, py::arg("minute") = 0,
                         py::arg("seconds") = 0.0L)
                    .def_property("year", &tba::DateTime::getYear,
                                  &tba::DateTime::setYear,
"")
                    .def_property("month", &tba::DateTime::getMonth,
                                  &tba::DateTime::setMonth,
"")
                    .def_property("day", &tba::DateTime::getDay,
                                  &tba::DateTime::setDay,
"")
                    .def_property("hour", &tba::DateTime::getHour,
                                  &tba::DateTime::setHour,
"")
                    .def_property("minute", &tba::DateTime::getMinute,
                                  &tba::DateTime::setMinute,
"")
                    .def_property("seconds", &tba::DateTime::getSeconds,
                                  &tba::DateTime::setSeconds,
"")
                    .def("iso_string", &tba::DateTime::isoString,
                         py::arg("add_T") = false,
                         py::arg("number_of_digits_seconds") = 15,
"")
                    .def("day_of_year", &tba::DateTime::dayOfYear,
"")
                    .def("epoch", &tba::DateTime::epoch<TIME_TYPE>,
"")
                    .def("julian_day", &tba::DateTime::julianDay<double>,
"")
                    .def("modified_julian_day",
                         &tba::DateTime::modifiedJulianDay<double>,
"");


                m.def("epoch_from_date_time_components",
                      &tba::timeFromDecomposedDateTime<TIME_TYPE>,
                      py::arg("year"), py::arg("month"), py::arg("day"),
                      py::arg("hour"), py::arg("minute"), py::arg("seconds"),
"");

                m.def("epoch_from_date_time_iso_string",
                      &tba::timeFromIsoString<TIME_TYPE>,
                      py::arg("iso_datetime"),
"");


                //    m.def("epoch_from_julian_day",
                //          &tudat::timeFromJulianDay< TIME_TYPE >,
                //          py::arg("julian_day"),

                //    m.def("epoch_from_modified_julian_day",
                //          &tudat::timeFromModifiedJulianDay< TIME_TYPE >,
                //          py::arg("modified_julian_day"),
                //          );

                m.def("date_time_from_epoch",
                      &tba::getCalendarDateFromTime<TIME_TYPE>,
                      py::arg("epoch"),
"");

                m.def("date_time_from_iso_string",
                      &tba::getCalendarDateFromTime<TIME_TYPE>,
                      py::arg("iso_datetime"),
"");


                /////////////// DEPRECATED

                m.def("calendar_date_to_julian_day_since_epoch",
                      &convertCalendarDateToJulianDaySinceEpochPy<double>,
                      py::arg("calendar_date"),
                      py::arg("days_since_julian_day_zero") =
                          tba::JULIAN_DAY_ON_J2000);
            }
        }  // namespace time_conversion
    }  // namespace astro
}  // namespace tudatpy