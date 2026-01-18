/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "../../wasm_module.h"
#include "../../shared_ptr_wasm.h"

#include <tudat/astro/basic_astro/timeConversions.h>
#include <tudat/astro/basic_astro/dateTime.h>

namespace tba = tudat::basic_astrodynamics;

WASM_MODULE_PATH("astro_time_representation")

EMSCRIPTEN_BINDINGS(tudatpy_astro_time_representation) {
    using namespace emscripten;

    // TimeScales enum
    enum_<tba::TimeScales>("astro_time_representation_TimeScales")
        .value("tai_scale", tba::tai_scale)
        .value("tt_scale", tba::tt_scale)
        .value("tdb_scale", tba::tdb_scale)
        .value("utc_scale", tba::utc_scale)
        .value("ut1_scale", tba::ut1_scale);

    // Time class
    class_<tudat::Time>("astro_time_representation_Time")
        .constructor<double>()
        .constructor<int, long double>()
        .function("toFloat", +[](const tudat::Time& t) { return static_cast<double>(t); });

    // DateTime class
    class_<tba::DateTime>("astro_time_representation_DateTime")
        .constructor<int, int, int, int, int, double>()
        .function("getYear", &tba::DateTime::getYear)
        .function("getMonth", &tba::DateTime::getMonth)
        .function("getDay", &tba::DateTime::getDay)
        .function("getHour", &tba::DateTime::getHour)
        .function("getMinute", &tba::DateTime::getMinute)
        .function("getSeconds", &tba::DateTime::getSeconds)
        .function("epoch", &tba::DateTime::epoch<double>)
        .function("julianDay", &tba::DateTime::julianDay<double>)
        .function("modifiedJulianDay", &tba::DateTime::modifiedJulianDay<double>)
        .function("isoString", &tba::DateTime::isoString)
        .class_function("fromJulianDay", &tba::DateTime::fromJulianDay)
        .class_function("fromModifiedJulianDay", &tba::DateTime::fromModifiedJulianDay);

    // Time conversion functions
    function("astro_time_representation_julian_day_to_seconds_since_epoch",
        &tba::convertJulianDayToSecondsSinceEpoch<double>);

    function("astro_time_representation_seconds_since_epoch_to_julian_day",
        &tba::convertSecondsSinceEpochToJulianDay<double>);

    function("astro_time_representation_julian_day_to_modified_julian_day",
        &tba::convertJulianDayToModifiedJulianDay<double>);

    function("astro_time_representation_modified_julian_day_to_julian_day",
        &tba::convertModifiedJulianDayToJulianDay<double>);

    function("astro_time_representation_is_leap_year",
        &tba::isLeapYear);

    function("astro_time_representation_get_days_in_month",
        &tba::getDaysInMonth);

    function("astro_time_representation_TAI_to_TT",
        &tba::convertTAItoTT<double>);

    function("astro_time_representation_TT_to_TAI",
        &tba::convertTTtoTAI<double>);

    function("astro_time_representation_TT_to_TDB_approximate",
        &tba::approximateConvertTTtoTDB);

    function("astro_time_representation_date_time_to_epoch",
        &tba::convertCalendarDateToJulianDaysSinceEpoch<double>);
}

#endif
