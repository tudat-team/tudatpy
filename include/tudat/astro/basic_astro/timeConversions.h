/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_TIME_CONVERSIONS_H
#define TUDAT_TIME_CONVERSIONS_H

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! List of time scales available
/*!
 *  List of time scales available. Two types of scales are included, earth-based time scales (which can be handled, in part, by SOFA),
 *  which represent a unique specific scale, and three relativistic scales, of which only barycentric is unique. The bodycentric and
 *  topocentric scales require additional identifiers to fully determine.
 */
enum TimeScales
{
    dummy_scale = -1,
    tai_scale = 0,
    tt_scale = 1,
    tdb_scale = 2,
    utc_scale = 3,
    ut1_scale = 4,
    body_centered_coordinate_time_scale = 5,
    barycentric_coordinate_time_scale = 6,
    local_proper_time_scale = 7
};

const static int JULIAN_DAY_ON_J2000_INT = 2451545;

//! Julian day at J2000, i.e. 01-01-2000, at 12:00 (in TT).
const static double JULIAN_DAY_ON_J2000 = mathematical_constants::getFloatingInteger< double >( JULIAN_DAY_ON_J2000_INT );

//! Julian day at J2000, i.e. 01-01-2000, at 12:00 (in TT),  in long double precision.
const static long double JULIAN_DAY_ON_J2000_LONG = mathematical_constants::getFloatingInteger< long double >( JULIAN_DAY_ON_J2000_INT );

//! Function to get the Julian day on J2000
/*!
 *  Function to get the Julian day on J2000, in the requested time representation type
 *  \return  Julian day on J2000
 */
template< typename TimeType >
TimeType getJulianDayOnJ2000( );

template< typename TimeType >
constexpr TimeType getModifiedJulianDayOnJ2000( )
{
    return static_cast< TimeType >( 51544 ) + mathematical_constants::getFloatingFraction< TimeType >( 1, 2 );
}

//! Julian day at Modified Julain Date 0, i.e. Nov 17, 1858, 00:00.
const static double JULIAN_DAY_AT_0_MJD = 2400000.5;

//! Julian day at Modified Julain Date 0, i.e. Nov 17, 1858, 00:00, in long double precision.
const static long double JULIAN_DAY_AT_0_MJD_LONG = 2400000.5L;

//! Function to get the Julian day on zero modified Julian day.
/*!
 *  Function to get the Julian day on on zero modified Julian day, in the requested time representation type
 *  \return Julian day on zero modified Julian day.
 */
template< typename TimeType >
TimeType getJulianDayOnMjd0( );

//! Julian day at which TT, TCG, and TCB all show exact same time (1977 January 1, 00:00:32.184)
const static double TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION = 2443144.5003725;

//! Julian day at which TT, TCG, and TCB all show exact same time (1977 January 1, 00:00:32.184), in long double precision.
const static long double TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION_LONG = 2443144.5003725L;

//! Function to get the synchronization Julian day of TT, TCG, and TCB.
/*!
 *  Function to get the synchronization Julian day of TT, TCG, and TCB, in the requested time representation type
 *  \return Synchronization Julian day of TT, TCG, and TCB.
 */
template< typename TimeType >
TimeType getTimeOfTaiSynchronizationJulianDay( );

//! Function to get the synchronization time of TT, TCG, and TCB, in seconds since J2000.
/*!
 *  Function to get the synchronization time of TT, TCG, and TCB, in seconds since J2000 (constant by definition).
 *  \return Synchronization time of TT, TCG, and TCB, in seconds since J2000.
 */
template< typename TimeType >
TimeType getTimeOfTaiSynchronizationSinceJ2000( )
{
    return ( getTimeOfTaiSynchronizationJulianDay< TimeType >( ) - getJulianDayOnJ2000< TimeType >( ) ) *
             physical_constants::getJulianDay< TimeType >( );
}

//! Difference between TDB and (TT, TCB and TCB) at TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION
const static double TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION = -6.55E-5;

//! Difference between TDB and (TT, TCB and TCB) at TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION, in the requested time type
const static double TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION_LONG = -6.55E-5L;

//! Function to get the difference between TDB and (TT, TCB and TCB) at TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION
/*!
 *  Function to get the difference between TDB and (TT, TCB and TCB) at TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION,
 *  in the requested time representation type
 *  \return Difference between TDB and (TT, TCB and TCB) at TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION (constant by definition).
 */
template< typename TimeType >
TimeType getTdbSecondsOffsetAtSynchronization( );

//! Offset of TT from TAI (constant by definition).
const static double TT_MINUS_TAI = 32.184;

//! Offset of TT from TAI (constant by definition), in long double precision.
const static double TT_MINUS_TAI_LONG = 32.184L;

//! Function to get the offset of TT from TAI (constant by definition)
/*!
 *  Function to get the offset of TT from TAI (constant by definition), in the requested time representation type
 *  \return Offset of TT from TAI (constant by definition).
 */
template< typename TimeType >
TimeType getTTMinusTai( );

//! Compute number of seconds since a reference Julian day.
/*!
 * Computes the number of seconds since a reference Julian day from
 * Julian day.
 * \param julianDay Julian day at which number of seconds since epoch is to be determined.
 * \param epochSinceJulianDayZero Epoch in Julian day.
 * \return Number of seconds since epoch.
 */
template< typename TimeScalarType = double >
TimeScalarType convertJulianDayToSecondsSinceEpoch(
        const TimeScalarType julianDay,
        const TimeScalarType epochSinceJulianDayZero = getJulianDayOnJ2000< TimeScalarType >( ) )
{
    return ( julianDay - epochSinceJulianDayZero ) * physical_constants::getJulianDay< TimeScalarType >( );
}

//! Compute Julian day from seconds since reference Julian day epoch.
/*!
 * Computes the Julian day bsaed on seconds since reference Julian day epoch provided.
 * \param secondsSinceEpoch Seconds since epoch.
 * \param epochSinceJulianDayZero Epoch in Julian day.
 * \return Number of Julian days since epoch.
 */
template< typename TimeScalarType = double >
TimeScalarType convertSecondsSinceEpochToJulianDay(
        const TimeScalarType secondsSinceEpoch,
        const TimeScalarType epochSinceJulianDayZero = getJulianDayOnJ2000< TimeScalarType >( ) )
{
    return ( secondsSinceEpoch / physical_constants::getJulianDay< TimeScalarType >( ) + epochSinceJulianDayZero );
}

inline int julianDayNumberFromDate(
    int year, int month, int day )
{
    int a, m, y, leap_days;
    a = (14 - month) / 12;
    m = (month - 3) + (12 * a);
    y = year + 4800 - a;
    leap_days = (y / 4) - (y / 100) + (y / 400);
    return day + (((153 * m) + 2) / 5) + (365 * y) + leap_days - 32045;

}

//! Compute Julian day from given date and time
/*!
  * Computes the Julian day from given year, month, day, hour, minutes, seconds as used in everyday
  * life. The function uses the internal calcualtions of the boost::date_time::gregorian class.
  *
  * \param calendarYear Year of the standard calendar in years.
  * \param calendarMonth Month of the standard calendar in months.
  * \param calendarDay Day of the standard calendar in days.
  * \param calendarHour Hour of the time of this day in hours.
  * \param calendarMinutes Minutes of the time of this day in minutes.
  * \param calendarSeconds Seconds of the time of this day in seconds.
  * \param referenceJulianDay Reference Julian day (i.e. t=0) for input time.
  * \return Days since referenceJulianDay for input date/time.
  */
template< typename TimeScalarType = double >
TimeScalarType convertCalendarDateToJulianDaysSinceEpoch( const int calendarYear,
                                                          const int calendarMonth,
                                                          const int calendarDay,
                                                          const int calendarHour,
                                                          const int calendarMinutes,
                                                          const TimeScalarType calendarSeconds,
                                                          const TimeScalarType referenceJulianDay )
{
    // Calculate julian day of calendar date.
    TimeScalarType julianDay =
            static_cast< TimeScalarType >( julianDayNumberFromDate( calendarYear, calendarMonth, calendarDay ) - referenceJulianDay );

    //Compute day fraction
    const TimeScalarType dayFraction =
            static_cast< TimeScalarType >( calendarHour ) /
            mathematical_constants::getFloatingInteger< TimeScalarType >( 24 )  +
            static_cast< TimeScalarType >( calendarMinutes ) /
            mathematical_constants::getFloatingInteger< TimeScalarType >( 24 * 60 ) +
            calendarSeconds / mathematical_constants::getFloatingInteger< TimeScalarType >( 24 * 3600 );


    // Compute Julian day by adding day fraction and subtracting 0.5 to reference to midnight instead of noon..
    return julianDay + dayFraction - mathematical_constants::getFloatingFraction< TimeScalarType >( 1, 2 );
}

//! Compute Julian day from given date and time
/*!
  * Computes the Julian day from given year, month, day, hour, minutes, seconds as used in everyday
  * life. The function uses the internal calcualtions of the boost::date_time::gregorian class.
  *
  * \param calendarYear Year of the standard calendar in years.
  * \param calendarMonth Month of the standard calendar in months.
  * \param calendarDay Day of the standard calendar in days.
  * \param calendarHour Hour of the time of this day in hours.
  * \param calendarMinutes Minutes of the time of this day in minutes.
  * \param calendarSeconds Seconds of the time of this day in seconds.
  */
template< typename TimeScalarType = double >
TimeScalarType convertCalendarDateToJulianDay( const int calendarYear,
                                               const int calendarMonth,
                                               const int calendarDay,
                                               const int calendarHour,
                                               const int calendarMinutes,
                                               const TimeScalarType calendarSeconds )
{
    return convertCalendarDateToJulianDaysSinceEpoch< TimeScalarType >
            ( calendarYear, calendarMonth, calendarDay, calendarHour, calendarMinutes,
              calendarSeconds, mathematical_constants::getFloatingInteger< TimeScalarType >( 0.0 ) );
}

//! Compute Julian day since J2000 from given date and time
/*!
  * Computes the Julian day from given year, month, day, hour, minutes, seconds as used in everyday
  * life. The function uses the internal calcualtions of the boost::date_time::gregorian class.
  *
  * \param calendarYear Year of the standard calendar in years.
  * \param calendarMonth Month of the standard calendar in months.
  * \param calendarDay Day of the standard calendar in days.
  * \param calendarHour Hour of the time of this day in hours.
  * \param calendarMinutes Minutes of the time of this day in minutes.
  * \param calendarSeconds Seconds of the time of this day in seconds.
  */
template< typename TimeScalarType = double >
TimeScalarType convertCalendarDateToJulianDaySinceJ2000(const int calendarYear,
                                                        const int calendarMonth,
                                                        const int calendarDay,
                                                        const int calendarHour,
                                                        const int calendarMinutes,
                                                        const TimeScalarType calendarSeconds)
{
  return basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch<TimeScalarType>
      (calendarYear, calendarMonth, calendarDay, calendarHour, calendarMinutes,
       calendarSeconds, basic_astrodynamics::JULIAN_DAY_ON_J2000);
}

//! Function to convert julian day to modified julian day.
/*!
 *  Function to convert julian day to modified julian day.
 *  \param julianDay Julian day to convert
 *  \return Modified julian day as obtained from julian day.
 */
template< typename TimeScalarType >
TimeScalarType convertJulianDayToModifiedJulianDay( const TimeScalarType julianDay )
{
    return julianDay - getJulianDayOnMjd0< TimeScalarType >( );
}

//! Function to convert modified julian dat to julian day.
/*!
 *  Function to convert modified julian dat to julian day.
 *  \param modifiedJulianDay Modified julian day to convert
 *  \return Julian day as obtained from modified julian day.
 */
template< typename TimeScalarType >
TimeScalarType convertModifiedJulianDayToJulianDay( const TimeScalarType modifiedJulianDay )
{
    return modifiedJulianDay + getJulianDayOnMjd0< TimeScalarType >( );
}

//! Function to convert the number of seconds since some reference julian day to a julian day since that epoch.
/*!
 *  Function to convert the number of seconds since some reference julian day to a julian year since that epoch.
 *  \param secondsSinceEpoch Seconds since some reference epoch to convert to julian year.
 *  \return Number of julian years since epoch
 */
template< typename TimeScalarType = double >
TimeScalarType convertSecondsSinceEpochToJulianYearsSinceEpoch( const TimeScalarType secondsSinceEpoch )
{
    return secondsSinceEpoch / ( physical_constants::getJulianDay< TimeScalarType >( ) *
                                 physical_constants::getJulianYearInDays< TimeScalarType >( ) );
}

//! Function to convert the number of seconds since some reference julian day to a julian century since that epoch.
/*!
 *  Function to convert the number of seconds since some reference julian day to a julian century since that epoch.
 *  \param secondsSinceEpoch Seconds since some reference epoch to convert to julian year.
 *  \return Number of julian centuries since epoch
 */
template< typename TimeScalarType = double >
TimeScalarType convertSecondsSinceEpochToJulianCenturiesSinceEpoch( const TimeScalarType secondsSinceEpoch )
{
    return convertSecondsSinceEpochToJulianYearsSinceEpoch( secondsSinceEpoch ) / 100.0;
}


template< typename ScalarType >
void convertShiftedJulianDayToCalendarDate(
    ScalarType shiftedJulianDay, int& day, int& month, int& year)
{
    // Declare temporary variables.
    int L, M, N, P, Q;

    // Execute algorithm.
    if( shiftedJulianDay > 2299160 )    // after Oct 4, 1582
    {
        L = shiftedJulianDay + 68569;
        M = (4 * L) / 146097;
        L = L - ((146097 * M + 3) / 4);
        N = (4000 * (L + 1)) / 1461001;
        L = L - ((1461 * N) / 4) + 31;
        P = (80 * L) / 2447;
        day = int(L - (2447 * P) / 80);
        L = P / 11;
        month = int(P + 2 - 12 * L);
        year = int(100 * (M - 49) + N + L);
    }
    else
    {
        P = shiftedJulianDay + 1402;
        Q = (P - 1) / 1461;
        L = P - 1461 * Q;
        M = (L - 1) / 365 - L / 1461;
        N = L - 365 * M + 30;
        P = (80 * N) / 2447;
        day = int(N - (2447 * P) / 80);
        N = P / 11;
        month = int( P + 2 - 12 * N );
        year = int( 4 * Q + M + N - 4716 );
        if(year <= 0)
        {
            --year;
        }
    }
    // catch century/non-400 non-leap years
    if( year > 1599 && !( year % 100 )
        && ( year % 400 ) && month == 2 && day == 29 )
    {
        month = 3;
        day = 1;
    }
}


template< typename ScalarType >
void convertJulianDayToCalendarDate( ScalarType julianDay,
                                     int& day, int& month, int& year)
{
    // Execute algorithm.
    ScalarType shiftedJulianDay = julianDay + mathematical_constants::getFloatingFraction< ScalarType >( 1, 2 );
    convertShiftedJulianDayToCalendarDate( shiftedJulianDay, day, month, year );
}

template< typename ScalarType >
void convertModifiedJulianDayToCalendarDate(
    ScalarType modifiedJulianDay,
    int& day, int& month, int& year)
{
    // Execute algorithm.
    ScalarType shiftedJulianDay = modifiedJulianDay + mathematical_constants::getFloatingInteger< ScalarType >( 2400001 );
    convertShiftedJulianDayToCalendarDate( shiftedJulianDay, day, month, year );
}

//! Function to convert julian day to gregorian calendar date.
/*!
 *  Function to convert julian day to gregorian calendar date. Algorithm from (Wertz, 2009, Eq. 4.3)
 *  \param julianDay Julian day to convert to gregorian calendar date.
 *  \return Gregorian calendar date at given julian day.
 */
//! Function to convert julian day to gregorian calendar date.
template< typename ScalarType >
boost::gregorian::date convertJulianDayToCalendarDate( ScalarType julianDay )
{

    // Declare date variables.
    int day, month, year;
    convertJulianDayToCalendarDate( julianDay, day, month, year );

    return boost::gregorian::date( year, month, static_cast< double >( day ) );
}

//! Function to convert gregorian calendar date and fraction of day to a number of julian days since a reference epoch.
/*!
 *  Function to convert gregorian calendar date and fraction of day to a number of julian days since a reference epoch.
 *  \param calendarDate Calendar date at which conversion is to be performed.
 *  \param fractionOfDay Fraction of day in current gregorian date.
 *  \param epochSinceJulianDayZero Reference epoch since which number of julian days is to be determined.
 *  \return Number of julian days since a reference epoch.
 */
template< typename TimeScalarType >
TimeScalarType calculateJulianDaySinceEpoch( const boost::gregorian::date calendarDate,
                                             const TimeScalarType fractionOfDay,
                                             const TimeScalarType epochSinceJulianDayZero = JULIAN_DAY_ON_J2000 )
{
    TimeScalarType julianDayAtMiddleOfDay = static_cast< TimeScalarType >( calendarDate.julian_day( ) );
    return julianDayAtMiddleOfDay +
            ( fractionOfDay - mathematical_constants::getFloatingFraction< TimeScalarType >( 1, 2 ) ) -
            epochSinceJulianDayZero;
}

//! Function to determine whether the given year is a leap year (i.e. has 366 days)
/*!
 *  Function to determine whether the given year is a leap year (i.e. has 366 days)
 *  \param year Year for which it is to be determined whether it is a leap year.
 *  \return True if year is a leap year, false otherwise.
 */
bool isLeapYear( const int year );

//! Function that returns number of days in given month number
/*!
 *  Function that returns number of days in given month number (January = 1, December = 12)
 *  \param month Month number
 *  \param year Year
 *  \return Number of days in given month
 */
int getDaysInMonth( const int month, const int year );

//! Determine number of full days that have passed in current year
/*!
 *  Determine number of full days that have passed in current year (i.e. 0 to 364 or 0 to 365 for leap year)
 *  This function creates a boost::gregorian::date object and calls the overloaded function.
 *  \param day Day in current year (1 is first day)
 *  \param month Month in current year (1 is first month)
 *  \param year Current year (for leap year purposes)
 *  \return Number of full days that have passed in current year
 */
int convertDayMonthYearToDayOfYear( const int day,
                                    const int month,
                                    const int year );

//! Determine number of full days that have passed in current year
/*!
 *  Determine number of full days that have passed in current year (i.e. 0 to 364 or 0 to 365 for leap year)
 *  \param calendarDate Gregorian calendat date object containing current date information.
 *  \return Number of full days that have passed in current year
 */
int convertDayMonthYearToDayOfYear( const boost::gregorian::date calendarDate );

//! Determine number of seconds into current day of given time
/*!
 *  Determine number of seconds into current day of given time
 *  \param julianDay Current julian day
 *  \return Number of seconds into current day.
 */
double calculateSecondsInCurrentJulianDay( const double julianDay );

//! Function to create the calendar date from the year and the number of days in the year
/*!
 *  Function to create the calendar date from the year and the number of days in the year, where Jan. 1 is day in year 0.
 *  \param year Year in which days in year are given
 *  \param daysInYear Number of days in current year, note that day must be in current year (i.e <= 365 of 366 for leap year)
 *  \return Gregorian date object as generated from input
 */
boost::gregorian::date convertYearAndDaysInYearToDate( const int year, const int daysInYear );

//! Function to convert TCB to TDB times scale
/*!
 *  Function to convert TCB to TDB times scale, with both input and output referenced to the J2000 reference time.
 *  \param tcbTime Input time in TCB scale
 *  \return Converted time in TDB scale
 */
template< typename TimeType >
TimeType convertTcbToTdb( const TimeType tcbTime )
{
    return tcbTime + getTdbSecondsOffsetAtSynchronization< TimeType >( ) -
            physical_constants::getLbTimeRateTerm< TimeType >( ) *
            ( tcbTime - getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) );
}

//! Function to convert TDB to TCB times scale
/*!
 *  Function to convert TDB to TCB times scale, with both input and output referenced to the J2000 reference time.
 *  \param tdbTime Input time in TDB scale
 *  \return Converted time in TCB scale
 */
template< typename TimeType >
TimeType convertTdbToTcb( const TimeType tdbTime )
{
    return ( tdbTime - getTdbSecondsOffsetAtSynchronization< TimeType >( ) -
             physical_constants::getLbTimeRateTerm< TimeType >( ) *
             getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) ) /
            ( mathematical_constants::getFloatingInteger< TimeType >( 1 ) -
              physical_constants::getLbTimeRateTerm< TimeType >( ) );
}

//! Function to convert TCG to TT times scale
/*!
 *  Function to convert TCG to TT times scale, with both input and output referenced to the J2000 reference time.
 *  \param tcgTime Input time in TCG scale
 *  \return Converted time in TT scale
 */
template< typename TimeType >
TimeType convertTcgToTt( const TimeType tcgTime  )
{
    return tcgTime - physical_constants::getLgTimeRateTerm< TimeType >( ) *
            ( tcgTime - getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) );
}

//! Function to convert TT to TCG times scale
/*!
 *  Function to convert TT to TCG times scale, with both input and output referenced to the J2000 reference time.
 *  \param ttTime Input time in TT scale
 *  \return Converted time in TCG scale
 */
template< typename TimeType >
TimeType convertTtToTcg( const TimeType ttTime  )
{
    return ttTime + physical_constants::getLgTimeRateTerm< TimeType >( )
            / ( mathematical_constants::getFloatingInteger< TimeType >( 1 ) -
                physical_constants::getLgTimeRateTerm< TimeType >( ) ) *
            ( ttTime - getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) );
}

//! Function to convert TAI to TT
/*!
 *  Function to convert TAI (International Atomic Time) to TT (Terrestrial Time) by adding bias as defined by Sofa.
 *  \param taiTime Time in TAI
 *  \return Time in TT
 */
template< typename TimeType >
TimeType convertTAItoTT( const TimeType taiTime )
{
    return taiTime + getTTMinusTai< TimeType >( );
}


//! Function to convert TT to TAI
/*!
 *  Function to convert TT (Terrestrial Time) to TAI (International Atomic Time) by subtracting bias as defined by Sofa.
 *  \param ttTime Time in TT
 *  \return Time in TAI
 */
template< typename TimeType >
TimeType convertTTtoTAI( const TimeType ttTime )
{
    return ttTime - getTTMinusTai< TimeType >( );
}


//! Perform apprixmate conversion of TT to TDB
/*!
 * Perform apprixmate conversion of TT to TDB, in which only the once-per-orbit sinusoidal effect of O(e) is taken into
 * account. Accurate conversions are calculated usinf Sofa
 * \param ttSecondsSinceJ2000 Terrestrial time in seconds since J2000
 * \return TDB in seconds since J2000
 */
double approximateConvertTTtoTDB( const double ttSecondsSinceJ2000);



} // namespace basic_astrodynamics
} // tudat

#endif // TUDAT_TIME_CONVERSIONS_H
