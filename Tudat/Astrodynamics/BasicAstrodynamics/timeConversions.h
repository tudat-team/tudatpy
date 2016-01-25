/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130218    D. Dirkx          File created from personal application.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_TIME_CONVERSIONS_H
#define TUDAT_TIME_CONVERSIONS_H

#include "boost/date_time/gregorian/gregorian.hpp"

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

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

//! Julian day at J2000.
/*!
 * Julian day at J2000, i.e. 01-01-2000, at 12:00 (in TT).
 */
const static double JULIAN_DAY_ON_J2000 = 2451545.0;

const static long double JULIAN_DAY_ON_J2000_LONG = static_cast< long double >( 2451545.0 );

template< typename TimeType >
TimeType getJulianDayOnJ2000( );

//! Julian day at Modified Julain Date 0.
/*!
 * Julian day at Modified Julain Date 0, i.e. Nov 17, 1858, 00:00.
 */
const static double JULIAN_DAY_AT_0_MJD = 2400000.5;

const static long double JULIAN_DAY_AT_0_MJD_LONG = static_cast< long double >( 2400000.5 );

template< typename TimeType >
TimeType getJulianDayOnMjd0( );

const static double TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION = 2443144.5003725;

const static long double TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION_LONG = static_cast< long double >( 2443144.5003725 );

template< typename TimeType >
TimeType getTimeOfTaiSynchronizationJulianDay( );

const static double TAI_JULIAN_DAY_SINCE_J2000_AT_TIME_SYNCHRONIZATION = -8.4004996275E3;

const static long double TAI_JULIAN_DAY_SINCE_J2000_AT_TIME_SYNCHRONIZATION_LONG = static_cast< long double >( -8400.4996275 );

template< typename TimeType >
TimeType getTimeOfTaiSynchronizationSinceJ2000( );

const static double TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION = -6.55E-5;

const static double TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION_LONG = static_cast< long double >( -6.55E-5 );

//! Compute number of seconds since a reference Julian day.
/*!
 * Computes the number of seconds since a reference Julian day from
 * Julian day.
 * \param julianDay Julian day at which number of seconds since epoch is to be determined.
 * \param epochSinceJulianDayZero Epoch in Julian day.
 * \return Number of seconds since epoch.
 */
double convertJulianDayToSecondsSinceEpoch( const double julianDay,
                                            const double epochSinceJulianDayZero );

//! Compute Julian day from seconds since reference Julian day epoch.
/*!
 * Computes the Julian day bsaed on seconds since reference Julian day epoch provided.
 * \param secondsSinceEpoch Seconds since epoch.
 * \param epochSinceJulianDayZero Epoch in Julian day.
 * \return Number of Julian days since epoch.
 */
double convertSecondsSinceEpochToJulianDay(
                                      const double secondsSinceEpoch,
                                      const double epochSinceJulianDayZero = JULIAN_DAY_ON_J2000 );

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
double convertCalendarDateToJulianDay( const int calendarYear,
                                       const int calendarMonth,
                                       const int calendarDay,
                                       const int calendarHour,
                                       const int calendarMinutes,
                                       const double calendarSeconds );

//! Function to determine whether a time scale is a general, relativistic time scale.
/*!
 *  Function to determine whether a time scale is a general, relativistic time scale. Available time scales are barycentric (TCB),
 *  bodycentric, of which TCG is a specific (Earth-centered) case, and topocentric, which represents the proper time at a
 *  body-fixed point.
 *  \param timeScale Time scale which is to be checked.
 *  \return True if timeScale is one of the general, relativistic scales, false otherwise.
 */
bool isTimeScaleRelativistic( const TimeScales timeScale );

//! Function to determine whether a time scale is one of the specific Earth-based scales
/*!
 *  Function to determine whether a time scale is one of the specific Earth-based scales, which are (in part) handled by sofa.
 *  These are the UTC and UT1 sidereal times, the TT and TAI atomic times and the TDB (rescaled TCB)
 *  \param timeScale Time scale which is to be checked.
 *  \return True if timeScale is one of the specific Earth-based time scales, false otherwise.
 */
bool isTimeScaleFromSofa( const TimeScales timeScale );

namespace basic_astrodynamics
{


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
 *  \param julianDay Modified julian day to convert
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
double convertSecondsSinceEpochToJulianYearsSinceEpoch( const double secondsSinceEpoch );

//! Function to convert the number of seconds since some reference julian day to a julian century since that epoch.
/*!
 *  Function to convert the number of seconds since some reference julian day to a julian century since that epoch.
 *  \param secondsSinceEpoch Seconds since some reference epoch to convert to julian year.
 *  \return Number of julian centuries since epoch
 */
double convertSecondsSinceEpochToJulianCenturiesSinceEpoch( const double secondsSinceEpoch );

//! Function to convert the number of seconds since a reference julian day to a number of full julian days + fraction of day since that epoch.
/*!
 *  Function to convert the number of seconds since a reference julian day to a number of full julian days + fraction of day since that epoch.
 *  Values are returned by passing to this function by reference.
 *  \param secondsSinceEpoch Seconds since epoch to convert to number of full julian days + fraction of day since that epoch.
 *  \param fullDays Full julian days since epochs (return by reference)
 *  \param fractionOfDay Fraction of day in current day since reference epoch (return by reference)
 *  \return Number of julian days since epoch
 */
void convertSecondsSinceEpochToFullJulianDaysAndFractionOfDay(
        const double secondsSinceEpoch, double& fullDays, double& fractionOfDay );

//! Function to convert julian day to gregorian calendar date.
/*!
 *  Function to convert julian day to gregorian calendar date. Algorithm from (Wertz, 2009, Eq. 4.3)
 *  \param julianDay Julian day to convert to gregorian calendar date.
 *  \return Gregorian calendar date at given julian day.
 */
boost::gregorian::date convertJulianDayToCalendarDate( const double julianDay );

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
    return julianDayAtMiddleOfDay + ( fractionOfDay - mathematical_constants< TimeScalarType >( ) ) -
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

//! Calculates the number of full days between 01-Jan of two given years.
/*!
 *  Calculates the number of full days between 01-Jan of two given years. Start year may be smaller than end year, in which case a negative
 *  result is produced. If the two input years are equal, result is zero.
 *  \param startYear Year from which to start counting number of days
 *  \param endYear Year up to which to start counting number of days.
 *  \param Number of days from first to second year (01 Jan)
 */
int calculateDaysBetweenYears( int startYear, int endYear );

//! Function to create the calendar date from the year and the number of days in the year
/*!
 *  Function to create the calendar date from the year and the number of days in the year, where Jan. 1 is day in year 0.
 *  \param year Year in which days in year are given
 *  \param daysInYear Number of days in current year, note that day must be in current year (i.e <= 365 of 366 for leap year)
 *  \return Gregorian date object as generated from input
 */
boost::gregorian::date convertYearAndDaysInYearToDate( const int year, const int daysInYear );

template< typename TimeType >
TimeType doDummyTimeConversion( const TimeType inputTime )
{
    return inputTime;
}

//! Function to convert TCB to TDB times scale
/*!
 *  Function to convert TCB to TDB times scale, with both input and output referenced to the J2000 reference time.
 *  \param inputTime Input time in TCB scale
 *  \param Converted time in TDB scale
 */
template< typename TimeType >
TimeType convertTcbToTdb( const TimeType inputTime )
{
    return inputTime + TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION - physical_constants::LB_TIME_RATE_TERM_LONG *
            ( inputTime - getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) );
}

//! Function to convert TDB to TCB times scale
/*!
 *  Function to convert TDB to TCB times scale, with both input and output referenced to the J2000 reference time.
 *  \param inputTime Input time in TDB scale
 *  \param Converted time in TCB scale
 */
template< typename TimeType >
TimeType convertTdbToTcb( const TimeType inputTime )
{
    return ( inputTime - TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION - physical_constants::LB_TIME_RATE_TERM_LONG *
             getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) ) /
            ( 1.0L - physical_constants::LB_TIME_RATE_TERM_LONG );
}

//! Function to convert TCG to TT times scale
/*!
 *  Function to convert TCG to TT times scale, with both input and output referenced to the J2000 reference time.
 *  \param inputTime Input time in TCG scale
 *  \param Converted time in TT scale
 */
template< typename TimeType >
TimeType convertTcgToTt( const TimeType inputTime  )
{
    return ( inputTime + physical_constants::LG_TIME_RATE_TERM_LONG *
             getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) ) /
            ( 1.0L + physical_constants::LG_TIME_RATE_TERM_LONG );
}

//! Function to convert TT to TCG times scale
/*!
 *  Function to convert TT to TCG times scale, with both input and output referenced to the J2000 reference time.
 *  \param inputTime Input time in TT scale
 *  \param Converted time in TCG scale
 */
template< typename TimeType >
TimeType convertTtToTcg( const TimeType inputTime  )
{
    return inputTime + physical_constants::LG_TIME_RATE_TERM_LONG *
            ( inputTime - getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) );
}

//! Function to calculate the difference between TCG and TT time scales from TT input
/*!
 *  Function to calculate the difference between TCG and TT time scales from TT input, with input value referenced to the J2000 reference time.
 *  \param terrestrialTime Input time in TT scale
 *  \param TCG minus TT at requested TT
 */
template< typename TimeType >
TimeType getTcgMinusTT( const TimeType terrestrialTime )
{
    return physical_constants::LG_TIME_RATE_TERM_LONG / ( 1.0L - physical_constants::LG_TIME_RATE_TERM_LONG ) * (
                terrestrialTime - getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) );
}

//! Function to calculate the difference between TCB and TDB time scales from TDB input
/*!
 *  Function to calculate the difference between TCB and TDB time scales from TDB input, with input value referenced to the J2000 reference time.
 *  \param tdbTime Input time in TDB scale
 *  \param TCB minus TDB at requested TDB
 */
template< typename TimeType >
TimeType getTcbMinusTdb( const TimeType tdbTime )
{
    return ( physical_constants::LB_TIME_RATE_TERM_LONG * (
                 tdbTime - getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) ) -
             TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION ) / ( 1.0L - physical_constants::LB_TIME_RATE_TERM_LONG );
}

//! Function to calculate the difference between TDB and TCB time scales from TCB input
/*!
 *  Function to calculate the difference between TDB and TCB time scales from TCB input, with input value referenced to the J2000 reference time.
 *  \param tdbTime Input time in TCB scale
 *  \param TDB minus TCB at requested TCB
 */
template< typename TimeType >
TimeType getTdbMinusTcb( const TimeType tcbTime )
{
    return TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION - physical_constants::LB_TIME_RATE_TERM_LONG *
            ( tcbTime - getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) );
}

} // namespace basic_astrodynamics
} // tudat

#endif // TUDAT_TIME_CONVERSIONS_H
