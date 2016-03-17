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


//! Julian day at J2000, i.e. 01-01-2000, at 12:00 (in TT).
const static double JULIAN_DAY_ON_J2000 = 2451545.0;

//! Julian day at J2000, i.e. 01-01-2000, at 12:00 (in TT),  in long double precision.
const static long double JULIAN_DAY_ON_J2000_LONG = 2451545.0L;

//! Function to get the Julian day on J2000
/*!
 *  Function to get the Julian day on J2000, in the requested time representation type
 *  \return  Julian day on J2000
 */
template< typename TimeType >
TimeType getJulianDayOnJ2000( );

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
            static_cast< TimeScalarType >( boost::gregorian::date( calendarYear, calendarMonth, calendarDay ).
                                           julian_day( ) ) - referenceJulianDay;

    //Compute day fraction
    const TimeScalarType dayFraction =
            static_cast< TimeScalarType >( calendarHour ) / static_cast< TimeScalarType >( 24.0 ) +
            static_cast< TimeScalarType >( calendarMinutes ) / static_cast< TimeScalarType >( 24.0 * 60.0 ) +
            calendarSeconds / static_cast< TimeScalarType >( 24.0 * 3600.0 );

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
    return secondsSinceEpoch / ( physical_constants::JULIAN_DAY * physical_constants::JULIAN_YEAR_IN_DAYS );
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

template< typename TimeType >
TimeType doDummyTimeConversion( const TimeType inputTime )
{
    return inputTime;
}

//! Function to convert TCB to TDB times scale
/*!
 *  Function to convert TCB to TDB times scale, with both input and output referenced to the J2000 reference time.
 *  \param inputTime Input time in TCB scale
 *  \return Converted time in TDB scale
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
 *  \return Converted time in TCB scale
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
 *  \return Converted time in TT scale
 */
template< typename TimeType >
TimeType convertTcgToTt( const TimeType tcgTime  )
{
    return tcgTime - physical_constants::LG_TIME_RATE_TERM_LONG *
            ( tcgTime - getTimeOfTaiSynchronizationSinceJ2000< TimeType >( ) );
}

//! Function to convert TT to TCG times scale
/*!
 *  Function to convert TT to TCG times scale, with both input and output referenced to the J2000 reference time.
 *  \param inputTime Input time in TT scale
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
 *  \param tai Time in TAI
 *  \return Time in TT
 */
template< typename TimeType >
TimeType convertTAItoTT( const TimeType tai )
{
    return tai + getTTMinusTai< TimeType >( );
}


//! Function to convert TT to TAI
/*!
 *  Function to convert TT (Terrestrial Time) to TAI (International Atomic Time) by subtracting bias as defined by Sofa.
 *  \param tt Time in TT
 *  \return Time in TAI
 */
template< typename TimeType >
TimeType convertTTtoTAI( const TimeType tt )
{
    return tt - getTTMinusTai< TimeType >( );
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
