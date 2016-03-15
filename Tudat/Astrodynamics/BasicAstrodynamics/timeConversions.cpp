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

#include <boost/date_time/gregorian/gregorian.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{
namespace basic_astrodynamics
{

template< >
double getJulianDayOnJ2000< double >( )
{
    return JULIAN_DAY_ON_J2000;
}

template< >
long double getJulianDayOnJ2000< long double >( )
{
    return JULIAN_DAY_ON_J2000_LONG;
}

template< >
double getJulianDayOnMjd0< double >( )
{
    return JULIAN_DAY_AT_0_MJD;
}

template< >
long double getJulianDayOnMjd0< long double >( )
{
    return JULIAN_DAY_AT_0_MJD_LONG;
}

template< >
double getTimeOfTaiSynchronizationJulianDay< double >( )
{
    return TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION;
}

template< >
long double getTimeOfTaiSynchronizationJulianDay< long double >( )
{
    return TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION_LONG;
}

template< >
double getTimeOfTaiSynchronizationSinceJ2000< double >( )
{
    return TAI_JULIAN_DAY_SINCE_J2000_AT_TIME_SYNCHRONIZATION * physical_constants::getJulianDay< double >( );
}

template< >
long double getTimeOfTaiSynchronizationSinceJ2000< long double >( )
{
    return TAI_JULIAN_DAY_SINCE_J2000_AT_TIME_SYNCHRONIZATION_LONG * physical_constants::getJulianDay< long double >( );
}

//! Function to convert julian day to gregorian calendar date.
boost::gregorian::date convertJulianDayToCalendarDate( const double julianDay )
{
    // Declare temporary variables.
    int L, M, N, P, Q;

    // Declare date variables.
    int day, month, year;

    // Execute algorithm.
    double shiftedJulianDay = julianDay + 0.5;
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

    // Create and return date object

    return boost::gregorian::date( year, month, static_cast< double >( day ) );
}

//! Compute the Julian day from the calendar date and time.
double convertCalendarDateToJulianDay( const int calendarYear,
                                       const int calendarMonth,
                                       const int calendarDay,
                                       const int calendarHour,
                                       const int calendarMinutes,
                                       const double calendarSeconds )
{
    return convertCalendarDateToJulianDaysSinceEpoch( calendarYear, calendarMonth, calendarDay, calendarHour, calendarMinutes,
                                           calendarSeconds, 0.0 );
}

double convertCalendarDateToJulianDaysSinceEpoch( const int calendarYear,
                                                  const int calendarMonth,
                                                  const int calendarDay,
                                                  const int calendarHour,
                                                  const int calendarMinutes,
                                                  const double calendarSeconds,
                                                  const double referenceJulianDay )
{
    // Calculate julian day of calendar date.
    double julianDay =
            static_cast< double >( boost::gregorian::date( calendarYear, calendarMonth, calendarDay ).julian_day( ) ) -
            referenceJulianDay;

    //Compute day fraction
    const double dayFraction = static_cast< double >( calendarHour ) / 24.0 +
                               static_cast< double >( calendarMinutes ) / ( 24.0 * 60.0 ) +
                               calendarSeconds / ( 24.0 * 3600.0 );

    // Compute Julian day by adding day fraction and subtracting 0.5 to reference to midnight
    // instead of noon..
    return julianDay + dayFraction - 0.5;
}


//! Function to determine whether the given year is a leap year (i.e. has 366 days)
bool isLeapYear( const int year )
{
    bool isLeapYear = 0;
    if( ( year % 4 == 0 ) && !( ( year % 100 == 0 ) && !( year % 400 == 0 ) ) )
    {
        isLeapYear = 1;
    }
    return isLeapYear;
}

//! Function that returns number of days in given month number
int getDaysInMonth( const int month,
                    const int year )
{
    // Declare number of days per month
    static const int daysPerMonth[ 12 ] =
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    int numberOfDays = 0;

    // Check input consistency
    if( month < 1 || month > 12 )
    {
        std::cerr<<"Month number "<<month<<" does not exist, value must be gretaer than 0 and smaller than 13"<<std::endl;
    }
    else
    {
        numberOfDays = daysPerMonth[ month - 1 ];

        // Check for leap year.
        if( month == 2 && isLeapYear( year ) )
        {
            numberOfDays++;
        }
    }
    return numberOfDays;
}

//! Determine number of full days that have passed in current year
int convertDayMonthYearToDayOfYear( const int day,
                                    const int month,
                                    const int year )
{
    return convertDayMonthYearToDayOfYear( boost::gregorian::date( year, month, day ) );
}

//! Determine number of full days that have passed in current year
int convertDayMonthYearToDayOfYear( const boost::gregorian::date calendarDate )
{
    // Calculate and return result (taking into account that this function should return 0 for first day, not 1)
    return ( calendarDate.day_of_year( ) - 1 );
}

//! Determine number of seconds into current day of given time
double calculateSecondsInCurrentJulianDay( const double julianDay )
{
    // Calculate and return result, taking into accoun the fact that julian day 0.0 is at 12:00 and result starts counting at 00:00
    return ( julianDay + 0.5 - std::floor( julianDay + 0.5 ) ) * physical_constants::JULIAN_DAY;
}

//! Function to create the calendar date from the year and the number of days in the year
boost::gregorian::date convertYearAndDaysInYearToDate( const int year, const int daysInYear )
{
    // Go to day 1 at 01-01 convention
    int daysLeft = daysInYear + 1;

    // Loop over all months (starting at January until month in which daysinYear is situated is found)
    int currentMonth = 1;
    int daysInCurrentMonth;
    bool isConverged = 0;
    while( !isConverged )
    {
        // Determine days in current month
        daysInCurrentMonth = getDaysInMonth( currentMonth, year );
        if( daysInCurrentMonth < daysLeft )
        {
            // Subtract days in current month from days left.
            daysLeft -= daysInCurrentMonth;

            // Update month and check consistency
            currentMonth++;
            if( currentMonth > 12 )
            {
                std::cerr<<"Error when converting year and days in year to date, month number has exceeded 12"<<std::endl;
            }
        }
        else
        {
            isConverged = 1;
        }
    }

    // Create date and return
    boost::gregorian::date date( year, currentMonth, daysLeft );

    if( date.day_of_year( ) != daysInYear + 1 )
    {
        std::cerr<<"Error when converting year and days in year to date, inconsistent output"<<std::endl;
    }

    return date;
}

} // namespace basic_astrodynamics
} // namespace tudat
