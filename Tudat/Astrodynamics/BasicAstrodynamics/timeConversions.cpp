/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <boost/date_time/gregorian/gregorian.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Basics/timeType.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Function to get the Julian day on J2000, in double precision.
template< >
double getJulianDayOnJ2000< double >( )
{
    return JULIAN_DAY_ON_J2000;
}

//! Function to get the Julian day on J2000, in long double precision.
template< >
long double getJulianDayOnJ2000< long double >( )
{
    return JULIAN_DAY_ON_J2000_LONG;
}

//! Function to get the Julian day on zero modified Julian day, in double precision.
template< >
double getJulianDayOnMjd0< double >( )
{
    return JULIAN_DAY_AT_0_MJD;
}

//! Function to get the Julian day on zero modified Julian day, in long double precision.
template< >
long double getJulianDayOnMjd0< long double >( )
{
    return JULIAN_DAY_AT_0_MJD_LONG;
}

//! Function to get the synchronization Julian day of TT, TCG, and TCB, in double precision.
template< >
double getTimeOfTaiSynchronizationJulianDay< double >( )
{
    return TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION;
}

//! Function to get the synchronization Julian day of TT, TCG, and TCB, in long double precision.
template< >
long double getTimeOfTaiSynchronizationJulianDay< long double >( )
{
    return TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION_LONG;
}

//! Function to get the difference between TDB and (TT, TCB and TCB) at TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION,
//! in double precision.
template< >
double getTdbSecondsOffsetAtSynchronization< double >( )
{
    return TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION;
}

//! Function to get the difference between TDB and (TT, TCB and TCB) at TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION,
//! in long double precision.
template< >
long double getTdbSecondsOffsetAtSynchronization< long double >( )
{
    return TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION_LONG;
}

//! Function to get the offset of TT from TAI (constant by definition), in double precision.
template< >
double getTTMinusTai< double >( )
{
    return TT_MINUS_TAI;
}

//! Function to get the offset of TT from TAI (constant by definition), in long double precision.
template< >
long double getTTMinusTai< long double >( )
{
    return TT_MINUS_TAI_LONG;
}

//! Function to get the offset of TT from TAI (constant by definition), in Time format.
template< >
Time getTTMinusTai< Time>( )
{
    return Time( TT_MINUS_TAI_LONG );
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
        throw std::runtime_error( "Error, month number " + std::to_string( month ) +
                                  " does not exist, value must be gretaer than 0 and smaller than 13" );
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
                throw std::runtime_error(
                            "Error when converting year and days in year to date, month number has exceeded 12" );
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
        throw std::runtime_error(
                    "Error when converting year and days in year to date, inconsistent output" );
    }

    return date;
}

//! Perform apprixmate conversion of TT to TDB
double approximateConvertTTtoTDB( const double ttSecondsSinceJ2000 )
{
    double ttCenturiesSinceJ2000 = ttSecondsSinceJ2000 / ( 100.0 * physical_constants::JULIAN_YEAR );
    return ttSecondsSinceJ2000 + 0.001657  * std::sin( 628.3076 * ttCenturiesSinceJ2000 + 6.2401 );
}


} // namespace basic_astrodynamics
} // namespace tudat
