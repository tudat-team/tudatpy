/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DATETIME_H
#define TUDAT_DATETIME_H

#include <chrono>
#include <ctime>

#include "tudat/basics/timeType.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/deprecationWarnings.h"

namespace tudat
{

namespace basic_astrodynamics
{

template< typename TimeType >
TimeType timeFromDecomposedDateTime( const int year,
                                     const int month,
                                     const int day,
                                     const int hour,
                                     const int minutes,
                                     const long double seconds )
{
    int daysSinceJ2000 = julianDayNumberFromDate( year, month, day ) - basic_astrodynamics::JULIAN_DAY_ON_J2000_INT;

    if( TIME_NORMALIZATION_INTEGER_TERM != 3600 )
    {
        throw std::runtime_error( "Error when getting time from decomposed date time; normalization period has been changed!" );
    }
    long double secondsIntoCurrentHour = static_cast< long double >( minutes ) * 60.0L + seconds;

    int fullPeriods = daysSinceJ2000 * TIME_NORMALIZATION_TERMS_PER_DAY + ( hour - 12 );
    long double secondsIntoFullPeriod = secondsIntoCurrentHour;
    return static_cast< TimeType >( Time( fullPeriods, secondsIntoFullPeriod ) );
}

//! Function to get Time from an ISO time string
inline void decomposedDateTimeFromIsoString( const std::string &isoTime,
                                             int &year,
                                             int &month,
                                             int &days,
                                             int &hours,
                                             int &minutes,
                                             long double &seconds )
{
    try
    {
        // Get year and month
        std::vector< std::string > splitTime;
        boost::algorithm::split( splitTime, isoTime, boost::is_any_of( "-‐" ), boost::algorithm::token_compress_on );
        year = boost::lexical_cast< int >( splitTime.at( 0 ) );
        month = boost::lexical_cast< int >( splitTime.at( 1 ) );

        // Get day
        std::string remainingString = splitTime.at( 2 );
        splitTime.clear( );
        boost::algorithm::split( splitTime, remainingString, boost::is_any_of( "T " ), boost::algorithm::token_compress_on );

        days = boost::lexical_cast< int >( splitTime.at( 0 ) );

        // Get hours, minutes, seconds
        remainingString = splitTime.at( 1 );
        splitTime.clear( );
        boost::algorithm::split( splitTime, remainingString, boost::is_any_of( ":" ), boost::algorithm::token_compress_on );
        hours = boost::lexical_cast< int >( splitTime.at( 0 ) );
        minutes = boost::lexical_cast< int >( splitTime.at( 1 ) );
        seconds = boost::lexical_cast< long double >( splitTime.at( 2 ) );
    }
    catch( std::runtime_error &caughtException )
    {
        throw std::runtime_error( "Error when parsing iso datetime string " + isoTime +
                                  ". Caught exception is: " + caughtException.what( ) );
    }
}

//! Function to get Time from an ISO time string
template< typename TimeType >
TimeType timeFromIsoString( const std::string &isoTime )
{
    int year, month, days, hours, minutes;
    long double seconds;

    decomposedDateTimeFromIsoString( isoTime, year, month, days, hours, minutes, seconds );
    return timeFromDecomposedDateTime< TimeType >( year, month, days, hours, minutes, seconds );
}

struct DateTime {
public:
    DateTime( int year, int month, int day, int hour, int minute, long double seconds ):
        year_( year ), month_( month ), day_( day ), hour_( hour ), minute_( minute ), seconds_( seconds )
    {
        verifyMonth( );
        verifyDay( );
        verifyHour( );
        verifyMinute( );
        verifySeconds( );
    }

    int getYear( ) const
    {
        return year_;
    }

    int getMonth( ) const
    {
        return month_;
    }

    int getDay( ) const
    {
        return day_;
    }

    int getHour( ) const
    {
        return hour_;
    }

    int getMinute( ) const
    {
        return minute_;
    }

    long double getSeconds( ) const
    {
        return seconds_;
    }

    void setYear( const int year )
    {
        year_ = year;
        verifyDay( );
    }

    void setMonth( const int month )
    {
        month_ = month;
        verifyMonth( );
        verifyDay( );
    }

    void setDay( const int day )
    {
        day_ = day;
        verifyDay( );
    }

    void setHour( const int hour )
    {
        hour_ = hour;
        verifyHour( );
    }

    void setMinute( const int minute )
    {
        minute_ = minute;
        verifyMinute( );
    }

    void setSeconds( const long double seconds )
    {
        seconds_ = seconds;
        verifySeconds( );
    }

    std::string isoString( const bool addT = false, const int numberOfFractionalSecondDigits = 15 ) const
    {
        std::string yearString = std::to_string( year_ );
        std::string monthString = utilities::paddedZeroIntString( month_, 2 );
        std::string dayString = utilities::paddedZeroIntString( day_, 2 );
        std::string hourString = utilities::paddedZeroIntString( hour_, 2 );
        std::string minuteString = utilities::paddedZeroIntString( minute_, 2 );
        std::string secondString = utilities::paddedZeroIntString( static_cast< int >( seconds_ ), 2 );
        long double fractionalSeconds =
                seconds_ - mathematical_constants::getFloatingInteger< long double >( static_cast< int >( seconds_ ) );
        std::string fractionalSecondString =
                utilities::to_string_with_precision< long double >( fractionalSeconds, numberOfFractionalSecondDigits );

        secondString += fractionalSecondString.substr( 1, fractionalSecondString.length( ) - 1 );
        std::string separationCharacter = addT ? "T" : " ";
        return yearString + "-" + monthString + "-" + dayString + separationCharacter + hourString + ":" + minuteString + ":" +
                secondString;
    }

    template< typename TimeType >
    TimeType julianDay( )
    {
        return julianDayFromTime< TimeType >( epoch< Time >( ) );
    }

    template< typename TimeType >
    TimeType modifiedJulianDay( )
    {
        return modifiedJulianDayFromTime< TimeType >( epoch< Time >( ) );
    }

    template< typename TimeType >
    TimeType epoch( ) const
    {
        return timeFromDecomposedDateTime< TimeType >( year_, month_, day_, hour_, minute_, seconds_ );
    }

    int dayOfYear( )
    {
        return basic_astrodynamics::convertDayMonthYearToDayOfYear( day_, month_, year_ );
    }

    static double minimumChronoRepresentableEpoch( )
    {
        // std::chrono::system_clock uses the Unix epoch (1970-01-01 00:00:00 UTC) as reference point
        DateTime referenceDateTime( 1970, 1, 1, 0, 0, 0.0L );

        constexpr std::chrono::system_clock::duration::rep minTickCount =
                std::numeric_limits< std::chrono::system_clock::duration::rep >::min( );

        std::chrono::system_clock::duration minDuration( minTickCount );

        std::chrono::duration< double > secondsDuration = minDuration;

        double lowerRepresentationCount = secondsDuration.count( );
        // const double lowerRepresentationCount =
        //         std::chrono::duration< double >( std::chrono::system_clock::time_point::min( ).time_since_epoch( ) ).count( );

        double minimumRepresentableEpoch = referenceDateTime.epoch< double >( ) + lowerRepresentationCount;

        std::cout << "Lower representation count" << lowerRepresentationCount << std::endl;
        std::cout << "Minimum representable epoch" << minimumRepresentableEpoch << std::endl;

        return minimumRepresentableEpoch;
    }

    static double maximumChronoRepresentableEpoch( )
    {
        // std::chrono::system_clock uses the Unix epoch (1970-01-01 00:00:00 UTC) as reference point
        DateTime referenceDateTime( 1970, 1, 1, 0, 0, 0.0L );

        constexpr std::chrono::system_clock::duration::rep maxTickCount =
                std::numeric_limits< std::chrono::system_clock::duration::rep >::max( );

        std::chrono::system_clock::duration maxDuration( maxTickCount );

        std::chrono::duration< double > secondsDuration = maxDuration;

        double upperRepresentationCount = secondsDuration.count( );
        // const double upperRepresentationCount =
        //         std::chrono::duration< double >( std::chrono::system_clock::time_point::max( ).time_since_epoch( ) ).count( );

        double maximumRepresentableEpoch = referenceDateTime.epoch< double >( ) + upperRepresentationCount;

        std::cout << "Upper representation count" << upperRepresentationCount << std::endl;
        std::cout << "Maximum representable epoch" << maximumRepresentableEpoch << std::endl;

        return maximumRepresentableEpoch;
    }

    std::chrono::system_clock::time_point timePoint( ) const
    {
        double minimumChronoEpoch = minimumChronoRepresentableEpoch( );
        double maximumChronoEpoch = maximumChronoRepresentableEpoch( );

        if( this->epoch< double >( ) > maximumChronoEpoch || this->epoch< double >( ) < minimumChronoEpoch )
        {
            throw std::runtime_error( " Date " + this->isoString( false, 3 ) +
                                      " is out of range for conversion to time point. Lower limit (in seconds from J2000) is: " +
                                      std::to_string( minimumChronoEpoch ) + ", upper limit: " + std::to_string( maximumChronoEpoch ) );
        }

        std::tm tm = { };
        tm.tm_sec = static_cast< int >( this->getSeconds( ) );
        tm.tm_min = this->getMinute( );
        tm.tm_hour = this->getHour( );
        tm.tm_mday = this->getDay( );
        tm.tm_mon = this->getMonth( ) - 1;
        tm.tm_year = this->getYear( ) - 1900;

        tm.tm_isdst = -1;

        std::cout << "Tudat IsoString: " << this->isoString( false, 3 ) << std::endl;

        std::cout << "tm from tudat: " << local_tm.tm_year << " " << local_tm.tm_mon << " " << local_tm.tm_mday << " " << local_tm.tm_hour
                  << " " << local_tm.tm_min << " " << local_tm.tm_sec << std::endl;

        std::time_t tt = std::mktime( &tm );
        std::cout << "time_t from mktime: " << std::to_string( tt ) << std::endl;

        std::chrono::system_clock::time_point timePoint = std::chrono::system_clock::from_time_t( tt );
        return timePoint +
                std::chrono::microseconds(
                        static_cast< int >( std::round( ( this->getSeconds( ) - static_cast< long double >( tm.tm_sec ) ) *
                                                        tudat::mathematical_constants::getFloatingInteger< long double >( 1E6 ) ) ) );
    }

    static DateTime fromTimePoint( const std::chrono::system_clock::time_point datetime )
    {
        std::time_t tt = std::chrono::system_clock::to_time_t( datetime );

        std::cout << "time_t from time_point: " << std::to_string( tt ) << std::endl;

        std::tm local_tm = *localtime( &tt );

        using namespace std::chrono;
        microseconds timeInMicroSeconds = duration_cast< microseconds >( datetime.time_since_epoch( ) );
        long long fractional_seconds = timeInMicroSeconds.count( ) % 1000000LL;

        std::cout << "Local tm: " << local_tm.tm_year << " " << local_tm.tm_mon << " " << local_tm.tm_mday << " " << local_tm.tm_hour << " "
                  << local_tm.tm_min << " " << local_tm.tm_sec << std::endl;
        std::cout << "Time in microseconds: " << timeInMicroSeconds.count( ) << std::endl;
        std::cout << "Fractional seconds: " << fractional_seconds << std::endl;

        std::cout << "Input seconds to DateTime: "
                  << static_cast< long double >( local_tm.tm_sec ) +
                        static_cast< long double >( fractional_seconds ) /
                                tudat::mathematical_constants::getFloatingInteger< long double >( 1000000LL )
                  << std::endl;

        return DateTime( local_tm.tm_year + 1900,
                         local_tm.tm_mon + 1,
                         local_tm.tm_mday,
                         local_tm.tm_hour,
                         local_tm.tm_min,
                         static_cast< long double >( local_tm.tm_sec ) +
                                 static_cast< long double >( fractional_seconds ) /
                                         tudat::mathematical_constants::getFloatingInteger< long double >( 1000000LL ) );
    }
    static DateTime fromYearAndDaysInYear( const int year, const int daysInYear )
    {
        boost::gregorian::date boostDateTime = convertYearAndDaysInYearToDate( year, daysInYear );
        return DateTime( boostDateTime.year( ), boostDateTime.month( ), boostDateTime.day( ), 0, 0, 0.0 );
    }

    template< typename TimeType >
    static DateTime fromTime( const TimeType &timeInput )
    {
        std::cout << "Creating DateTime from time " << timeInput << std::endl;
        Time time = Time( timeInput );
        std::cout << "Creating DateTime from time " << time.getFullPeriods( ) << std::endl;

        int fullPeriodsSinceMidnightJD0 = time.getFullPeriods( ) +
                basic_astrodynamics::JULIAN_DAY_ON_J2000_INT * TIME_NORMALIZATION_TERMS_PER_DAY + TIME_NORMALIZATION_TERMS_PER_HALF_DAY;
        if( fullPeriodsSinceMidnightJD0 < 0 )
        {
            throw std::runtime_error( "Error calendar date from Time not implemented for negative Julian days" );
        }
        int fullDaysSinceMidnightJD0 = fullPeriodsSinceMidnightJD0 / TIME_NORMALIZATION_TERMS_PER_DAY;
        int day, month, year;

        basic_astrodynamics::convertShiftedJulianDayToCalendarDate< int >( fullDaysSinceMidnightJD0, day, month, year );

        if( TIME_NORMALIZATION_INTEGER_TERM != 3600 )
        {
            throw std::runtime_error( "Error, Time to calendar date only implemented for normalization term of 3600 s" );
        }
        int hour = time.fullPeriodsSinceMidnight( );
        int minute = std::floor( time.getSecondsIntoFullPeriod( ) / 60.0L );

        std::cout << "Test output A" << minute << " " << time.getSecondsIntoFullPeriod( ) << " " << 60.0L << std::endl;
        long double seconds = time.getSecondsIntoFullPeriod( ) - static_cast< long double >( 60 * minute );
        std::cout << "Test output B" << seconds << " " << time.getSecondsIntoFullPeriod( ) << " "
                  << static_cast< long double >( 60 * minute ) << " " << 60 * minute << std::endl;

        return DateTime( year, month, day, hour, minute, seconds );
    }

    static DateTime fromIsoString( const std::string &isoTime )
    {
        int year, month, days, hours, minutes;
        long double seconds;

        decomposedDateTimeFromIsoString( isoTime, year, month, days, hours, minutes, seconds );
        return DateTime( year, month, days, hours, minutes, seconds );
    }

    template< typename TimeType >
    DateTime addSecondsToDateTime( const TimeType timeToAdd )
    {
        return fromTime< Time >( this->epoch< Time >( ) + timeToAdd );
    }

    template< typename TimeType >
    DateTime addDaysToDateTime( const TimeType daysToAdd )
    {
        return fromTime< Time >( this->epoch< Time >( ) + daysToAdd * mathematical_constants::getFloatingInteger< long double >( 86400 ) );
    }

protected:
    int year_;
    int month_;
    int day_;
    int hour_;
    int minute_;
    long double seconds_;

    void verifyMonth( )
    {
        if( month_ > 12 || month_ < 1 )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input month was " + std::to_string( month_ ) );
        }
    }

    void verifyDay( )
    {
        std::cout << "Verifying day " << day_ << " " << month_ << " " << year_ << " "
                  << basic_astrodynamics::getDaysInMonth( month_, year_ ) << std::endl;
        if( day_ > basic_astrodynamics::getDaysInMonth( month_, year_ ) )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input date was " + std::to_string( day_ ) + "-" +
                                      std::to_string( month_ ) + "-" + std::to_string( year_ ) );
        }
    }

    void verifyHour( )
    {
        if( hour_ > 23 || hour_ < 0 )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input hour was " + std::to_string( hour_ ) );
        }
    }

    void verifyMinute( )
    {
        if( minute_ > 59 || minute_ < 0 )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input minute was " + std::to_string( minute_ ) );
        }
    }

    void verifySeconds( )
    {
        if( seconds_ > 60.0L || seconds_ < 0.0L || ( seconds_ != seconds_ ) )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input seconds was " + std::to_string( seconds_ ) +
                                      ", full date time was " + std::to_string( year_ ) + ", " + std::to_string( month_ ) + ", " +
                                      std::to_string( day_ ) + ", " + std::to_string( hour_ ) + ", " + std::to_string( minute_ ) + ", " +
                                      std::to_string( seconds_ ) );
        }
    }
};

template< typename TimeType >
DateTime addSecondsToDateTime( const DateTime &dateTime, const TimeType timeToAdd )
{
    utilities::printDeprecationWarning( "add_seconds_to_datetime", "DateTime.add_seconds" );

    return DateTime::fromTime< Time >( dateTime.epoch< Time >( ) + timeToAdd );
}

template< typename TimeType >
DateTime addDaysToDateTime( const DateTime &dateTime, const TimeType daysToAdd )
{
    utilities::printDeprecationWarning( "add_days_to_datetime", "DateTime.add_days" );
    return DateTime::fromTime< Time >( dateTime.epoch< Time >( ) +
                                       daysToAdd * mathematical_constants::getFloatingInteger< long double >( 86400 ) );
}

template< typename TimeType >
TimeType getTimeDifferenceBetweenDateTimes( const DateTime &firstDateTime, const DateTime &secondDateTime )
{
    return firstDateTime.epoch< TimeType >( ) - secondDateTime.epoch< TimeType >( );
}

}  // namespace basic_astrodynamics

}  // namespace tudat
#endif  // TUDAT_DATETIME_H
