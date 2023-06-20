/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TIMETYPE_H
#define TUDAT_TIMETYPE_H

#include <cmath>
#include <algorithm>

#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/timeConversions.h"

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"

namespace tudat
{

static constexpr int TIME_NORMALIZATION_INTEGER_TERM = 3600;

static constexpr long double TIME_NORMALIZATION_TERM = mathematical_constants::getFloatingInteger< long double >( TIME_NORMALIZATION_INTEGER_TERM );

static constexpr int TIME_NORMALIZATION_TERMS_PER_DAY = physical_constants::JULIAN_DAY_INT / TIME_NORMALIZATION_INTEGER_TERM;

static constexpr int TIME_NORMALIZATION_TERMS_PER_HALF_DAY = TIME_NORMALIZATION_TERMS_PER_DAY / 2;

static constexpr int J2000_JULIAN_DAY_IN_FULL_PERIODS = physical_constants::JULIAN_DAY_INT * TIME_NORMALIZATION_TERMS_PER_DAY;

struct DateTime
{
    DateTime( int year, int month, int day, int hour, int minute, long double seconds ):
        year_( year ), month_( month ), day_ (day ), hour_( hour ), minute_( minute ), seconds_( seconds )
    {
        if( month > 12 || month < 1 )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input month was " + std::to_string( month ) );
        }

        if(  day > basic_astrodynamics::getDaysInMonth( month, year ) )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input date was " +
                std::to_string( day ) + "-" + std::to_string( month ) + "-" + std::to_string( year ) );
        }

        if( hour > 23 || hour < 0 )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input hour was " + std::to_string( hour ) );
        }

        if( minute > 59 || minute < 0 )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input minute was " + std::to_string( minute ) );
        }

        if( seconds > 60.0L || seconds < 0.0L || (seconds != seconds ) )
        {
            throw std::runtime_error( "Error when creating Tudat DateTime, input seconds was " + std::to_string( seconds ) );
        }
    }

    int year_;
    int month_;
    int day_;
    int hour_;
    int minute_;
    long double seconds_;
};

//! Class for defining time with a resolution that is sub-fs for very long periods of time.
/*!
 *  Class for defining time with a resolution that is sub-fs for very long periods of time. Using double or long double
 *  precision as a representation of time, the issue of reduced quality will occur that over long time-period. For instance,
 *  over a period of 10^8 seconds (about 3 years), double and long double representations have resolution of about 10^-8 and
 *  10^-11 s respectively, which is insufficient for various applications. This type uses an int to represent the number of
 *  hours since an epoch, and long double to represent the number of seconds into the present hour. This provides a
 *  resulution of < 1 femtosecond, over a range of 2147483647 hours (about 300,000 years), which is more than sufficient for
 *  practical applications.
 */
class Time
{
public:

    constexpr Time( ):fullPeriods_( 0 ), secondsIntoFullPeriod_( 0.0L ), daysToAdd( 0 ){ }

    //! Constructor, sets current hour and time into current hour directly
    /*!
     * Constructor, sets current hour and time into current hour directly
     * \param fullPeriods Number of full hours since epoch
     * \param secondsIntoFullPeriod Number of seconds into current hour. Note that this value need not be in the range
     * between 0 and 3600: the time representation is normalized upon construction to ensure that the internal representation
     * is in this range.
     */
    constexpr Time( const int fullPeriods, const long double secondsIntoFullPeriod ):
        fullPeriods_( fullPeriods ), secondsIntoFullPeriod_( secondsIntoFullPeriod ), daysToAdd( 0 )
    {
        normalizeMembers( );
    }

    //! Constructor, sets number of seconds since epoch (with long double representation as input)
    /*!
     * Constructor, sets number of seconds since epoch (with long double representation as input)
     * \param numberOfSeconds Number of seconds since epoch.
     */
    constexpr Time( const long double numberOfSeconds ):
        fullPeriods_( 0 ), secondsIntoFullPeriod_( numberOfSeconds ), daysToAdd( 0 )
    {
        normalizeMembers( );
    }

    //! Constructor, sets number of seconds since epoch (with double representation as input)
    /*!
     * Constructor, sets number of seconds since epoch (with double representation as input)
     * \param secondsIntoFullPeriod Number of seconds since epoch.
     */
    constexpr Time( const double secondsIntoFullPeriod ):
        fullPeriods_( 0 ), secondsIntoFullPeriod_( static_cast< long double >( secondsIntoFullPeriod ) ), daysToAdd( 0 )
    {
        normalizeMembers( );
    }

    //! Constructor, sets number of seconds since epoch (with int representation as input)
    /*!
     * Constructor, sets number of seconds since epoch (with int representation as input)
     * \param secondsIntoFullPeriod Number of seconds since epoch.
     */
    constexpr Time( const int secondsIntoFullPeriod ):
        fullPeriods_( 0 ), secondsIntoFullPeriod_( static_cast< long double >( secondsIntoFullPeriod ) ), daysToAdd( 0 )
    {
        normalizeMembers( );
    }

    //! Copy constructor
    /*!
     * Copy constructor
     * \param otherTime Time that is to be copied.
     */
    constexpr Time( const Time& otherTime ):
        fullPeriods_( otherTime.fullPeriods_ ), secondsIntoFullPeriod_( otherTime.secondsIntoFullPeriod_ ), daysToAdd( 0 )
    {
        normalizeMembers( );
    }

    //! Definition of = operator for Time type
    /*!
     * Definition of = operator for Time type
     * \param timeToCopy Time that is to be copied by operator
     * \return Assigned Time object
     */
    Time& operator=( const Time& timeToCopy )
    {
        if( this == &timeToCopy )
        {
            return *this;
        }
        else
        {
            secondsIntoFullPeriod_ = timeToCopy.secondsIntoFullPeriod_;
            fullPeriods_ = timeToCopy.fullPeriods_;
            normalizeMembers( );
            return *this;
        }
    }


    //! Addition operator for two Time objects
    /*!
     * Addition operator for two Time objects
     * \param timeToAdd1 First time that is to be added.
     * \param timeToAdd2 Second time that is to be added.
     * \return Input arguments, added together
     */
    friend Time operator+( const Time& timeToAdd1, const Time& timeToAdd2 )
    {
        return Time( timeToAdd1.fullPeriods_ + timeToAdd2.fullPeriods_,
                     timeToAdd1.secondsIntoFullPeriod_ + timeToAdd2.secondsIntoFullPeriod_ );
    }

    //! Addition operator for double variable with Time object.
    /*!
     * Addition operator for two Time objects
     * \param timeToAdd1 First time that is to be added (as a double).
     * \param timeToAdd2 Second time that is to be added (as a Time object).
     * \return Input arguments, added together as Time object.
     */
    friend Time operator+( const double& timeToAdd1, const Time& timeToAdd2 )
    {
        return Time( timeToAdd2.fullPeriods_, timeToAdd2.secondsIntoFullPeriod_ + static_cast< long double >( timeToAdd1 ) );
    }

    //! Addition operator for long double variable with Time object.
    /*!
     * Addition operator for two Time objects
     * \param timeToAdd1 First time that is to be added (as a long double).
     * \param timeToAdd2 Second time that is to be added (as a Time object).
     * \return Input arguments, added together as Time object.
     */
    friend Time operator+( const long double& timeToAdd1, const Time& timeToAdd2 )
    {
        return Time( timeToAdd2.fullPeriods_, timeToAdd2.secondsIntoFullPeriod_ + timeToAdd1 );
    }

    //! Addition operator for Time object with double variable
    /*!
     * Addition operator for two Time objects
     * \param timeToAdd1 First time that is to be added (as Time object).
     * \param timeToAdd2 Second time that is to be added (as a double).
     * \return Input arguments, added together as Time object.
     */
    friend Time operator+( const Time& timeToAdd2, const double& timeToAdd1 )
    {
        return timeToAdd1 + timeToAdd2;
    }

    //! Addition operator for Time object with long double variable
    /*!
     * Addition operator for two Time objects
     * \param timeToAdd1 First time that is to be added (as Time object).
     * \param timeToAdd2 Second time that is to be added (as a long double).
     * \return Input arguments, added together as Time object.
     */
    friend Time operator+( const Time& timeToAdd2, const long double& timeToAdd1 )
    {
        return timeToAdd1 + timeToAdd2;
    }



    //! Subtraction operator for two Time objects
    /*!
     * Addition operator for two Time objects
     * \param timeToSubtract1 Time from which second time is to be subtracted
     * \param timeToSubtract2 Time that is to be subtracted from first input
     * \return Input arguments, subtracted from one another
     */
    friend Time operator-( const Time& timeToSubtract1, const Time& timeToSubtract2 )
    {

        return Time( timeToSubtract1.fullPeriods_ - timeToSubtract2.fullPeriods_,
                     timeToSubtract1.secondsIntoFullPeriod_ - timeToSubtract2.secondsIntoFullPeriod_ );
    }

    //! Subtraction operator for double from Time object
    /*!
     * Addition operator for two Time objects
     * \param timeToSubtract1 Time from which second time is to be subtracted (as a Time object)
     * \param timeToSubtract2 Time that is to be subtracted from first input (as a double)
     * \return Input arguments, subtracted from one another
     */
    friend Time operator-( const Time& timeToSubtract1, const double timeToSubtract2 )
    {
        return Time( timeToSubtract1.fullPeriods_,
                     timeToSubtract1.secondsIntoFullPeriod_ - static_cast< long double >( timeToSubtract2 ) );
    }

    //! Subtraction operator for double from Time object
    /*!
     * Addition operator for two Time objects
     * \param timeToSubtract1 Time from which second time is to be subtracted (as a Time object)
     * \param timeToSubtract2 Time that is to be subtracted from first input (as a long double)
     * \return Input arguments, subtracted from one another
     */
    friend Time operator-( const Time& timeToSubtract1, const long double timeToSubtract2 )
    {
        return Time( timeToSubtract1.fullPeriods_, timeToSubtract1.secondsIntoFullPeriod_ - timeToSubtract2 );
    }

    //! Subtraction operator for Time object from double
    /*!
     * Addition operator for two Time objects
     * \param timeToSubtract1 Time from which second time is to be subtracted (as a double)
     * \param timeToSubtract2 Time that is to be subtracted from first input (as a Time object)
     * \return Input arguments, subtracted from one another
     */
    friend Time operator-( const double timeToSubtract1, const Time& timeToSubtract2 )
    {
        return Time( -timeToSubtract2.fullPeriods_,
                     static_cast< long double >( timeToSubtract1 ) - timeToSubtract2.secondsIntoFullPeriod_ );
    }

    //! Subtraction operator for Time object from long double
    /*!
     * Addition operator for two Time objects
     * \param timeToSubtract1 Time from which second time is to be subtracted (as a long double)
     * \param timeToSubtract2 Time that is to be subtracted from first input (as a Time object)
     * \return Input arguments, subtracted from one another
     */
    friend Time operator-( const long double timeToSubtract1, const Time& timeToSubtract2 )
    {
        return Time( -timeToSubtract2.fullPeriods_, timeToSubtract1 - timeToSubtract2.secondsIntoFullPeriod_ );
    }


    //! Multiplication operator of a long double with a Time object (i.e. to rescale time)
    /*!
     * Multiplication operator of a long double with a Time object (i.e. to rescale time)
     * \param timeToMultiply1 Value by which Time is to be multiplied
     * \param timeToMultiply2 Time that is to be multiplied by first input argument
     * \return Multiplied Time object.
     */
    friend Time operator*( const long double timeToMultiply1, const Time& timeToMultiply2 )
    {
        long double newPeriods = timeToMultiply1 * static_cast< long double >( timeToMultiply2.fullPeriods_ );
        long double roundedNewPeriods = static_cast< long double >( std::floor( newPeriods ) );

        int newfullPeriods = static_cast< int >( std::round( roundedNewPeriods ) );
        long double newSecondsIntoFullPeriod_ = timeToMultiply2.secondsIntoFullPeriod_ * timeToMultiply1;
        newSecondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods ) * TIME_NORMALIZATION_TERM;

        return Time( newfullPeriods, newSecondsIntoFullPeriod_ );
    }

    //! Multiplication operator of a long double with a Time object (i.e. to rescale time)
    /*!
     * Multiplication operator of a long double with a Time object (i.e. to rescale time)
     * \param timeToMultiply1 Time that is to be multiplied by second input argument
     * \param timeToMultiply2 Value by which Time is to be multiplied
     * \return Multiplied Time object.
     */
    friend Time operator*( const Time& timeToMultiply1, const long double timeToMultiply2 )
    {
        return timeToMultiply2 * timeToMultiply1;
    }

    //! Multiplication operator of a double with a Time object (i.e. to rescale time)
    /*!
     * Multiplication operator of a double with a Time object (i.e. to rescale time)
     * \param timeToMultiply1 Value by which Time is to be multiplied
     * \param timeToMultiply2 Time that is to be multiplied by first input argument
     * \return Multiplied Time object.
     */
    friend Time operator*( const double timeToMultiply1, const Time& timeToMultiply2 )
    {
        long double newPeriods = static_cast< long double >( timeToMultiply1 ) *
                static_cast< long double >( timeToMultiply2.fullPeriods_ );
        long double roundedNewPeriods = std::floor( newPeriods );

        int newfullPeriods = static_cast< int >( std::round( roundedNewPeriods ) );
        long double newSecondsIntoFullPeriod_ = timeToMultiply2.secondsIntoFullPeriod_ *
                static_cast< long double >( timeToMultiply1 );
        newSecondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods ) * TIME_NORMALIZATION_TERM;

        return Time( newfullPeriods, newSecondsIntoFullPeriod_ );
    }

    //! Multiplication operator of a double with a Time object (i.e. to rescale time)
    /*!
     * Multiplication operator of a double with a Time object (i.e. to rescale time)
     * \param timeToMultiply1 Time that is to be multiplied by second input argument
     * \param timeToMultiply2 Value by which Time is to be multiplied
     * \return Multiplied Time object.
     */
    friend Time operator*( const Time& timeToMultiply1, const double timeToMultiply2 )
    {
        return timeToMultiply2 * timeToMultiply1;
    }

    //! Division operator of a Time object with a double (i.e. to rescale time)
    /*!
     * Division operator of a Time object with a double (i.e. to rescale time)
     * \param original Time that is to be divided by second input argument
     * \param doubleToDivideBy Value by which first argument is to be divided.
     * \return Divided Time object.
     */
    friend const Time operator/( const Time& original, const double doubleToDivideBy )
    {
        long double newPeriods =  static_cast< long double >( original.fullPeriods_ ) /
                static_cast< long double >( doubleToDivideBy );

        long double roundedNewPeriods = std::floor( newPeriods );

        int newfullPeriods = static_cast< int >( std::round( roundedNewPeriods ) );
        long double newSecondsIntoFullPeriod_ = original.secondsIntoFullPeriod_ /
                static_cast< long double >( doubleToDivideBy );
        newSecondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods ) * TIME_NORMALIZATION_TERM;

        return Time( newfullPeriods, newSecondsIntoFullPeriod_ );
    }


    //! Division operator of a Time object with a long double (i.e. to rescale time)
    /*!
     * Division operator of a Time object with a long double (i.e. to rescale time)
     * \param original Time that is to be divided by second input argument
     * \param doubleToDivideBy Value by which first argument is to be divided.
     * \return Divided Time object.
     */
    friend const Time operator/( const Time& original, const long double doubleToDivideBy )
    {
        long double newPeriods =  static_cast< long double >( original.fullPeriods_ ) / doubleToDivideBy;

        long double roundedNewPeriods = std::floor( newPeriods );

        int newfullPeriods = static_cast< int >( std::round( roundedNewPeriods ) );
        long double newSecondsIntoFullPeriod_ = original.secondsIntoFullPeriod_ /  doubleToDivideBy;
        newSecondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods ) * TIME_NORMALIZATION_TERM;

        return Time( newfullPeriods, newSecondsIntoFullPeriod_ );
    }


    //! Add and assign operator for adding a Time
    /*!
     *  Add and assign operator for adding a Time
     *  \param timeToAdd Time that is to be added
     */
    void operator+=( const Time& timeToAdd )
    {
        fullPeriods_ += timeToAdd.fullPeriods_;
        secondsIntoFullPeriod_ += timeToAdd.secondsIntoFullPeriod_;
        normalizeMembers( );
    }

    //! Add and assign operator for adding a double
    /*!
     *  Add and assign operator for adding a double
     *  \param timeToAdd Time that is to be added (in double representation)
     */
    void operator+=( const double timeToAdd )
    {
        secondsIntoFullPeriod_ += static_cast< long double >( timeToAdd );
        normalizeMembers( );
    }

    //! Add and assign operator for adding a double
    /*!
     *  Add and assign operator for adding a double
     *  \param timeToAdd Time that is to be added (in long double representation)
     */
    void operator+=( const long double timeToAdd )
    {
        secondsIntoFullPeriod_ += timeToAdd;
        normalizeMembers( );
    }

    //! Subtract and assign operator for adding a Time
    /*!
     *  Subtract and assign operator for subtracting a Time
     *  \param timeToSubtract Time that is to be subtracted
     */
    void operator-=( const Time& timeToSubtract )
    {
        fullPeriods_ -= timeToSubtract.fullPeriods_;
        secondsIntoFullPeriod_ -= timeToSubtract.secondsIntoFullPeriod_;
        normalizeMembers( );
    }

    //! Subtract and assign operator for adding a double
    /*!
     *  Subtract and assign operator for subtracting a double
     *  \param timeToSubtract Time that is to be subtracted (in double representation)
     */
    void operator-=( const double timeToSubtract )
    {
        secondsIntoFullPeriod_ -= static_cast< long double >( timeToSubtract );
        normalizeMembers( );
    }

    //! Subtract and assign operator for adding a long double
    /*!
     *  Subtract and assign operator for subtracting a long double
     *  \param timeToSubtract Time that is to be subtracted (in long double representation)
     */
    void operator-=( const long double timeToSubtract )
    {
        secondsIntoFullPeriod_ -= timeToSubtract;
        normalizeMembers( );
    }

    //! Multiply and assign operator for multiplying by double
    /*!
     *  Multiply and assign operator for multiplying by a double
     *  \param timeToMultiply Value by which time is to be multiplied (in double representation)
     */
    void operator*=( const double timeToMultiply )
    {
        long double newPeriods = static_cast< long double >( timeToMultiply ) * static_cast< long double >( fullPeriods_ );
        long double roundedNewPeriods = std::floor( newPeriods );

        fullPeriods_ = static_cast< int >( std::round( roundedNewPeriods ) );
        secondsIntoFullPeriod_ *= static_cast< long double >( timeToMultiply );
        secondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods ) * TIME_NORMALIZATION_TERM;

        normalizeMembers( );
    }

    //! Multiply and assign operator for multiplying by long double
    /*!
     *  Multiply and assign operator for multiplying by a long double
     *  \param timeToMultiply Value by which time is to be multiplied (in long double representation)
     */
    void operator*=( const long double timeToMultiply )
    {
        long double newPeriods = timeToMultiply * static_cast< long double >( fullPeriods_ );
        long double roundedNewPeriods = std::floor( newPeriods );

        fullPeriods_ = static_cast< int >( std::round( roundedNewPeriods ) );
        secondsIntoFullPeriod_ *= timeToMultiply;
        secondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods ) * TIME_NORMALIZATION_TERM;

        normalizeMembers( );
    }

    //! Divided and assign operator for dividing by double
    /*!
     *  Divided and assign operator for dividing by double
     *  \param timeToDivide Value by which time is to be divided (in double representation)
     */
    void operator/=( const double timeToDivide )
    {
        long double newPeriods = static_cast< long double >( timeToDivide ) / static_cast< long double >( fullPeriods_ );
        long double roundedNewPeriods = std::floor( newPeriods );

        fullPeriods_ = static_cast< int >( std::round( roundedNewPeriods ) );
        secondsIntoFullPeriod_ *= timeToDivide;
        secondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods );

        normalizeMembers( );
    }

    //! Divided and assign operator for dividing by long double
    /*!
     *  Divided and assign operator for dividing by long double
     *  \param timeToDivide Value by which time is to be divided (in long  double representation)
     */
    void operator/=( const long double timeToDivide )
    {
        long double newPeriods = timeToDivide / static_cast< long double >( fullPeriods_ );
        long double roundedNewPeriods = std::floor( newPeriods );

        fullPeriods_ = static_cast< int >( std::round( roundedNewPeriods ) );
        secondsIntoFullPeriod_ *= timeToDivide;
        secondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods );

        normalizeMembers( );
    }


    //! Equality operator for two Time objects
    /*!
     * Equality operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare
     * \return True if two times are fully equal; false if not.
     */
    friend bool operator==( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        return ( ( timeToCompare1.getFullPeriods( ) == timeToCompare2.getFullPeriods( ) ) &&
                 ( timeToCompare1.getSecondsIntoFullPeriod( ) == timeToCompare2.getSecondsIntoFullPeriod( ) ) );
    }

    //! Inequality operator for two Time objects
    /*!
     * Inqquality operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare
     * \return False if two times are fully equal; true if not.
     */
    friend bool operator!=( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        return !operator==( timeToCompare1, timeToCompare2 );
    }

    //! Equality operator for a Time object with an integer.
    /*!
     *  Equality operator for a Time object with an integer. Comparison is performed at double precision (i.e. integer is
     *  cast to double and compared).
     *  \param timeToCompare1 First time to compare
     *  \param timeToCompare2 Second time to compare (in integer precision)
     *  \return True if two times are fully equal; false if not.
     */
    friend bool operator==( const Time& timeToCompare1, const int timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< double >( ) == static_cast< double >( timeToCompare2 ) );
    }

    //! Equality operator for a Time object with a double.
    /*!
     *  Equality operator for a Time object with a double. Comparison is performed at double precision (i.e. Time is
     *  cast to double and compared)
     *  \param timeToCompare1 First time to compare
     *  \param timeToCompare2 Second time to compare (in double precision)
     *  \return True if two times are fully equal; false if not.
     */
    friend bool operator==( const Time& timeToCompare1, const double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< double >( ) == timeToCompare2 );
    }

    //! Equality operator for a Time object with a double
    /*!
     *  Equality operator for a Time object with a double. Comparison is performed at double precision (i.e. Time is
     *  cast to double and compared)
     *  \param timeToCompare1 First time to compare (in double precision)
     *  \param timeToCompare2 Second time to compare
     *  \return True if two times are fully equal; false if not.
     */
    friend bool operator==( const double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare2.getSeconds< double >( ) == timeToCompare1 );
    }

    //! Inequality operator for a Time object with a double
    /*!
     *  Inequality operator for a Time object with a double. Comparison is performed at double precision (i.e. Time is
     *  cast to double and compared)
     *  \param timeToCompare1 First time to compare
     *  \param timeToCompare2 Second time to compare (in double precision)
     *  \return True if two times are fully equal; false if not.
     */
    friend bool operator!=( const Time& timeToCompare1, const double timeToCompare2 )
    {
        return !operator==( timeToCompare1, timeToCompare2 );
    }

    //! Inequality operator for a Time object with a double
    /*!
     *  Inequality operator for a Time object with a double. Comparison is performed at double precision (i.e. Time is
     *  cast to double and compared)
     *  \param timeToCompare1 First time to compare (in double precision)
     *  \param timeToCompare2 Second time to compare
     *  \return True if two times are fully equal; false if not.
     */
    friend bool operator!=( const double timeToCompare1, const Time& timeToCompare2 )
    {
        return !operator==( timeToCompare1, timeToCompare2 );
    }

    //! Equality operator for a Time object with a long double
    /*!
     *  Equality operator for a Time object with a long double. Comparison is performed at lgon double precision (i.e. Time
     *  is cast to long double and compared)
     *  \param timeToCompare1 First time to compare
     *  \param timeToCompare2 Second time to compare (in long double precision)
     *  \return True if two times are fully equal; false if not.
     */
    friend bool operator==( const Time& timeToCompare1, const long double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< long double >( ) == timeToCompare2 );
    }

    //! Equality operator for a Time object with a long double
    /*!
     *  Equality operator for a Time object with a long double. Comparison is performed at lgon double precision (i.e. Time
     *  is cast to long double and compared)
     *  \param timeToCompare1 First time to compare (in long double precision)
     *  \param timeToCompare2 Second time to compare
     *  \return True if two times are fully equal; false if not.
     */
    friend bool operator==( const long double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare2.getSeconds< long double >( ) == timeToCompare1 );
    }

    //! Inequality operator for a Time object with a long double
    /*!
     *  Inequality operator for a Time object with a long double. Comparison is performed at lgon double precision (i.e. Time
     *  is cast to long double and compared)
     *  \param timeToCompare1 First time to compare
     *  \param timeToCompare2 Second time to compare (in long double precision)
     *  \return True if two times are fully equal; false if not.
     */
    friend bool operator!=( const Time& timeToCompare1, const long double timeToCompare2 )
    {
        return !operator==( timeToCompare1, timeToCompare2 );
    }

    //! Inequality operator for a Time object with a long double
    /*!
     *  Inequality operator for a Time object with a long double. Comparison is performed at lgon double precision (i.e. Time
     *  is cast to long double and compared)
     *  \param timeToCompare1 First time to compare (in long double precision)
     *  \param timeToCompare2 Second time to compare
     *  \return True if two times are fully equal; false if not.
     */
    friend bool operator!=( const long double timeToCompare1, const Time& timeToCompare2 )
    {
        return !operator==( timeToCompare1, timeToCompare2 );
    }

    //! Greater-than operator for two Time objects
    /*!
     * Greater-than operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is larger than timeToCompare2, false otherwise.
     */
    friend bool operator> ( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        if( timeToCompare1.getFullPeriods( ) > timeToCompare2.getFullPeriods( ) )
        {
            return true;
        }
        else if( ( timeToCompare1.getFullPeriods( ) == timeToCompare2.getFullPeriods( ) ) &&
                 ( timeToCompare1.getSecondsIntoFullPeriod( ) > timeToCompare2.getSecondsIntoFullPeriod( ) ) )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    //! Greater-than-or-equal-to operator for two Time objects
    /*!
     * Greater-than-or-equal-to-than operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is larger than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator>= ( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        if( timeToCompare1.getFullPeriods( ) >= timeToCompare2.getFullPeriods( ) )
        {
            return true;
        }
        else if( ( timeToCompare1.getFullPeriods( ) == timeToCompare2.getFullPeriods( ) ) &&
          ( timeToCompare1.getSecondsIntoFullPeriod( ) >= timeToCompare2.getSecondsIntoFullPeriod( ) ) )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    //! Smaller-than operator for two Time objects
    /*!
     * Smaller-than operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is smaller than timeToCompare2, false otherwise.
     */
    friend bool operator< ( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        if( timeToCompare1.getFullPeriods( ) < timeToCompare2.getFullPeriods( ) )
        {
            return true;
        }
        else if( ( timeToCompare1.getFullPeriods( ) == timeToCompare2.getFullPeriods( ) ) &&
                 ( timeToCompare1.getSecondsIntoFullPeriod( ) < timeToCompare2.getSecondsIntoFullPeriod( ) ) )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    //! Smaller-than-or-equal-to operator for two Time objects
    /*!
     * Smaller-than-or-equal-to-than operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is smaller than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator<= ( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        if( timeToCompare1.getFullPeriods( ) < timeToCompare2.getFullPeriods( ) )
        {
            return true;
        }
        else if( ( timeToCompare1.getFullPeriods( ) == timeToCompare2.getFullPeriods( ) ) &&
                 ( timeToCompare1.getSecondsIntoFullPeriod( ) <= timeToCompare2.getSecondsIntoFullPeriod( ) ) )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    //! Smaller-than operator for Time object with double
    /*!
     * Smaller-than operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare (as double)
     * \return True if timeToCompare1 is smaller than timeToCompare2, false otherwise.
     */
    friend bool operator< ( const Time& timeToCompare1, const double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< double >( ) < timeToCompare2 );
    }

    //! Smaller-than operator for Time object with long double
    /*!
     * Smaller-than operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare (as double)
     * \return True if timeToCompare1 is smaller than timeToCompare2, false otherwise.
     */
    friend bool operator< ( const Time& timeToCompare1, const long double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< long double >( ) < timeToCompare2 );
    }

    //! Smaller-than-or-equal operator for Time object with double
    /*!
     * Smaller-than-or-equal operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare (as double)
     * \return True if timeToCompare1 is smaller than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator<= ( const Time& timeToCompare1, const double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< double >( ) <= timeToCompare2 );
    }

    //! Smaller-than-or-equal operator for Time object with long double
    /*!
     * Smaller-than-or-equal operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare (as double)
     * \return True if timeToCompare1 is smaller than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator<= ( const Time& timeToCompare1, const long double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< long double >( ) <= timeToCompare2 );
    }

    //! Greater-than operator for Time object with double
    /*!
     * Greater-than operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare (as double)
     * \return True if timeToCompare1 is larger than timeToCompare2, false otherwise.
     */
    friend bool operator> ( const Time& timeToCompare1, const double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< double >( ) > timeToCompare2 );
    }

    //! Greater-than operator for Time object with long double
    /*!
     * Greater-than operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare (as double)
     * \return True if timeToCompare1 is larger than timeToCompare2, false otherwise.
     */
    friend bool operator> ( const Time& timeToCompare1, const long double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< long double >( ) > timeToCompare2 );
    }

    //! Greater-than-or-equal operator for Time object with double
    /*!
     * Greater-than-or-equal operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare (as double)
     * \return True if timeToCompare1 is larger than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator>= ( const Time& timeToCompare1, const double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< double >( ) >= timeToCompare2 );
    }

    //! Greater-than-or-equal operator for Time object with long double
    /*!
     * Greater-than-or-equal operator for two Time objects
     * \param timeToCompare1 First time to compare
     * \param timeToCompare2 Second time to compare (as double)
     * \return True if timeToCompare1 is larger than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator>= ( const Time& timeToCompare1, const long double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< long double >( ) >= timeToCompare2 );
    }

    //! Smaller-than operator for Time object with double
    /*!
     * Smaller-than operator for two Time objects
     * \param timeToCompare1 First time to compare (as double)
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is smaller than timeToCompare2, false otherwise.
     */
    friend bool operator< ( const double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare1 < timeToCompare2.getSeconds< double >( ) );
    }

    //! Smaller-than operator for Time object with long double
    /*!
     * Smaller-than operator for two Time objects
     * \param timeToCompare1 First time to compare (as long double)
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is smaller than timeToCompare2, false otherwise.
     */
    friend bool operator< ( const long double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare1 < timeToCompare2.getSeconds< long double >( ) );
    }

    //! Smaller-than-or-equal operator for Time object with double
    /*!
     * Smaller-than-or-equal operator for two Time objects
     * \param timeToCompare1 First time to compare (as double)
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is smaller than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator<= ( const double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare1 <= timeToCompare2.getSeconds< double >( ) );
    }

    //! Smaller-than-or-equal operator for Time object with long double
    /*!
     * Smaller-than-or-equal operator for two Time objects
     * \param timeToCompare1 First time to compare (as long double)
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is smaller than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator<= ( const long double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare1 <= timeToCompare2.getSeconds< long double >( ) );
    }

    //! Greater-than operator for Time object with double
    /*!
     * Greater-than operator for two Time objects
     * \param timeToCompare1 First time to compare (as double)
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is larger than timeToCompare2, false otherwise.
     */
    friend bool operator> ( const double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare1 > timeToCompare2.getSeconds< double >( ) );
    }

    //! Greater-than operator for Time object with long double
    /*!
     * Greater-than operator for two Time objects
     * \param timeToCompare1 First time to compare (as long double)
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is larger than timeToCompare2, false otherwise.
     */
    friend bool operator> ( const long double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare1 > timeToCompare2.getSeconds< long double >( ) );
    }

    //! Greater-than-or-equal operator for Time object with double
    /*!
     * Greater-than-or-equal operator for two Time objects
     * \param timeToCompare1 First time to compare (as double)
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is larger than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator>= ( const double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare1 >= timeToCompare2.getSeconds< double >( ) );
    }

    //! Greater-than-or-equal operator for Time object with long double
    /*!
     * Greater-than-or-equal operator for two Time objects
     * \param timeToCompare1 First time to compare (as long double)
     * \param timeToCompare2 Second time to compare
     * \return True if timeToCompare1 is larger than or equal to timeToCompare2, false otherwise.
     */
    friend bool operator>= ( const long double timeToCompare1, const Time& timeToCompare2 )
    {
        return ( timeToCompare1 >= timeToCompare2.getSeconds< long double >( ) );
    }

    //!Output operator for Time object
    friend std::ostream& operator << ( std::ostream& stream, const Time& timeToPrint )
    {
        stream << "(" << timeToPrint.getFullPeriods( ) << ", " << timeToPrint.getSecondsIntoFullPeriod( ) << ") ";
        return stream;
    }

    //! Function to get the total seconds since epoch, in templated precision
    /*!
     *  Function to get the total seconds since epoch, in templated precision
     *  \return Total seconds since epoch.
     */
    template< typename ScalarType >
    ScalarType getSeconds( ) const
    {
        return static_cast< ScalarType >(
                    static_cast< long double >( fullPeriods_ ) * TIME_NORMALIZATION_TERM + secondsIntoFullPeriod_ );
    }

    //! Function to get the total seconds since epoch, in int precision (cast of Time to int)
    /*!
     *  Function to get the total seconds since epoch, in int precision (cast of Time to int)
     *  \return Total seconds int epoch.
     */
    operator int( ) const
    {
        return getSeconds< int >( );
    }

    //! Function to get the total seconds since epoch, in int precision (cast of Time to double)
    /*!
     *  Function to get the total seconds since epoch, in int precision (cast of Time to double)
     *  \return Total seconds int epoch.
     */
    operator double( ) const
    {
        return getSeconds< double >( );
    }

    //! Function to get the total seconds since epoch, in int precision (cast of Time to long sdouble)
    /*!
     *  Function to get the total seconds since epoch, in int precision (cast of Time to long double)
     *  \return Total seconds int epoch.
     */
    operator long double( ) const
    {
        return getSeconds< long double >( );
    }

    //! Function to get the number of full hours since epoch
    /*!
     * \brief Function to get the number of full hours since epoch
     * \return Number of full hours since epoch
     */
    int getFullPeriods( ) const
    {
        return fullPeriods_;
    }

    //! Function to get the number of seconds into current hour
    /*!
     * \brief Function to get the number of seconds into current hour
     * \return
     */
    long double getSecondsIntoFullPeriod( ) const
    {
        return secondsIntoFullPeriod_;
    }

    int fullDaysSinceEpoch( ) const
    {
        if( fullPeriods_ >= 0 || ( ( fullPeriods_ % TIME_NORMALIZATION_TERMS_PER_DAY ) == 0 ) )
        {
            return fullPeriods_ / TIME_NORMALIZATION_TERMS_PER_DAY;
        }
        else
        {
            return fullPeriods_ / TIME_NORMALIZATION_TERMS_PER_DAY - 1;
        }
    }

    int fullPeriodsIntoCurrentDay( ) const
    {
        return fullPeriods_ - fullDaysSinceEpoch( ) * TIME_NORMALIZATION_TERMS_PER_DAY;
    }

    long double secondsSinceNoon( ) const
    {
        return static_cast< long double >( fullPeriodsIntoCurrentDay( ) * TIME_NORMALIZATION_TERM ) + secondsIntoFullPeriod_;
    }

    int fullPeriodsSinceMidnight( ) const
    {
        int fullPeriodsIntoCurrentDay = this->fullPeriodsIntoCurrentDay( );
        if( fullPeriodsIntoCurrentDay < TIME_NORMALIZATION_TERMS_PER_HALF_DAY)
        {
            return fullPeriodsIntoCurrentDay + TIME_NORMALIZATION_TERMS_PER_HALF_DAY;
        }
        else
        {
            return fullPeriodsIntoCurrentDay - TIME_NORMALIZATION_TERMS_PER_HALF_DAY;
        }
    }

    long double secondsSinceMidnight( ) const
    {
        return static_cast< long double >( fullPeriodsSinceMidnight( ) * TIME_NORMALIZATION_TERM ) + secondsIntoFullPeriod_;
    }

protected:

    //! Function to renormalize the members of the Time object, so that secondsIntoFullPeriod_ is between 0 and 3600
    constexpr void normalizeMembers( )
    {
        if( secondsIntoFullPeriod_ < 0.0L || secondsIntoFullPeriod_ >= TIME_NORMALIZATION_TERM )
        {
            basic_mathematics::computeModuloAndRemainder< long double >(
                        secondsIntoFullPeriod_,TIME_NORMALIZATION_TERM, secondsIntoFullPeriod_, daysToAdd );
            fullPeriods_ += daysToAdd;
        }
    }

    //! Number of full hours since epoch
    int fullPeriods_;

    //! Number of seconds into current hour
    long double secondsIntoFullPeriod_;

    //! Pre-declared variable used in often-called normalizeMembers function
    int daysToAdd;

};

template< typename TimeType >
TimeType timeFromDecomposedDateTime(
    const int year, const int month, const int day,
    const int hour, const int minutes, const long double seconds )
{
    boost::gregorian::date currentDate( year, month, day );
    int daysSinceJ2000 = currentDate.julian_day( ) - basic_astrodynamics::JULIAN_DAY_ON_J2000_INT;

    if( TIME_NORMALIZATION_INTEGER_TERM != 3600 )
    {
        throw std::runtime_error( "Error when getting time from decomposed date time; normalization period has been changed!" );
    }
    long double secondsIntoCurrentHour =
        static_cast< long double >( minutes ) * 60.0L +
        seconds;

    int fullPeriods = daysSinceJ2000 * TIME_NORMALIZATION_TERMS_PER_DAY + ( hour - 12 );
    long double secondsIntoFullPeriod = secondsIntoCurrentHour;
    return static_cast< TimeType >( Time( fullPeriods, secondsIntoFullPeriod ) );
}


//! The Time at JD0
constexpr static Time TIME_AT_JD0 =  Time(
    -basic_astrodynamics::JULIAN_DAY_ON_J2000_INT *  tudat::TIME_NORMALIZATION_TERMS_PER_DAY,
    0.0L );

//! The Time at MJD0
constexpr static Time TIME_AT_MJD0 = Time(
    -51545 * tudat::TIME_NORMALIZATION_TERMS_PER_DAY,
    static_cast< long double >( TIME_NORMALIZATION_TERMS_PER_HALF_DAY ) * TIME_NORMALIZATION_TERM );

//! Function to get the Julian day from the current Time
template< typename ScalarType >
ScalarType julianDayFromTime( const Time& time )
{
    return ( time - TIME_AT_JD0 ).fullDaysSinceEpoch( ) + time.secondsSinceNoon( ) / physical_constants::getJulianDay< ScalarType >( );
}

//! Function to get the modified Julian day from the current Time
template< typename ScalarType >
ScalarType modifiedJulianDayFromTime( const Time& time )
{
    Time scaledTime = time - TIME_AT_MJD0;
    return scaledTime.fullDaysSinceEpoch( ) + scaledTime.secondsSinceNoon( ) / physical_constants::getJulianDay< ScalarType >( );
}

//! Function to get Time from the current Julian day
template< typename ScalarType >
Time timeFromJulianDay( const ScalarType julianDay )
{
    return Time( 0.0, ( julianDay - basic_astrodynamics::getJulianDayOnJ2000< ScalarType >( ) ) * physical_constants::getJulianDay< ScalarType >( ) );
}

//! Function to get Time from the current modified Julian day
template< typename ScalarType >
Time timeFromModifiedJulianDay( const ScalarType julianDay )
{
    return Time( 0.0, ( julianDay - basic_astrodynamics::getModifiedJulianDayOnJ2000< ScalarType >( ) ) * physical_constants::getJulianDay< ScalarType >( ) );
}

//! Function to get Time from an ISO time string
template< typename TimeType >
TimeType timeFromIsoString( const std::string& isoTime )
{

    try
    {
        // Get year and month
        std::vector<std::string> splitTime;
        boost::algorithm::split( splitTime, isoTime,
                                 boost::is_any_of( "-‐" ),
                                 boost::algorithm::token_compress_on );
        int year = boost::lexical_cast<int>( splitTime.at( 0 ));
        int month = boost::lexical_cast<int>( splitTime.at( 1 ));

        // Get day
        std::string remainingString = splitTime.at( 2 );
        splitTime.clear( );
        boost::algorithm::split( splitTime, remainingString,
                                 boost::is_any_of( "T " ),
                                 boost::algorithm::token_compress_on );

        int days = boost::lexical_cast<int>( splitTime.at( 0 ));

        // Get hours, minutes, seconds
        remainingString = splitTime.at( 1 );
        splitTime.clear( );
        boost::algorithm::split( splitTime, remainingString,
                                 boost::is_any_of( ":" ),
                                 boost::algorithm::token_compress_on );
        int hours = boost::lexical_cast<int>( splitTime.at( 0 ));
        int minutes = boost::lexical_cast<int>( splitTime.at( 1 ));
        long double seconds = boost::lexical_cast<long double>( splitTime.at( 2 ));

        return timeFromDecomposedDateTime<TimeType>(
            year, month, days, hours, minutes, seconds );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error when parsing iso datetime string " + isoTime + ". Caught exception is: " +
                                      caughtException.what( ) );
    }
}


//! Function to get Time from an ISO time string
template< typename TimeType >
TimeType timeFromDateTime( const DateTime& dateTime )
{
    return timeFromDecomposedDateTime< TimeType >(
        dateTime.year_, dateTime.month_, dateTime.day_, dateTime.hour_, dateTime.minute_, dateTime.seconds_ );
}

template< typename TimeType >
inline boost::gregorian::date convertTimeToCalendarDate( const TimeType timeInput )
{
    Time time = Time( timeInput );
    long double julianDay = julianDayFromTime< long double >( time );
    boost::gregorian::date currentDate = basic_astrodynamics::convertJulianDayToCalendarDate( julianDay );
    return currentDate;
}

template< typename TimeType >
inline DateTime getCalendarDateFromTime( const TimeType& timeInput )
{
    Time time = Time( timeInput );
    int fullPeriodsSinceMidnightJD0 = time.getFullPeriods( ) + basic_astrodynamics::JULIAN_DAY_ON_J2000_INT * TIME_NORMALIZATION_TERMS_PER_DAY +  TIME_NORMALIZATION_TERMS_PER_HALF_DAY;
    if( fullPeriodsSinceMidnightJD0 < 0 )
    {
        throw std::runtime_error( "Error calendar date from Time not implemented for negative Julian days" );
    }
    int fullDaysSinceMidnightJD0 = fullPeriodsSinceMidnightJD0 / TIME_NORMALIZATION_TERMS_PER_DAY;
    int day, month, year;
    basic_astrodynamics::convertShiftedJulianDayToCalendarDate< int >(
        fullDaysSinceMidnightJD0, day, month, year);

    if( TIME_NORMALIZATION_INTEGER_TERM != 3600 )
    {
        throw std::runtime_error( "Error, Time to calendar date only implemented for normalization term of 3600 s" );
    }
    int hour = time.fullPeriodsSinceMidnight( );
    int minute = std::floor( time.getSecondsIntoFullPeriod( ) / 60.0L );
    long double seconds = time.getSecondsIntoFullPeriod( ) - 60.0L * static_cast< long double >( minute );

    return DateTime( year, month, day, hour, minute, seconds );
}


} // namespace tudat

#endif // TUDAT_TIMETYPE_H
