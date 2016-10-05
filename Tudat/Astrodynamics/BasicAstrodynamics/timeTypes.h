#ifndef TIMETYPES_H
#define TIMETYPES_H

#include <cmath>
#include <algorithm>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

namespace tudat
{

static const long double TIME_NORMALIZATION_TERM = 3600.0L;

class Time
{
public:
    Time( ):fullPeriods_( 0 ), secondsIntoFullPeriod_( 0.0L ){ }

    Time( const int fullPeriods, const long double secondsIntoFullPeriod ):
        fullPeriods_( fullPeriods ), secondsIntoFullPeriod_( secondsIntoFullPeriod )
    {
        normalizeMembers( );
    }

    Time( const long double numberOfSeconds ):
        fullPeriods_( 0 ), secondsIntoFullPeriod_( numberOfSeconds )
    {
        normalizeMembers( );
    }

    Time( const double secondsIntoFullPeriod ):
        fullPeriods_( 0 ), secondsIntoFullPeriod_( static_cast< long double >( secondsIntoFullPeriod ) )
    {
        normalizeMembers( );
    }

    Time( const int secondsIntoFullPeriod ):
        fullPeriods_( 0 ), secondsIntoFullPeriod_( static_cast< long double >( secondsIntoFullPeriod ) )
    {
        normalizeMembers( );
    }

    Time( const Time& otherTime ):
        fullPeriods_( otherTime.fullPeriods_ ), secondsIntoFullPeriod_( otherTime.secondsIntoFullPeriod_ )
    {
        normalizeMembers( );
    }


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
            return *this;
        }
    }



    friend Time operator+( const Time& timeToAdd1, const Time& timeToAdd2 )
    {
        return Time( timeToAdd1.fullPeriods_ + timeToAdd2.fullPeriods_, timeToAdd1.secondsIntoFullPeriod_ + timeToAdd2.secondsIntoFullPeriod_ );
    }

    friend Time operator+( const double& timeToAdd1, const Time& timeToAdd2 )
    {
        return Time( timeToAdd2.fullPeriods_, timeToAdd2.secondsIntoFullPeriod_ + static_cast< long double >( timeToAdd1 ) );
    }

    friend  Time operator+( const long double& timeToAdd1, const Time& timeToAdd2 )
    {
        return Time( timeToAdd2.fullPeriods_, timeToAdd2.secondsIntoFullPeriod_ + timeToAdd1 );
    }

    friend Time operator+( const Time& timeToAdd2, const double& timeToAdd1 )
    {
        return timeToAdd1 + timeToAdd2;
    }

    friend  Time operator+( const Time& timeToAdd2, const long double& timeToAdd1 )
    {
        return timeToAdd1 + timeToAdd2;
    }




    friend Time operator-( const Time& timeToSubtract1, const Time& timeToSubtract2 )
    {

        return Time( timeToSubtract1.fullPeriods_ - timeToSubtract2.fullPeriods_,
                     timeToSubtract1.secondsIntoFullPeriod_ - timeToSubtract2.secondsIntoFullPeriod_ );

    }

    friend Time operator-( const Time& timeToSubtract1, const double timeToSubtract2 )
    {
        return Time( timeToSubtract1.fullPeriods_, timeToSubtract1.secondsIntoFullPeriod_ - static_cast< long double >( timeToSubtract2 ) );
    }

    friend Time operator-( const double timeToSubtract2, const Time& timeToSubtract1 )
    {
        return Time( -timeToSubtract1.fullPeriods_, static_cast< long double >( timeToSubtract2 ) - timeToSubtract1.secondsIntoFullPeriod_ );
    }

    friend Time operator-( const Time& timeToSubtract1, const long double timeToSubtract2 )
    {
        return Time( timeToSubtract1.fullPeriods_, timeToSubtract1.secondsIntoFullPeriod_ - timeToSubtract2 );
    }

    friend Time operator-( const long double timeToSubtract2, const Time& timeToSubtract1 )
    {
        return Time( -timeToSubtract1.fullPeriods_, timeToSubtract2 - timeToSubtract1.secondsIntoFullPeriod_ );
    }

    friend Time operator*( const long double timeToMultiply1, const Time& timeToMultiply2 )
    {
        long double newPeriods = timeToMultiply1 * static_cast< long double >( timeToMultiply2.fullPeriods_ );
        long double roundedNewPeriods = std::floor( newPeriods );

        int newfullPeriods = static_cast< int >( std::round( roundedNewPeriods ) );
        long double newSecondsIntoFullPeriod_ = timeToMultiply2.secondsIntoFullPeriod_ * timeToMultiply1;
        newSecondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods ) * TIME_NORMALIZATION_TERM;

        return Time( newfullPeriods, newSecondsIntoFullPeriod_ );
    }

    friend Time operator*( const double timeToMultiply1, const Time& timeToMultiply2 )
    {
        long double newPeriods = static_cast< long double >( timeToMultiply1 ) * static_cast< long double >( timeToMultiply2.fullPeriods_ );
        long double roundedNewPeriods = std::floor( newPeriods );

        int newfullPeriods = static_cast< int >( std::round( roundedNewPeriods ) );
        long double newSecondsIntoFullPeriod_ = timeToMultiply2.secondsIntoFullPeriod_ * static_cast< long double >( timeToMultiply1 );
        newSecondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods ) * TIME_NORMALIZATION_TERM;

        return Time( newfullPeriods, newSecondsIntoFullPeriod_ );
    }

    void operator+=( const Time& timeToAdd )
    {
        fullPeriods_ += timeToAdd.fullPeriods_;
        secondsIntoFullPeriod_ += timeToAdd.secondsIntoFullPeriod_;
        normalizeMembers( );
    }

    void operator+=( const double timeToAdd )
    {
        secondsIntoFullPeriod_ += static_cast< long double >( timeToAdd );
        normalizeMembers( );
    }

    void operator+=( const long double timeToAdd )
    {
        secondsIntoFullPeriod_ += timeToAdd;
        normalizeMembers( );
    }

    void operator-=( const Time& timeToSubtract )
    {
        fullPeriods_ -= timeToSubtract.fullPeriods_;
        secondsIntoFullPeriod_ -= timeToSubtract.secondsIntoFullPeriod_;
        normalizeMembers( );
    }

    void operator-=( const double timeToSubtract )
    {
        secondsIntoFullPeriod_ -= static_cast< long double >( timeToSubtract );
        normalizeMembers( );
    }

    void operator-=( const long double timeToSubtract )
    {
        secondsIntoFullPeriod_ -= timeToSubtract;
        normalizeMembers( );
    }

    //    void operator*=( const double timeToMultiply )
    //    {
    //        secondsIntoFullPeriod_ *= timeToMultiply;
    //        normalizeMembers( );
    //    }

    void operator*=( const long double timeToMultiply )
    {
        long double newPeriods = timeToMultiply * static_cast< long double >( fullPeriods_ );
        long double roundedNewPeriods = std::floor( newPeriods );

        fullPeriods_ = static_cast< int >( std::round( roundedNewPeriods ) );
        secondsIntoFullPeriod_ *= timeToMultiply;
        secondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods ) * TIME_NORMALIZATION_TERM;

        normalizeMembers( );
    }

//    void operator/=( const long double timeToMultiply )
//    {
//        long double newPeriods = timeToMultiply * static_cast< long double >( fullPeriods_ );
//        long double roundedNewPeriods = std::floor( newPeriods );

//        fullPeriods_ = static_cast< int >( std::round( roundedNewPeriods ) );
//        secondsIntoFullPeriod_ *= timeToMultiply;
//        secondsIntoFullPeriod_ += ( newPeriods - roundedNewPeriods );

//        normalizeMembers( );
//    }


    friend bool operator==( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        return ( ( timeToCompare1.getFullPeriods( ) == timeToCompare2.getFullPeriods( ) ) &&
                 ( timeToCompare1.getSecondsIntoFullPeriod( ) == timeToCompare2.getSecondsIntoFullPeriod( ) ) ) ? true : false ;
    }

    friend bool operator!=( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        return !operator==( timeToCompare1, timeToCompare2 );
    }

    friend bool operator==( const Time& timeToCompare1, const double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< double >( ) == timeToCompare2 ) ? true : false ;
    }

    friend bool operator!=( const Time& timeToCompare1, const double timeToCompare2 )
    {
        return !operator==( timeToCompare1, timeToCompare2 );
    }


    friend bool operator==( const Time& timeToCompare1, const long double timeToCompare2 )
    {
        return ( timeToCompare1.getSeconds< long double >( ) == timeToCompare2 ) ? true : false ;
    }

    friend bool operator!=( const Time& timeToCompare1, const long double timeToCompare2 )
    {
        return !operator==( timeToCompare1, timeToCompare2 );
    }


    friend bool operator> ( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        return timeToCompare1.getSeconds< long double >( ) > timeToCompare2.getSeconds< long double >( );
    }

    friend bool operator>= ( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        return timeToCompare1.getSeconds< long double >( ) >= timeToCompare2.getSeconds< long double >( );
    }

    friend bool operator< ( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        return timeToCompare1.getSeconds< long double >( ) < timeToCompare2.getSeconds< long double >( );
    }

    friend bool operator<= ( const Time& timeToCompare1, const Time& timeToCompare2 )
    {
        return timeToCompare1.getSeconds< long double >( ) <= timeToCompare2.getSeconds< long double >( );
    }


    //    friend const double operator/( const double& original, const Time& timeToDivide )
    //    {
    //        return original / timeToDivide.getSeconds< double >( );
    //    }

    //    friend const long double operator/( const long double& original, const Time& timeToDivide )
    //    {
    //        return original / timeToDivide.getSeconds< long double >( );
    //    }


    friend const Time operator/( const Time& original, const double doubleToDivideBy )
    {
        return Time( static_cast< int >( static_cast< double >( original.fullPeriods_ ) / doubleToDivideBy ),
                     original.secondsIntoFullPeriod_ / static_cast< long double >( doubleToDivideBy ) +
                     basic_mathematics::computeModulo< long double >( static_cast< long double >( original.fullPeriods_ ),
                                                                               static_cast< long double >( doubleToDivideBy ) ) );
    }


    friend const Time operator/( const Time& original, const long double doubleToDivideBy )
    {
        return Time( static_cast< int >( static_cast< long double >( original.fullPeriods_ ) / doubleToDivideBy ),
                     original.secondsIntoFullPeriod_ / doubleToDivideBy + basic_mathematics::computeModulo< long double >(
                         static_cast< long double >( original.fullPeriods_ ), doubleToDivideBy ) );
    }

    friend const Time operator/( const Time& original, const Time& doubleToDivideBy )
    {
        return original / static_cast< long double >( doubleToDivideBy );
    }

    friend std::ostream& operator<<( std::ostream& stream, const Time& timeToPrint )
    {
        stream<<"("<<timeToPrint.getFullPeriods( )<<", "<<timeToPrint.getSecondsIntoFullPeriod( )<<") ";
        return stream;
    }

    template< typename ScalarType >
    ScalarType getSeconds( ) const
    {
        return static_cast< ScalarType >(
                    static_cast< long double >( fullPeriods_ ) * TIME_NORMALIZATION_TERM + secondsIntoFullPeriod_ );
    }

    operator int( ) const { return getSeconds< int >( ); }

    operator double( ) const { return getSeconds< double >( ); }

    operator long double( ) const { return getSeconds< long double >( ); }

    int getFullPeriods( ) const { return fullPeriods_; }

    long double getSecondsIntoFullPeriod( ) const { return secondsIntoFullPeriod_; }

protected:
    void normalizeMembers( )
    {
        if( secondsIntoFullPeriod_ < 0.0L || secondsIntoFullPeriod_ >= TIME_NORMALIZATION_TERM )
        {
            basic_mathematics::computeModuloAndRemainder< long double >(
                        secondsIntoFullPeriod_,TIME_NORMALIZATION_TERM, secondsIntoFullPeriod_, daysToAdd );
            fullPeriods_ += daysToAdd;
        }
    }

    int daysToAdd;

    int fullPeriods_;

    long double secondsIntoFullPeriod_;

};

}

#endif // TUDAT_TIMETYPES_H
