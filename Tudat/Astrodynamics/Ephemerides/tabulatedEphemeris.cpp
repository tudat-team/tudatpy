/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//! Get cartesian state from ephemeris (in double precision), for double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< double, double >::getCartesianState(
        const double ephemerisTime)
{
    return interpolator_->interpolate( ephemerisTime );
}

//! Get cartesian state from ephemeris (in long double precision), for double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< double, double >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    return interpolator_->interpolate( secondsSinceEpoch ).cast< long double >( );
}

//! Get cartesian state from ephemeris (in double precision from Time input), for double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< double, double >::getCartesianStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time.getSeconds< double >( ) );
}

//! Get cartesian state from ephemeris (in long double precision from Time input), for double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< double, double >::getCartesianLongStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time.getSeconds< double >( ) ).cast< long double >( );
}





//! Get cartesian state from ephemeris (in double precision), for long double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< long double, double >::getCartesianState(
        const double ephemerisTime )
{
    return interpolator_->interpolate( ephemerisTime ).cast< double >( );
}

//! Get cartesian state from ephemeris (in long double precision), for long double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< long double, double >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    return interpolator_->interpolate( secondsSinceEpoch );
}

//! Get cartesian state from ephemeris (in double precision from Time input), for double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< long double, double >::getCartesianStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time.getSeconds< double >( ) ).cast< double >( );
}

//! Get cartesian state from ephemeris (in long double precision from Time input), for double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< long double, double >::getCartesianLongStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time.getSeconds< double >( ) );
}


//! Get cartesian state from ephemeris (in double precision), double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< double, Time >::getCartesianState(
        const double ephemerisTime )
{
    return interpolator_->interpolate( Time( ephemerisTime ) );
}

//! Get cartesian state from ephemeris (in long double precision), for double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< double, Time >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    return interpolator_->interpolate( Time( secondsSinceEpoch ) ).cast< long double >( );
}

//! Get cartesian state from ephemeris (in double precision from Time input).
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< double, Time >::getCartesianStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time );
}

//! Get cartesian state from ephemeris (in long double precision from Time input).
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< double, Time >::getCartesianLongStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time ).cast< long double >( );
}



//! Get cartesian state from ephemeris (in double precision), for long double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< long double, Time >::getCartesianState(
        const double ephemerisTime )
{
    return interpolator_->interpolate( Time( ephemerisTime ) ).cast< double >( );
}

//! Get cartesian state from ephemeris (in long double precision), for long double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< long double, Time >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    return interpolator_->interpolate( Time( secondsSinceEpoch ) );
}

//! Get cartesian state from ephemeris (in double precision from Time input).
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< long double, Time >::getCartesianStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time ).cast< double >( );
}

//! Get cartesian state from ephemeris (in long double precision from Time input).
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< long double, Time >::getCartesianLongStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time );
}


//! Function to check whether an ephemeris is a (type of) tabulated ephemeris
bool isTabulatedEphemeris( const std::shared_ptr< Ephemeris > ephemeris )
{
    bool objectIsTabulated = 0;
    if( ( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >( ephemeris ) != nullptr ) ||
            ( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, double > >( ephemeris ) != nullptr ) ||
            ( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, Time > >( ephemeris ) != nullptr ) ||
            ( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >( ephemeris ) != nullptr ) )
    {
        objectIsTabulated = 1;
    }
    return objectIsTabulated;
}

//! Function that retrieves the time interval at which a tabulated ephemeris can be safely interrogated
std::pair< double, double > getTabulatedEphemerisSafeInterval( const std::shared_ptr< Ephemeris > ephemeris )
{
    // Initialize return pair
    std::pair< double, double > safeInterval = std::make_pair( TUDAT_NAN, TUDAT_NAN );

    // Check input consistency
    if( !isTabulatedEphemeris( ephemeris ) )
    {
        throw std::runtime_error( "Error wgen getting tabulated ephemeris safe interval, input is not a tabulated ephemeris" );
    }
    // Identify type of tabulated ephemeris, and call associated safe interval function
    else if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >( ephemeris ) != nullptr )
    {
        safeInterval = std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >(
                    ephemeris )->getSafeInterpolationInterval( );
    }
    else if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, double > >( ephemeris ) != nullptr )
    {
        safeInterval = std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, double > >(
                    ephemeris )->getSafeInterpolationInterval( );
    }
    else if ( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, Time > >( ephemeris ) != nullptr )
    {
        safeInterval = std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, Time > >(
                    ephemeris )->getSafeInterpolationInterval( );
    }
    else if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >( ephemeris ) != nullptr )
     {
        safeInterval = std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >(
                    ephemeris )->getSafeInterpolationInterval( );
    }
    return safeInterval;
}


} // namespace ephemerides

} // namespace tudat

