/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//! Get cartesian state from ephemeris (in double precision), for double StateScalarType
template<>
Eigen::Vector6d TabulatedCartesianEphemeris< double, double >::getCartesianState( const double ephemerisTime )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( ephemerisTime );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state using the configured high-precision scalar, for double StateScalarType
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 > TabulatedCartesianEphemeris< double, double >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( secondsSinceEpoch ).cast< HighPrecisionStateScalar >( );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state from ephemeris (in double precision from Time input), for double StateScalarType
template<>
Eigen::Vector6d TabulatedCartesianEphemeris< double, double >::getCartesianStateFromExtendedTime( const Time& time )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( time.getSeconds< double >( ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state using the configured high-precision scalar from Time input, for double StateScalarType
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 > TabulatedCartesianEphemeris< double, double >::getCartesianLongStateFromExtendedTime(
        const Time& time )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( time.getSeconds< double >( ) ).cast< HighPrecisionStateScalar >( );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

#if TUDAT_BUILD_WITH_HIGH_PRECISION_STATE_SCALAR
//! Get cartesian state in double precision, for HighPrecisionStateScalar StateScalarType
template<>
Eigen::Vector6d TabulatedCartesianEphemeris< HighPrecisionStateScalar, double >::getCartesianState( const double ephemerisTime )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( ephemerisTime ).cast< double >( );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state using the configured high-precision scalar
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 > TabulatedCartesianEphemeris< HighPrecisionStateScalar, double >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( secondsSinceEpoch );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state from ephemeris (in double precision from Time input), for double StateScalarType
template<>
Eigen::Vector6d TabulatedCartesianEphemeris< HighPrecisionStateScalar, double >::getCartesianStateFromExtendedTime( const Time& time )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( time.getSeconds< double >( ) ).cast< double >( );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state using the configured high-precision scalar from Time input
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 >
TabulatedCartesianEphemeris< HighPrecisionStateScalar, double >::getCartesianLongStateFromExtendedTime( const Time& time )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( time.getSeconds< double >( ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}
#endif

//! Get cartesian state from ephemeris (in double precision), double StateScalarType
template<>
Eigen::Vector6d TabulatedCartesianEphemeris< double, Time >::getCartesianState( const double ephemerisTime )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( Time( ephemerisTime ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state using the configured high-precision scalar, for double StateScalarType
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 > TabulatedCartesianEphemeris< double, Time >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( Time( secondsSinceEpoch ) ).cast< HighPrecisionStateScalar >( );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state from ephemeris (in double precision from Time input).
template<>
Eigen::Vector6d TabulatedCartesianEphemeris< double, Time >::getCartesianStateFromExtendedTime( const Time& time )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( time );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state using the configured high-precision scalar from Time input.
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 > TabulatedCartesianEphemeris< double, Time >::getCartesianLongStateFromExtendedTime(
        const Time& time )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( time ).cast< HighPrecisionStateScalar >( );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

#if TUDAT_BUILD_WITH_HIGH_PRECISION_STATE_SCALAR
//! Get cartesian state in double precision, for HighPrecisionStateScalar StateScalarType
template<>
Eigen::Vector6d TabulatedCartesianEphemeris< HighPrecisionStateScalar, Time >::getCartesianState( const double ephemerisTime )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( Time( ephemerisTime ) ).cast< double >( );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state using the configured high-precision scalar
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 > TabulatedCartesianEphemeris< HighPrecisionStateScalar, Time >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( Time( secondsSinceEpoch ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state from ephemeris (in double precision from Time input).
template<>
Eigen::Vector6d TabulatedCartesianEphemeris< HighPrecisionStateScalar, Time >::getCartesianStateFromExtendedTime( const Time& time )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( time ).cast< double >( );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}

//! Get cartesian state using the configured high-precision scalar from Time input.
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 >
TabulatedCartesianEphemeris< HighPrecisionStateScalar, Time >::getCartesianLongStateFromExtendedTime( const Time& time )
{
    if( interpolator_ == nullptr )
    {
        throw std::runtime_error( "Error when calling TabulatedCartesianEphemeris, no state interpolator defined" );
    }
    try
    {
        return interpolator_->interpolate( time );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in tabulated ephemeris.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }
}
#endif

//! Function to check whether an ephemeris is a (type of) tabulated ephemeris
bool isTabulatedEphemeris( const std::shared_ptr< Ephemeris > ephemeris )
{
    bool objectIsTabulated = 0;
    if( ( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >( ephemeris ) != nullptr ) ||
#if TUDAT_BUILD_WITH_HIGH_PRECISION_STATE_SCALAR
        ( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< HighPrecisionStateScalar, double > >( ephemeris ) != nullptr ) ||
        ( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< HighPrecisionStateScalar, Time > >( ephemeris ) != nullptr ) ||
#endif
        ( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >( ephemeris ) != nullptr ) )
    {
        objectIsTabulated = 1;
    }
    return objectIsTabulated;
}

//! Function that retrieves the time interval at which a tabulated ephemeris can be safely interrogated
std::pair< double, double > getTabulatedEphemerisSafeInterval( const std::shared_ptr< Ephemeris > ephemeris,
                                                               const bool acceptUserDefinedRisk )
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
        safeInterval = std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >( ephemeris )
                               ->getSafeInterpolationInterval( acceptUserDefinedRisk );
    }
#if TUDAT_BUILD_WITH_HIGH_PRECISION_STATE_SCALAR
    else if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< HighPrecisionStateScalar, double > >( ephemeris ) != nullptr )
    {
        safeInterval = std::dynamic_pointer_cast< TabulatedCartesianEphemeris< HighPrecisionStateScalar, double > >( ephemeris )
                               ->getSafeInterpolationInterval( acceptUserDefinedRisk );
    }
    else if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< HighPrecisionStateScalar, Time > >( ephemeris ) != nullptr )
    {
        safeInterval = std::dynamic_pointer_cast< TabulatedCartesianEphemeris< HighPrecisionStateScalar, Time > >( ephemeris )
                               ->getSafeInterpolationInterval( acceptUserDefinedRisk );
    }
#endif
    else if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >( ephemeris ) != nullptr )
    {
        safeInterval = std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >( ephemeris )
                               ->getSafeInterpolationInterval( acceptUserDefinedRisk );
    }
    else
    {
        throw std::runtime_error( "Error when getting tabulated ephemeris safe interval, tabulated ephemeris not recognized" );
    }
    return safeInterval;
}

}  // namespace ephemerides

}  // namespace tudat
