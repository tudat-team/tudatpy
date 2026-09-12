/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cmath>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"

namespace tudat
{

namespace sofa_interface
{

//! Function to calculate number of leap seconds from UTC input
double getDeltaAtFromUtc( const double utcInJulianDays )
{
    // Declare return variables (by reference) from sofa calendar function.
    int year, month, day;
    double fractionOfDay;

    // Get calendar date and check feasibility of input.
    if( iauJd2cal( basic_astrodynamics::JULIAN_DAY_ON_J2000, utcInJulianDays, &year, &month, &day, &fractionOfDay ) != 0 )
    {
        throw std::runtime_error( "Provided julian date too small to convert to calendar date" + std::to_string( year ) + " " +
                                  std::to_string( month ) + " " + std::to_string( day ) );
    }

    // Get number of leap seconds and check feasibility of calculation
    double deltaAt;
    int deltaAtReturn = iauDat( year, month, day, fractionOfDay, &deltaAt );
    if( deltaAtReturn != 0 )
    {
        throw std::runtime_error( "Provided caledar date cannot properly give Delta AT" + std::to_string( year ) + " " +
                                  std::to_string( month ) + " " + std::to_string( day ) );
    }

    return deltaAt;
}

//! Function to calculate number of leap seconds from TAI input
double getDeltaAtFromTai( const double taiInJulianDays )
{
    // Declare return variables (by reference) from sofa calendar function.
    int year, month, day;
    double fractionOfDay;

    // Get calendar date and check feasibility of input.
    if( iauJd2cal( basic_astrodynamics::JULIAN_DAY_ON_J2000, taiInJulianDays, &year, &month, &day, &fractionOfDay ) != 0 )
    {
        throw std::runtime_error( "Provided julian date too small to convert to calendar date" );
    }

    // Estimate number of leap seconds (by assuming TAI = UTC for Sofa input) and check feasibility of calculation
    double deltaAt;
    int deltaAtReturn = iauDat( year, month, day, fractionOfDay, &deltaAt );
    if( deltaAtReturn != 0 )
    {
        throw std::runtime_error( "Provided caledar date cannot properly give Delta AT" );
    }

    // Reperform calculation with converted utc time and check consistency with previous calculation.
    double utcInJulianDays = taiInJulianDays - deltaAt / physical_constants::JULIAN_DAY;
    double deltaAtCheck = getDeltaAtFromUtc( utcInJulianDays );
    if( deltaAt != deltaAtCheck )
    {
        deltaAt--;
        throw std::runtime_error( "Warning, Delta TAI calculation encountered error, iteration not yet implemented" );
    }

    return deltaAt;
}

// double convertUTCtoUT1( const double utc )
//{
//	return 0.0; // iauUtcut1()
// }

//! Function to calculate difference between TDB and TT
double getTDBminusTT( const double tdbTime,
                      const double universalTimeFractionOfDay,
                      const double stationLongitude,
                      const double distanceFromSpinAxis,
                      const double distanceFromEquatorialPlane )
{
    return iauDtdb( basic_astrodynamics::JULIAN_DAY_ON_J2000,
                    tdbTime / physical_constants::JULIAN_DAY,
                    universalTimeFractionOfDay,
                    stationLongitude,
                    distanceFromSpinAxis / 1000.0,
                    distanceFromEquatorialPlane / 1000.0 );
}

//! Function to calculate difference between TDB and TT.
double getTDBminusTT( const double tdbTime, const double universalTimeFractionOfDay, const Eigen::Vector3d& stationCartesianPosition )
{
    double stationLongitude = std::atan2( stationCartesianPosition( 1 ), stationCartesianPosition( 0 ) );
    return getTDBminusTT( tdbTime,
                          universalTimeFractionOfDay,
                          stationLongitude,
                          std::sqrt( stationCartesianPosition.x( ) * stationCartesianPosition.x( ) +
                                     stationCartesianPosition.y( ) * stationCartesianPosition.y( ) ),
                          stationCartesianPosition.z( ) );
}

//! Function to calculate difference between TDB and TT.
double getTDBminusTT( const double ttOrTdbSinceJ2000,
                      const double stationLongitude,
                      const double distanceFromSpinAxis,
                      const double distanceFromEquatorialPlane )
{
    // Calculate current TAI (approximately if input is in TDB)
    double tai = basic_astrodynamics::convertTTtoTAI< double >( ttOrTdbSinceJ2000 );

    // Calculate current UT1 (by assuming it equal to UTC)
    double ut1 = TUDAT_NAN;

    // Conversion is only valid from 1961 onwards
    if( static_cast< double >( tai ) / physical_constants::JULIAN_DAY >
        ( basic_astrodynamics::JULIAN_DAY_OF_UTC_INTRODUCTION - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) )
    {
        ut1 = static_cast< double >( convertTAItoUTC< double >( tai ) );
    }
    else
    {
        // Rough approximation, but suficient for TT<->TDB computation
        ut1 = tai;
    }
    double ut1FractionOfDay = std::fmod(
            ( ut1 / physical_constants::JULIAN_DAY ) - static_cast< double >( std::floor( ut1 / physical_constants::JULIAN_DAY ) ) + 0.5,
            1.0 );

    // Calculate and return difference (introducing addition approximation if input is in TT, by assuming TDB is equal to TT)
    double tdbMinusTT =
            getTDBminusTT( ttOrTdbSinceJ2000, ut1FractionOfDay, stationLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane );

    return tdbMinusTT;
}

//! Function to calculate difference between TDB and TT.
double getTDBminusTT( const double ttOrTdbSinceJ2000, const Eigen::Vector3d& earthFixedPosition )
{
    double siteLongitude = std::atan2( earthFixedPosition.y( ), earthFixedPosition.x( ) );
    double distanceFromSpinAxis =
            std::sqrt( earthFixedPosition.x( ) * earthFixedPosition.x( ) + earthFixedPosition.y( ) * earthFixedPosition.y( ) );
    double distanceFromEquatorialPlane = earthFixedPosition.z( );

    return getTDBminusTT( ttOrTdbSinceJ2000, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane );
}

#if TUDAT_HIGH_PRECISION_STATE_SCALAR_IS_CPP_BIN_FLOAT_QUAD
// Evaluation adapted from ERFA dtdb.c; the coefficient include carries the upstream license.
HighPrecisionStateScalar getHighPrecisionTDBminusTT( const HighPrecisionStateScalar& seconds, const Eigen::Vector3d& position )
{
    using Scalar = HighPrecisionStateScalar;
    using std::atan2;
    using std::cos;
    using std::floor;
    using std::sin;
    using std::sqrt;

#include "fairheadBretagnonCoefficients.inc"

    const Scalar t = seconds / Scalar( 31557600000LL );
    const Scalar tai = seconds - Scalar( 32184 ) / Scalar( 1000 );
    Scalar utc = tai;
    if( tai / Scalar( 86400 ) > Scalar( basic_astrodynamics::JULIAN_DAY_OF_UTC_INTRODUCTION - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) )
    {
        // Leap offsets are exact integers in the supported modern UTC era. Only calendar lookup uses double.
        utc = convertTAItoUTC< Scalar >( tai );
    }
    const Scalar pi = mathematical_constants::getPi< Scalar >( );
    const Scalar ut = utc / Scalar( 86400 ) + Scalar( 1 ) / Scalar( 2 );
    const Scalar tsol = ( ut - floor( ut ) ) * Scalar( 2 ) * pi + atan2( Scalar( position.y( ) ), Scalar( position.x( ) ) );
    const Scalar u =
            sqrt( Scalar( position.x( ) ) * Scalar( position.x( ) ) + Scalar( position.y( ) ) * Scalar( position.y( ) ) ) / Scalar( 1000 );
    const Scalar v = Scalar( position.z( ) ) / Scalar( 1000 );
    const Scalar w = t / Scalar( 3600 );
    const Scalar degreesToRadians = pi / Scalar( 180 );
    const Scalar elsun = ( Scalar( 280.46645683 ) + Scalar( 1296027711.03429 ) * w ) * degreesToRadians;
    const Scalar emsun = ( Scalar( 357.52910918 ) + Scalar( 1295965810.481 ) * w ) * degreesToRadians;
    const Scalar d = ( Scalar( 297.85019547 ) + Scalar( 16029616012.090 ) * w ) * degreesToRadians;
    const Scalar elj = ( Scalar( 34.35151874 ) + Scalar( 109306899.89453 ) * w ) * degreesToRadians;
    const Scalar els = ( Scalar( 50.07744430 ) + Scalar( 44046398.47038 ) * w ) * degreesToRadians;
    const Scalar topocentric = Scalar( 0.00029e-10 ) * u * sin( tsol + elsun - els ) +
            Scalar( 0.00100e-10 ) * u * sin( tsol - Scalar( 2 ) * emsun ) + Scalar( 0.00133e-10 ) * u * sin( tsol - d ) +
            Scalar( 0.00133e-10 ) * u * sin( tsol + elsun - elj ) - Scalar( 0.00229e-10 ) * u * sin( tsol + Scalar( 2 ) * elsun + emsun ) -
            Scalar( 0.02200e-10 ) * v * cos( elsun + emsun ) + Scalar( 0.05312e-10 ) * u * sin( tsol - emsun ) -
            Scalar( 0.13677e-10 ) * u * sin( tsol + Scalar( 2 ) * elsun ) - Scalar( 1.31840e-10 ) * v * cos( elsun ) +
            Scalar( 3.17679e-10 ) * u * sin( tsol );

    const int edges[] = { 0, 474, 679, 764, 784, 787 };
    Scalar fairhead = 0;
    for( int power = 4; power >= 0; --power )
    {
        Scalar sum = 0;
        for( int j = edges[ power + 1 ] - 1; j >= edges[ power ]; --j )
        {
            sum += Scalar( fairhd[ j ][ 0 ] ) * sin( Scalar( fairhd[ j ][ 1 ] ) * t + Scalar( fairhd[ j ][ 2 ] ) );
        }
        fairhead = fairhead * t + sum;
    }
    const Scalar massAdjustment = Scalar( 0.00065e-6 ) * sin( Scalar( 6069.776754 ) * t + Scalar( 4.021194 ) ) +
            Scalar( 0.00033e-6 ) * sin( Scalar( 213.299095 ) * t + Scalar( 5.543132 ) ) -
            Scalar( 0.00196e-6 ) * sin( Scalar( 6208.294251 ) * t + Scalar( 5.696701 ) ) -
            Scalar( 0.00173e-6 ) * sin( Scalar( 74.781599 ) * t + Scalar( 2.435900 ) ) + Scalar( 0.03638e-6 ) * t * t;
    return topocentric + fairhead + massAdjustment;
}
#endif

}  // namespace sofa_interface

}  // namespace tudat
