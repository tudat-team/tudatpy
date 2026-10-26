/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/basics/utilities.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/interface/sofa/earthOrientation.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "boost/algorithm/string.hpp"

#include <cmath>
#include <cstdio>

namespace tudat
{

namespace ephemerides
{

Eigen::Matrix3d getRotationMatrixFromTemeToJ2000( const double epochSinceJ2000 )
{
    // First, rotate to the True Of Date (TOD) frame.
    double equationOfEquinoxes = sofa_interface::calculateEquationOfEquinoxes( epochSinceJ2000 );

    // Rotate around pole (z-axis)
    Eigen::AngleAxisd rotationObject = Eigen::AngleAxisd( equationOfEquinoxes, Eigen::Vector3d::UnitZ( ) );
    Eigen::Matrix3d rotationMatrix1 = rotationObject.toRotationMatrix( );

    // Now that we have our state vector in the TOD frame, we need to obtain the combined precession + nutation matrix from Sofa
    // (according to the 1976/1980 model)
    Eigen::Matrix3d precessionNutationMatrix = sofa_interface::getPrecessionNutationMatrix( epochSinceJ2000 );
    Eigen::Matrix3d rotationMatrix2 = precessionNutationMatrix.transpose( );

    return rotationMatrix2 * rotationMatrix1;
}

Eigen::Matrix3d getRotationMatrixFromJ2000ToTeme( const double epochSinceJ2000 )
{
    return getRotationMatrixFromTemeToJ2000( epochSinceJ2000 ).transpose( );
}

Eigen::Matrix3d getRotationMatrixFromTemeToFrame( const double epochSinceJ2000, const std::string& frameOrientation )
{
    const Eigen::Matrix3d temeToJ2000 = getRotationMatrixFromTemeToJ2000( epochSinceJ2000 );
    if( frameOrientation == "J2000" )
    {
        return temeToJ2000;
    }
    if( frameOrientation == "ECLIPJ2000" )
    {
        return spice_interface::getRotationFromJ2000ToEclipJ2000( ) * temeToJ2000;
    }
    throw std::runtime_error( "TLE state conversion to target frame " + frameOrientation + " is currently unsupported." );
}

Eigen::Matrix3d getRotationMatrixFromFrameToTeme( const double epochSinceJ2000, const std::string& frameOrientation )
{
    return getRotationMatrixFromTemeToFrame( epochSinceJ2000, frameOrientation ).transpose( );
}

TleEphemeris::TleEphemeris( const std::string& referenceFrameOrigin,
                            const std::string& referenceFrameOrientation,
                            const std::shared_ptr< Tle > tle_ptr,
                            const bool useSDP ): Ephemeris( referenceFrameOrigin, referenceFrameOrientation )
{
    if( referenceFrameOrigin != "Earth" )
    {
        throw std::runtime_error( "Error: TleEphemeris only supports an Earth-centered reference frame." );
    }
    if( useSDP )
    {
        throw std::runtime_error( "TLE SDP propagator requested, which is not yet implemented." );
    }
    useSDP_ = useSDP;
    tle_ = tle_ptr;
    // Create table with ephemeris data
    // Create interpolator based on table
}

Eigen::Vector6d TleEphemeris::getCartesianStateInTemeFrame( double secondsSinceEpoch )
{
    // Call Spice interface to retrieve the spacecraft's state from the TLE in the True Equator, Mean Equinox frame (see
    // Vallado: Fundamentals of Astrodynamics and Applications 4th ed. (2013)). This frame is idiosyncratic in nature and
    // therefore needs to be converted to an intermediate standard reference frame.
    return spice_interface::getCartesianStateFromTleAtEpoch( secondsSinceEpoch, tle_ );
}

Eigen::Vector6d TleEphemeris::getCartesianState( double secondsSinceEpoch )
{
    // Compute TEME state
    Eigen::Vector6d cartesianStateAtEpochTEME = getCartesianStateInTemeFrame( secondsSinceEpoch );

    // Apply the same centralized TEME-to-output-frame rotation to position and velocity.
    const Eigen::Matrix3d rotationToOutputFrame = getRotationMatrixFromTemeToFrame( secondsSinceEpoch, referenceFrameOrientation_ );
    Eigen::Vector6d cartesianStateOutput = Eigen::Vector6d::Zero( );
    cartesianStateOutput.segment( 0, 3 ) = rotationToOutputFrame * cartesianStateAtEpochTEME.segment( 0, 3 );
    cartesianStateOutput.segment( 3, 3 ) = rotationToOutputFrame * cartesianStateAtEpochTEME.segment( 3, 3 );

    return cartesianStateOutput;
}

// Computes the TLE checksum (0-9) for a 69-character line (checksum digit is at index 68).
int computeTleChecksum( const std::string& line )
{
    if( line.size( ) != 69 )
    {
        throw std::runtime_error( "TLE line must be exactly 69 characters." );
    }

    int checksum = 0;
    for( std::size_t i = 0; i < 68; ++i )
    {
        const char c = line[ i ];
        if( c == '-' )
        {
            ++checksum;
        }
        else if( std::isdigit( static_cast< unsigned char >( c ) ) )
        {
            checksum += static_cast< int >( c ) - '0';
        }
    }

    return checksum % 10;
}

void validateTleChecksumOrThrow( const std::string& line )
{
    const int expected = computeTleChecksum( line );
    const int provided = static_cast< int >( line.at( 68 ) ) - '0';

    if( expected != provided )
    {
        throw std::runtime_error( "Invalid TLE checksum: expected " + std::to_string( expected ) + ", got " + std::to_string( provided ) );
    }
}

namespace
{

//! Format a value as a fixed-width, space-padded field with the given number of decimals (e.g. degrees fields).
std::string formatTleFixedField( const double value, const int totalWidth, const int fractionDigits )
{
    char buffer[ 32 ];
    std::snprintf( buffer, sizeof( buffer ), "%*.*f", totalWidth, fractionDigits, value );
    return std::string( buffer );
}

//! Format a value in (-1, 1) without the leading zero, e.g. " .00012345" or "-.00012345" (mean motion first derivative).
std::string formatTleSignedDecimalField( const double value, const int fractionDigits )
{
    const char sign = ( value < 0.0 ) ? '-' : ' ';
    char buffer[ 32 ];
    std::snprintf( buffer, sizeof( buffer ), "%.*f", fractionDigits, std::fabs( value ) );
    const std::string digitsString( buffer );
    return std::string( 1, sign ) + digitsString.substr( digitsString.find( '.' ) );
}

//! Format a value using the TLE mantissa/exponent convention: value = sign * mantissa * 10^(exponent - 5),
//! with mantissa a 5-digit integer (columns: sign + 5 digits) and exponent a signed single digit.
//! Used for the mean motion second derivative and B* fields.
std::string formatTleExponentialField( const double value )
{
    const char sign = ( value < 0.0 ) ? '-' : ' ';
    const double absoluteValue = std::fabs( value );

    int exponent = 0;
    int mantissa = 0;
    if( absoluteValue > 0.0 )
    {
        exponent = static_cast< int >( std::floor( std::log10( absoluteValue ) ) ) + 1;
        double mantissaValue = absoluteValue / std::pow( 10.0, exponent );
        mantissa = static_cast< int >( std::llround( mantissaValue * 1.0e5 ) );
        // Guard against rounding a mantissa such as 99999.6 up to 100000 (would overflow the 5-digit field).
        if( mantissa >= 100000 )
        {
            mantissa /= 10;
            exponent += 1;
        }
    }
    if( std::abs( exponent ) > 9 )
    {
        throw std::runtime_error( "Cannot represent value " + std::to_string( value ) +
                                  " in a TLE mantissa/exponent field: exponent magnitude exceeds the single allotted digit." );
    }

    char buffer[ 16 ];
    std::snprintf( buffer, sizeof( buffer ), "%c%05d%+d", sign, mantissa, exponent );
    return std::string( buffer );
}

}  // namespace

void Tle::computeRawLines( )
{
    // Convert the stored TDB epoch back to UTC to recover the calendar date used by the TLE format. This mirrors,
    // in reverse, the UTC-to-TT-to-(approximate)TDB conversion performed by the string-parsing constructor.
    const double ttEpoch = sofa_interface::convertTDBtoTT< double >( epoch_, Eigen::Vector3d::Zero( ) );
    const double utcEpoch = sofa_interface::convertTTtoUTC< double >( ttEpoch );

    const basic_astrodynamics::DateTime epochDateTime = basic_astrodynamics::DateTime::fromTime< double >( utcEpoch );
    const int epochYear = epochDateTime.getYear( );
    const basic_astrodynamics::DateTime startOfYear( epochYear, 1, 1, 0, 0, 0.0L );
    // TLE day numbering starts at 1 for January 1st, 00:00.
    const double epochDayFraction = ( utcEpoch - startOfYear.epoch< double >( ) ) / 86400.0 + 1.0;
    char dayFractionBuffer[ 16 ];
    std::snprintf( dayFractionBuffer, sizeof( dayFractionBuffer ), "%012.8f", epochDayFraction );

    std::string pieceField = internationalDesignatorPiece_;
    pieceField.resize( 3, ' ' );

    // Physical rates are stored in radians/minute^{1,2,3}; the TLE fields are in revolutions/day^{1,2,3}.
    const double revPerDayToRadPerMin = 2.0 * mathematical_constants::PI / 1440.0;

    // ---- Line 1: identification, epoch and drag/ballistic terms ----
    std::string line1 = "1 " + utilities::paddedZeroIntString( noradCatalogNumber_, 5 ) + std::string( 1, classification_ ) + " " +
            utilities::paddedZeroIntString( internationalDesignatorLaunchYear_ % 100, 2 ) +
            utilities::paddedZeroIntString( internationalDesignatorLaunchNumber_, 3 ) + pieceField + " " +
            utilities::paddedZeroIntString( epochYear % 100, 2 ) + std::string( dayFractionBuffer ) + " " +
            formatTleSignedDecimalField( meanMotionFirstDerivative_ / revPerDayToRadPerMin, 8 ) + " " +
            formatTleExponentialField( meanMotionSecondDerivative_ / ( revPerDayToRadPerMin / 1440.0 / 1440.0 ) ) + " " +
            formatTleExponentialField( bStar_ ) + " " + std::to_string( ephemerisType_ ) + " " +
            utilities::paddedZeroIntString( elementSetNumber_, 4 );
    line1 += std::to_string( computeTleChecksum( line1 + "0" ) );

    // ---- Line 2: mean orbital elements ----
    std::string line2 = "2 " + utilities::paddedZeroIntString( noradCatalogNumber_, 5 ) + " " +
            formatTleFixedField( unit_conversions::convertRadiansToDegrees( inclination_ ), 8, 4 ) + " " +
            formatTleFixedField( unit_conversions::convertRadiansToDegrees( rightAscension_ ), 8, 4 ) + " " +
            utilities::paddedZeroIntString( static_cast< int >( std::llround( eccentricity_ * 1.0e7 ) ), 7 ) + " " +
            formatTleFixedField( unit_conversions::convertRadiansToDegrees( argOfPerigee_ ), 8, 4 ) + " " +
            formatTleFixedField( unit_conversions::convertRadiansToDegrees( meanAnomaly_ ), 8, 4 ) + " " +
            formatTleFixedField( meanMotion_ / revPerDayToRadPerMin, 11, 8 ) +
            utilities::paddedZeroIntString( revolutionNumberAtEpoch_, 5 );
    line2 += std::to_string( computeTleChecksum( line2 + "0" ) );

    rawLine1_ = line1;
    rawLine2_ = line2;
}

Tle::Tle( const std::string& lines )
{
    // First of all, the TLE lines to be checked for validity: they shall not be empty, contain more than 69 characters
    // or have an invalid checksum.
    // Although unlikely, TLEs can become corrupted due to (un)intentional changes or data corruption during transmission.

    // Check if string is not empty
    if( lines.empty( ) )
    {
        throw std::runtime_error( "Error: TLE class was instantiated with string object, but string is empty." );
    }
    std::vector< std::string > tleLines;
    boost::algorithm::split( tleLines, lines, boost::is_any_of( "\n\r" ) );
    if( tleLines.size( ) != 2 )
    {
        throw std::runtime_error( "Error: TLE class was instantiated with string object, but string does not contain two 2 lines." );
    }
    // Check line length
    for( std::string line : tleLines )
    {
        if( line.length( ) != 69 )
        {
            std::cout << line << std::endl;
            std::string s( 1, line.at( 0 ) );
            throw std::runtime_error( "Error: TLE class was instantiated with string object, but line " + s +
                                      " contains an invalid number of characters (counted " + std::to_string( line.length( ) ) +
                                      " characters)" );
        }
        // Checksum
        int checksum = 0;
        for( unsigned int i = 0; i < line.length( ) - 1; i++ )
        {
            char c = line[ i ];
            if( c == '-' )
            {
                checksum++;
            }
            else if( std::isdigit( c ) )
            {
                checksum += static_cast< int >( c ) - 48;
            }
            // Else do nothing
        }
        checksum %= 10;
        if( checksum == static_cast< int >( line[ 68 ] ) - 48 )
        {
            // std::cout << "TLE checksum verified." << std::endl;
        }
        else
        {
            throw std::runtime_error( "TLE checksum invalid!" );
        }
    }

    // Parse TLE
    std::string line1 = tleLines.at( 0 );
    std::string line2 = tleLines.at( 1 );

    rawLine1_ = line1;
    rawLine2_ = line2;

    // Local-only metadata
    const int lineNumber1 = std::stoi( line1.substr( 0, 1 ) );
    const int lineNumber2 = std::stoi( line2.substr( 0, 1 ) );

    if( lineNumber1 != 1 || lineNumber2 != 2 )
    {
        throw std::runtime_error( "Error, line numbers incorrect in TLE: \n " + line1 + "\n" + line2 );
    }
    validateTleChecksumOrThrow( line1 );
    validateTleChecksumOrThrow( line2 );

    noradCatalogNumber_ = std::stoi( line1.substr( 2, 5 ) );
    classification_ = line1.at( 7 );

    // IF missing INTLDES, replace with mock designation
    std::string launchYearStr = boost::algorithm::trim_copy( line1.substr( 9, 2 ) );
    internationalDesignatorLaunchYear_ = launchYearStr.empty( ) ? 0 : std::stoi( launchYearStr );
    std::string launchNumberStr = boost::algorithm::trim_copy( line1.substr( 11, 3 ) );
    internationalDesignatorLaunchNumber_ = launchNumberStr.empty( ) ? 0 : std::stoi( launchNumberStr );
    internationalDesignatorPiece_ =
            boost::algorithm::trim_copy( line1.substr( 14, 3 ) ).empty( ) ? "A" : boost::algorithm::trim_copy( line1.substr( 14, 3 ) );

    int epochYear = std::stoi( line1.substr( 18, 2 ) );
    double epochDayFraction = std::stod( line1.substr( 20, 12 ) );

    // Convert to seconds since J2000
    // TLE day number starts with a 1, so a day fraction of 1.0 would mean Jan 1st, 0:00. Hence, we have to subtract 1.5 days
    // in seconds to obtain number of seconds since Jan 1st, noon (as dictated by J2000).
    if( epochYear < 57 )
    {
        epochYear += 2000;
    }
    else
    {
        epochYear += 1900;
    }
    // TLE day numbering starts with 1, whereas Tudat assumes January 1st to be number 0
    referenceDate_ = basic_astrodynamics::DateTime( epochYear, 1, 1, 0, 0, 0.0 );
    referenceDate_ = referenceDate_.addDaysToDateTime( epochDayFraction - 1.0 );
    epoch_ = referenceDate_.epoch< Time >( ).getSeconds< double >( );

    epoch_ = sofa_interface::convertUTCtoTT< double >( epoch_ );
    epoch_ += sofa_interface::getTDBminusTT( epoch_, Eigen::Vector3d::Zero( ) );

    // Mean motion first derivative (rev/day^2)
    meanMotionFirstDerivative_ = std::stod( line1.substr( 33, 10 ) );
    meanMotionFirstDerivative_ *= 2.0 * mathematical_constants::PI / ( ( 60.0 * 24.0 ) * ( 60.0 * 24.0 ) );
    // Mean motion second derivative (rev/day^3), mantissa+exponent format
    const double nDdotMantissa = std::stod( line1.substr( 44, 6 ) );
    const int nDdotExponent = std::stoi( line1.substr( 50, 2 ) );
    meanMotionSecondDerivative_ = nDdotMantissa * 1.0e-5 * std::pow( 10.0, nDdotExponent );
    meanMotionSecondDerivative_ *= 2.0 * mathematical_constants::PI / ( ( 60.0 * 24.0 ) * ( 60.0 * 24.0 ) * ( 60.0 * 24.0 ) );

    double bStar = std::stod( line1.substr( 53, 6 ) );
    double bStarExp = std::stod( line1.substr( 59, 2 ) );
    bStar_ = bStar * std::pow( 10, bStarExp - 5 );

    ephemerisType_ = std::stoi( line1.substr( 62, 1 ) );
    elementSetNumber_ = std::stoi( line1.substr( 64, 4 ) );

    // ---- Line 2: classical elements ----
    inclination_ = unit_conversions::convertDegreesToRadians( std::stod( line2.substr( 8, 8 ) ) );
    rightAscension_ = unit_conversions::convertDegreesToRadians( std::stod( line2.substr( 17, 8 ) ) );

    std::string eccentricityString = "0." + line2.substr( 26, 7 );
    eccentricity_ = std::stod( eccentricityString );

    argOfPerigee_ = unit_conversions::convertDegreesToRadians( std::stod( line2.substr( 34, 8 ) ) );
    meanAnomaly_ = unit_conversions::convertDegreesToRadians( std::stod( line2.substr( 43, 8 ) ) );

    // Convert mean motion to radians per min (as opposed to rev/day as given in TLE)
    double meanMotionRevDay = std::stod( line2.substr( 52, 11 ) );
    meanMotion_ = meanMotionRevDay * 2.0 * mathematical_constants::PI / ( 60.0 * 24.0 );

    revolutionNumberAtEpoch_ = std::stoi( line2.substr( 63, 5 ) );
}

Tle::Tle( const std::string& tleLine1, const std::string& tleLine2 ): Tle( std::string( tleLine1 + "\n" + tleLine2 ) ) {}

Tle::Tle( const double* spiceElements )
{
    for( double& spiceElement : spiceElements_ )
    {
        spiceElement = *spiceElements;
        spiceElements++;
    }
}

}  // namespace ephemerides
}  // namespace tudat
