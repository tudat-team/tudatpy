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

#ifndef TUDAT_READINPOPTIMEEPHEMERIS_H
#define TUDAT_READINPOPTIMEEPHEMERIS_H

#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <iomanip>

// #include "tudat/math/basic/linearAlgebraTypes.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/math/interpolators/chebyshevApproximation.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"

namespace tudat
{

namespace input_output
{

//! Utility function to trim leading/trailing whitespace from a string.
/*!
 *  \param str String to be trimmed in-place.
 */
inline void trimString( std::string& str )
{
    const std::string::size_type firstNonWhitespace = str.find_first_not_of( " \t\r\n" );
    if( firstNonWhitespace == std::string::npos )
    {
        str.clear( );
        return;
    }

    str.erase( 0, firstNonWhitespace );
    str.erase( str.find_last_not_of( " \t\r\n" ) + 1 );
}

//! Utility function to split a whitespace-separated line into tokens.
/*!
 *  \param str Input string.
 *  \return List of tokens.
 */
inline std::vector< std::string > splitString( const std::string& str )
{
    std::istringstream iss( str );
    std::vector< std::string > tokens;
    std::string token;
    while( iss >> token )
    {
        tokens.push_back( token );
    }
    return tokens;
}

//! Read INPOP time-ephemeris Chebyshev coefficients from an ASCII file.
/*!
 *  \param fileName Path to INPOP time-ephemeris coefficient file.
 *  \param referenceJulianDay Reference epoch in Julian day used to convert epochs to seconds.
 *  \return Map from segment start epoch (seconds since reference epoch) to Chebyshev coefficients.
 */
template< typename DifferenceScalarType >
std::map< double, std::vector< DifferenceScalarType > > readInpopTimeEphemeris(
        const std::string& fileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 )
{
    std::map< double, std::vector< DifferenceScalarType > > chebyshevCoefficients;
    std::ifstream stream( fileName );
    if( !stream.is_open( ) )
    {
        throw std::runtime_error( "Data file '" + fileName + "' could not be opened." );
    }

    std::string line;
    // Skip first line
    std::getline( stream, line );
    // Read header line with Chebyshev order
    std::getline( stream, line );
    trimString( line );
    std::vector< std::string > headerTokens = splitString( line );

    int order = std::stoi( headerTokens.at( 6 ) );

    while( std::getline( stream, line ) )
    {
        trimString( line );

        // Convert D-format to E-format
        std::replace( line.begin( ), line.end( ), 'D', 'E' );
        std::replace( line.begin( ), line.end( ), 'd', 'e' );

        std::vector< std::string > tokens = splitString( line );

        if( tokens.size( ) > 1 )
        {
            std::vector< DifferenceScalarType > currentCoefficients( order );
            for( int i = 0; i < order; ++i )
            {
                std::istringstream( tokens.at( i + 2 ) ) >> currentCoefficients[ i ];
            }

            double timeSinceJ2000 = ( std::stod( tokens.at( 0 ) ) - referenceJulianDay ) * physical_constants::JULIAN_DAY;
            chebyshevCoefficients[ timeSinceJ2000 ] = currentCoefficients;
        }
    }

    return chebyshevCoefficients;
}

//! Read INPOP Cartesian state ephemeris from position and velocity ASCII files.
/*!
 *  \param positionFileName Path to INPOP position data file.
 *  \param velocityFileName Path to INPOP velocity data file.
 *  \param referenceJulianDay Reference epoch in Julian day used to convert epochs to seconds.
 *  \return Map from epoch (seconds since reference epoch) to state vectors.
 */
std::map< double, std::vector< Eigen::Vector6d > > readInpopStateEphemeris(
        const std::string& positionFileName,
        const std::string& velocityFileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

//! Create a double-precision Chebyshev interpolator from INPOP time-ephemeris data.
/*!
 *  \param fileName Path to INPOP time-ephemeris coefficient file.
 *  \param referenceJulianDay Reference epoch in Julian day used to convert epochs to seconds.
 *  \return Interpolator for time difference values as double.
 */
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > createInpopTimeEphemerisInterpolator(
        const std::string& fileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

//! Create a long-double-precision Chebyshev interpolator from INPOP time-ephemeris data.
/*!
 *  \param fileName Path to INPOP time-ephemeris coefficient file.
 *  \param referenceJulianDay Reference epoch in Julian day used to convert epochs to seconds.
 *  \return Interpolator for time difference values as long double.
 */
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, long double > > createLongInpopTimeEphemerisInterpolator(
        const std::string& fileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

//! Construct tabulated translational ephemeris from INPOP position and velocity files.
/*!
 *  \param positionFileName Path to INPOP position data file.
 *  \param velocityFileName Path to INPOP velocity data file.
 *  \param referenceJulianDay Reference epoch in Julian day used to convert epochs to seconds.
 *  \param useGeocentricReference If true, generated ephemeris is interpreted in geocentric frame.
 *  \return Ephemeris object containing INPOP state history.
 */
std::shared_ptr< ephemerides::Ephemeris > createInpopEphemerisFromFiles(
        const std::string& positionFileName,
        const std::string& velocityFileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000,
        const bool useGeocentricReference = false );

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READINPOPTIMEEPHEMERIS_H
