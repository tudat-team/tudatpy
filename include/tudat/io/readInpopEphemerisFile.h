#ifndef TUDAT_READINPOPTIMEEPHEMERIS_H
#define TUDAT_READINPOPTIMEEPHEMERIS_H

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <iomanip>

//#include "tudat/math/basic/linearAlgebraTypes.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/math/interpolators/chebyshevApproximation.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"

namespace tudat
{

namespace input_output
{

// Utility: Trim whitespace from a string
inline void trimString( std::string& str )
{
    str.erase( 0, str.find_first_not_of( " \t\r\n" ) );
    str.erase( str.find_last_not_of( " \t\r\n" ) + 1 );
}

// Utility: Split string into tokens
inline std::vector< std::string > splitString( const std::string& str )
{
    std::istringstream iss( str );
    std::vector< std::string > tokens;
    std::string token;
    while ( iss >> token )
    {
        tokens.push_back( token );
    }
    return tokens;
}

// Templated function to read INPOP time ephemeris Chebyshev coefficients
template< typename DifferenceScalarType >
std::map< double, std::vector< DifferenceScalarType > > readInpopTimeEphemeris(
        const std::string& fileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 )
{
    std::map< double, std::vector< DifferenceScalarType > > chebyshevCoefficients;
    std::ifstream stream( fileName );
    if ( !stream.is_open( ) )
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

    while ( std::getline( stream, line ) )
    {
        trimString( line );

        // Convert D-format to E-format
        std::replace( line.begin( ), line.end( ), 'D', 'E' );
        std::replace( line.begin( ), line.end( ), 'd', 'e' );

        std::vector< std::string > tokens = splitString( line );

        if ( tokens.size( ) > 1 )
        {
            std::vector< DifferenceScalarType > currentCoefficients( order );
            for ( int i = 0; i < order; ++i )
            {
                std::istringstream( tokens.at( i + 2 ) ) >> currentCoefficients[ i ];
            }

            double timeSinceJ2000 = ( std::stod( tokens.at( 0 ) ) - referenceJulianDay ) *
                                    physical_constants::JULIAN_DAY;
            chebyshevCoefficients[ timeSinceJ2000 ] = currentCoefficients;
        }
    }


    return chebyshevCoefficients;
}

// Reads INPOP state ephemerides (defined in .cpp)
std::map< double, std::vector< Eigen::Vector6d > > readInpopStateEphemeris(
        const std::string& positionFileName,
        const std::string& velocityFileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

// Creates a double-precision interpolator from INPOP time ephemeris data
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > createInpopTimeEphemerisInterpolator(
        const std::string& fileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

// Creates a long double-precision interpolator from INPOP time ephemeris data
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, long double > > createLongInpopTimeEphemerisInterpolator(
        const std::string& fileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

// Constructs a full tabulated ephemeris from INPOP position and velocity data
std::shared_ptr< ephemerides::Ephemeris > createInpopEphemerisFromFiles(
        const std::string& positionFileName,
        const std::string& velocityFileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000,
        const bool useGeocentricReference = false );

} // namespace input_output

} // namespace tudat

#endif // TUDAT_READINPOPTIMEEPHEMERIS_H
