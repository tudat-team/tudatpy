#include "tudat/io/readOceanPoleTideData.h"
#include "tudat/io/textParsing.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
void requireValidOceanPoleTideGrid( const bool condition, const std::string& message )
{
    if( !condition )
    {
        throw std::runtime_error( "Invalid ocean pole-tide grid: " + message );
    }
}

double parseOceanPoleTideNumber( const std::string& token )
{
    const auto value = text_parsing::parseFiniteNumber( token );
    requireValidOceanPoleTideGrid( value.has_value( ), "nonfinite or malformed number" );
    return *value;
}
}  // namespace

void validateOceanPoleTideGrid( const OceanPoleTideGrid& grid )
{
    requireValidOceanPoleTideGrid( grid.longitudeCount >= 2 && grid.latitudeCount >= 2, "each axis requires at least two nodes" );
    requireValidOceanPoleTideGrid( grid.longitudeCount <= std::numeric_limits< std::size_t >::max( ) / grid.latitudeCount,
                                   "grid size overflow" );
    requireValidOceanPoleTideGrid( grid.coefficients.size( ) == grid.longitudeCount * grid.latitudeCount, "missing or excess nodes" );
    requireValidOceanPoleTideGrid( std::isfinite( grid.firstLongitudeDegrees ) && std::isfinite( grid.firstLatitudeDegrees ) &&
                                           std::isfinite( grid.longitudeStepDegrees ) && std::isfinite( grid.latitudeStepDegrees ),
                                   "nonfinite geometry" );
    requireValidOceanPoleTideGrid( grid.longitudeStepDegrees > 0.0 && grid.latitudeStepDegrees > 0.0, "steps must be positive" );
    const double lastLatitude = grid.firstLatitudeDegrees + ( grid.latitudeCount - 1 ) * grid.latitudeStepDegrees;
    const double longitudeSpan = ( grid.longitudeCount - 1 ) * grid.longitudeStepDegrees;
    requireValidOceanPoleTideGrid( grid.firstLatitudeDegrees >= -90.0 && lastLatitude <= 90.0, "latitude outside [-90,90]" );
    requireValidOceanPoleTideGrid( longitudeSpan < 360.0 && grid.longitudeCount * grid.longitudeStepDegrees <= 360.0 + 1.0e-9,
                                   "longitude axis exceeds one period or repeats an endpoint" );
    requireValidOceanPoleTideGrid( std::isfinite( grid.firstLongitudeDegrees + longitudeSpan ), "longitude geometry overflow" );
    for( const auto& value : grid.coefficients )
    {
        requireValidOceanPoleTideGrid( value.allFinite( ), "nonfinite coefficient" );
    }
    requireValidOceanPoleTideGrid( grid.referenceOrigin == OceanLoadingReferenceOrigin::unspecified ||
                                           grid.referenceOrigin == OceanLoadingReferenceOrigin::solidEarthCentre ||
                                           grid.referenceOrigin == OceanLoadingReferenceOrigin::wholeEarthCentreOfMass,
                                   "unknown origin" );
}

OceanPoleTideGrid readOceanPoleTideGrid( std::istream& stream, const std::string& source )
{
    std::size_t lineNumber = 0;
    try
    {
        OceanPoleTideGrid grid;
        grid.source = source;
        std::map< std::string, double > geometry;
        bool columnsSeen = false;
        bool unitsSeen = false;
        bool dataStarted = false;
        std::string line;
        while( std::getline( stream, line ) )
        {
            ++lineNumber;
            if( !line.empty( ) && line.back( ) == '\r' )
            {
                line.pop_back( );
            }
            std::istringstream words( line );
            std::string first;
            if( !( words >> first ) )
            {
                continue;
            }
            const bool numeric = first.front( ) == '+' || first.front( ) == '-' || first.front( ) == '.' ||
                    ( first.front( ) >= '0' && first.front( ) <= '9' );

            // The provider separator line starts with '-', but contains no number.
            if( numeric && first.find_first_not_of( '-' ) != std::string::npos )
            {
                requireValidOceanPoleTideGrid( columnsSeen && unitsSeen, "missing component order or units header" );
                if( !dataStarted )
                {
                    const std::array< std::string, 8 > keys = { "Number_longitude_Grid_Points", "First_longitude_degrees",
                                                                "Last_longitude_degrees",       "Longitude_step_degrees",
                                                                "Number_latitude_grid_points",  "First_latitude_degrees",
                                                                "Last_latitude_degrees",        "Latitude_step_degrees" };
                    for( const auto& key : keys )
                    {
                        requireValidOceanPoleTideGrid( geometry.count( key ) != 0, "missing geometry field " + key );
                    }
                    const double nlon = geometry.at( keys[ 0 ] ), nlat = geometry.at( keys[ 4 ] );
                    requireValidOceanPoleTideGrid( nlon >= 2 && nlat >= 2 && nlon <= 1000000 && nlat <= 1000000 &&
                                                           nlon == std::floor( nlon ) && nlat == std::floor( nlat ),
                                                   "invalid node count" );
                    grid.longitudeCount = static_cast< std::size_t >( nlon );
                    grid.latitudeCount = static_cast< std::size_t >( nlat );
                    grid.firstLongitudeDegrees = geometry.at( keys[ 1 ] );
                    grid.longitudeStepDegrees = geometry.at( keys[ 3 ] );
                    grid.firstLatitudeDegrees = geometry.at( keys[ 5 ] );
                    grid.latitudeStepDegrees = geometry.at( keys[ 7 ] );
                    requireValidOceanPoleTideGrid( std::abs( geometry.at( keys[ 2 ] ) - grid.firstLongitudeDegrees -
                                                             ( nlon - 1 ) * grid.longitudeStepDegrees ) < 1e-9 &&
                                                           std::abs( geometry.at( keys[ 6 ] ) - grid.firstLatitudeDegrees -
                                                                     ( nlat - 1 ) * grid.latitudeStepDegrees ) < 1e-9,
                                                   "last coordinates disagree with counts and steps" );
                    requireValidOceanPoleTideGrid( grid.longitudeStepDegrees > 0 && grid.latitudeStepDegrees > 0, "nonpositive steps" );
                    dataStarted = true;
                }
                requireValidOceanPoleTideGrid( grid.coefficients.size( ) < grid.longitudeCount * grid.latitudeCount, "excess nodes" );
                std::array< double, 8 > row;
                row[ 0 ] = parseOceanPoleTideNumber( first );
                std::string token;
                for( std::size_t i = 1; i < row.size( ); ++i )
                {
                    requireValidOceanPoleTideGrid( static_cast< bool >( words >> token ), "truncated eight-column row" );
                    row[ i ] = parseOceanPoleTideNumber( token );
                }
                requireValidOceanPoleTideGrid( !( words >> token ), "extra columns" );
                const std::size_t node = grid.coefficients.size( );
                const double expectedLongitude = grid.firstLongitudeDegrees + ( node % grid.longitudeCount ) * grid.longitudeStepDegrees;
                const double expectedLatitude = grid.firstLatitudeDegrees + ( node / grid.longitudeCount ) * grid.latitudeStepDegrees;
                requireValidOceanPoleTideGrid(
                        std::abs( row[ 0 ] - expectedLongitude ) < 1e-9 && std::abs( row[ 1 ] - expectedLatitude ) < 1e-9,
                        "duplicate, missing, reordered or off-grid node" );
                Eigen::Matrix< double, 3, 2 > coefficients;
                for( int component = 0; component < 3; ++component )
                {
                    coefficients( component, 0 ) = row[ 2 + 2 * component ];
                    coefficients( component, 1 ) = row[ 3 + 2 * component ];
                }
                grid.coefficients.push_back( coefficients );
            }
            else
            {
                requireValidOceanPoleTideGrid( !dataStarted, "unexpected text after grid data" );
                grid.sourceHeaders.push_back( line );
                if( line.find( "CM correction applied" ) != std::string::npos )
                {
                    grid.referenceOrigin = OceanLoadingReferenceOrigin::wholeEarthCentreOfMass;
                }
                if( first == "Longitude" )
                {
                    requireValidOceanPoleTideGrid( !columnsSeen, "duplicate component header" );
                    const std::array< std::string, 7 > names = { "Latitude", "u_r^R", "u_r^I", "u_n^R", "u_n^I", "u_e^R", "u_e^I" };
                    for( const auto& name : names )
                    {
                        std::string actual;
                        requireValidOceanPoleTideGrid( static_cast< bool >( words >> actual ) && actual == name, "wrong component order" );
                    }
                    std::string extra;
                    requireValidOceanPoleTideGrid( !( words >> extra ), "extra component header" );
                    columnsSeen = true;
                }
                else if( first == "(degrees)" )
                {
                    requireValidOceanPoleTideGrid( !unitsSeen, "duplicate units header" );
                    std::string token;
                    requireValidOceanPoleTideGrid( static_cast< bool >( words >> token ) && token == "(degrees)",
                                                   "wrong coordinate units" );
                    std::string remainder;
                    std::getline( words, remainder );
                    remainder.erase( std::remove_if( remainder.begin( ), remainder.end( ), []( char c ) { return c == ' ' || c == '\t'; } ),
                                     remainder.end( ) );
                    requireValidOceanPoleTideGrid( remainder == "()()()()()()", "six coefficient units must be dimensionless" );
                    unitsSeen = true;
                }
                else if( line.find( '=' ) != std::string::npos )
                {
                    std::string equals, token, extra;
                    requireValidOceanPoleTideGrid( static_cast< bool >( words >> equals >> token ) && equals == "=" && !( words >> extra ),
                                                   "malformed geometry header" );
                    requireValidOceanPoleTideGrid( geometry.emplace( first, parseOceanPoleTideNumber( token ) ).second,
                                                   "duplicate geometry header" );
                }
            }
        }
        requireValidOceanPoleTideGrid( !stream.bad( ) && !( stream.fail( ) && !stream.eof( ) ), "stream read failure" );
        validateOceanPoleTideGrid( grid );
        return grid;
    }
    catch( const std::exception& error )
    {
        throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": " + error.what( ) );
    }
}

OceanPoleTideGrid readOceanPoleTideGrid( const std::string& filePath )
{
    std::ifstream stream( filePath );
    if( !stream )
    {
        throw std::runtime_error( "Cannot open ocean pole-tide grid: " + filePath );
    }
    return readOceanPoleTideGrid( stream, filePath );
}
}  // namespace input_output
}  // namespace tudat
