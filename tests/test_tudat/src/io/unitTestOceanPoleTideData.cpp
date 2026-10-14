#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <sstream>
#include "tudat/io/readOceanPoleTideData.h"

namespace
{
std::string fixture( )
{
    std::ifstream file( std::string( __FILE__ ).substr( 0, std::string( __FILE__ ).find_last_of( "/\\" ) ) +
                        "/../../reference_files/ocean_pole_tide/iers_grid_extract.txt" );
    BOOST_REQUIRE( file.good( ) );
    return std::string( std::istreambuf_iterator< char >( file ), std::istreambuf_iterator< char >( ) );
}
void rejects( const std::string& text )
{
    std::istringstream stream( text );
    BOOST_CHECK_THROW( tudat::input_output::readOceanPoleTideGrid( stream ), std::exception );
}
std::string replace( std::string text, const std::string& from, const std::string& to )
{
    const auto position = text.find( from );
    BOOST_REQUIRE( position != std::string::npos );
    text.replace( position, from.size( ), to );
    return text;
}
}  // namespace

BOOST_AUTO_TEST_CASE( providerGridExtractAndMetadata )
{
    using namespace tudat::input_output;
    const auto text = fixture( );
    std::istringstream input( text );
    const auto grid = readOceanPoleTideGrid( input, "pinned official extract" );
    BOOST_CHECK_EQUAL( grid.longitudeCount, 2 );
    BOOST_CHECK_EQUAL( grid.latitudeCount, 2 );
    BOOST_CHECK_EQUAL( grid.source, "pinned official extract" );
    BOOST_CHECK( grid.referenceOrigin == OceanLoadingReferenceOrigin::wholeEarthCentreOfMass );
    BOOST_CHECK_EQUAL( grid.sourceHeaders.size( ), 14 );
    BOOST_CHECK_EQUAL( grid.coefficients[ 0 ]( 0, 0 ), 0.216192 );
    BOOST_CHECK_EQUAL( grid.coefficients[ 0 ]( 0, 1 ), 0.285652 );
    BOOST_CHECK_EQUAL( grid.coefficients[ 0 ]( 1, 0 ), 0.012697 );
    BOOST_CHECK_EQUAL( grid.coefficients[ 0 ]( 2, 1 ), 0.010389 );
    std::string crlf;
    for( char c : text )
    {
        if( c == '\n' )
        {
            crlf += '\r';
        }
        crlf += c;
    }
    crlf = replace( crlf, "0.216192", "2.16192D-1" );
    std::istringstream windows( crlf );
    BOOST_CHECK_SMALL( ( readOceanPoleTideGrid( windows ).coefficients[ 0 ] - grid.coefficients[ 0 ] ).norm( ), 1e-15 );
}

BOOST_AUTO_TEST_CASE( rejectsMalformedIncompleteDuplicateAndWrongConventions )
{
    const auto text = fixture( );
    rejects( text.substr( 0, text.find_last_of( '\n', text.size( ) - 2 ) ) );
    rejects( text + text.substr( text.find_last_of( '\n', text.size( ) - 2 ) + 1 ) );
    rejects( replace( text, "232.75     -43.75", "232.25     -43.75" ) );
    rejects( replace( text, "0.216192", "nan" ) );
    rejects( replace( text, "u_n^R", "u_e^R" ) );
    rejects( replace( text, "(        )", "(metres)" ) );
    rejects( replace( text, "Last_latitude_degrees        = -43.25", "Last_latitude_degrees        = -42.75" ) );
    rejects( replace( text, "Latitude_step_degrees        =      0.50", "Latitude_step_degrees        = -0.50" ) );
    rejects( text + "garbage\n" );
    rejects( "" );
    BOOST_CHECK_THROW( tudat::input_output::readOceanPoleTideGrid( "/nonexistent/ocean/pole/grid.txt" ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( validatesNeutralPreloadedRecords )
{
    std::istringstream input( fixture( ) );
    auto grid = tudat::input_output::readOceanPoleTideGrid( input );
    grid.firstLatitudeDegrees = 90;
    BOOST_CHECK_THROW( tudat::input_output::validateOceanPoleTideGrid( grid ), std::runtime_error );
    grid.firstLatitudeDegrees = -43.75;
    grid.coefficients.pop_back( );
    BOOST_CHECK_THROW( tudat::input_output::validateOceanPoleTideGrid( grid ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( parseErrorsIncludeProviderSourceAndLine )
{
    for( const auto& invalid : { replace( fixture( ), "0.216192", "broken" ), replace( fixture( ), "0.216192", "nan" ) } )
    {
        std::istringstream input( invalid );
        try
        {
            tudat::input_output::readOceanPoleTideGrid( input, "provider-grid.txt" );
            BOOST_FAIL( "Malformed coefficient was accepted" );
        }
        catch( const std::runtime_error& error )
        {
            BOOST_CHECK_NE( std::string( error.what( ) ).find( "provider-grid.txt:15:" ), std::string::npos );
        }
    }
    std::istringstream failed( fixture( ) );
    failed.setstate( std::ios::failbit );
    BOOST_CHECK_THROW( tudat::input_output::readOceanPoleTideGrid( failed, "failed-stream" ), std::runtime_error );
}
