#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <sstream>
#include "tudat/io/readHighFrequencyEopData.h"

using namespace tudat::input_output;
namespace
{
const std::string row = "M2 2 0 0 -2 0 -2 255.555 .51752505 1 -2 3 -4 5 -6 7 -8";
HighFrequencyEopData parse( const std::string& text )
{
    std::istringstream stream( text );
    return readHighFrequencyEopData( stream, "parser fixture" );
}
}  // namespace

BOOST_AUTO_TEST_CASE( neutral_reader_preserves_labels_headers_and_converts_all_components )
{
    const auto data = parse( "  # retained provenance\n\n" + row + "\r\n" );
    BOOST_CHECK_EQUAL( data.source, "parser fixture" );
    BOOST_REQUIRE_EQUAL( data.header.size( ), 1 );
    BOOST_REQUIRE_EQUAL( data.terms.size( ), 1 );
    const auto& term = data.terms.front( );
    BOOST_CHECK_EQUAL( term.name, "M2" );
    BOOST_CHECK_EQUAL( term.doodsonNumber, "255.555" );
    BOOST_CHECK_EQUAL( term.argumentMultipliers[ 0 ], 2 );
    const double uas = std::acos( -1.0 ) / ( 180.0 * 3600.0 * 1e6 );
    for( int i = 0; i < 4; ++i )
    {
        const double unit = i < 2 ? uas : 1e-6;
        BOOST_CHECK_SMALL( term.sineAmplitudes[ i ] - ( 2 * i + 1 ) * unit, 1e-20 );
        BOOST_CHECK_SMALL( term.cosineAmplitudes[ i ] + ( 2 * i + 2 ) * unit, 1e-20 );
    }
    BOOST_CHECK_EQUAL( parse( "- 1 -1 0 -2 -2 -2 1X9.E55 1 0 0 0 0 0 0 0 0" ).terms[ 0 ].doodsonNumber, "1X9.E55" );
}

BOOST_AUTO_TEST_CASE( malformed_empty_duplicate_and_nonfinite_tables_are_rejected )
{
    for( const auto& text : { "",
                              "# only header",
                              "M2 2 0",
                              "M2 2.5 0 0 -2 0 -2 255.555 1 0 0 0 0 0 0 0 0",
                              "M2 2 0 0 -2 0 -2 255.555 0 0 0 0 0 0 0 0 0",
                              "M2 2 0 0 -2 0 -2 255.555 -1 0 0 0 0 0 0 0 0",
                              "M2 2 0 0 -2 0 -2 255.555 1 nan 0 0 0 0 0 0 0",
                              "M2 2 0 0 -2 0 -2 255.555 1 1e999 0 0 0 0 0 0 0" } )
        BOOST_CHECK_THROW( parse( text ), std::runtime_error );
    BOOST_CHECK_THROW( parse( row + " extra" ), std::runtime_error );
    BOOST_CHECK_THROW( parse( row + "\n" + row ), std::runtime_error );
    BOOST_CHECK_THROW( readHighFrequencyEopData( "/missing/desai-fixture.txt" ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( builtin_table_is_the_complete_ocean_only_supplement )
{
    const auto data = getDesaiSibois2016EopData( );
    BOOST_REQUIRE_EQUAL( data.terms.size( ), 159 );
    int diurnal = 0, semidiurnal = 0;
    for( const auto& term : data.terms )
    {
        diurnal += term.argumentMultipliers[ 0 ] == 1;
        semidiurnal += term.argumentMultipliers[ 0 ] == 2;
    }
    BOOST_CHECK_EQUAL( diurnal, 86 );
    BOOST_CHECK_EQUAL( semidiurnal, 73 );
    BOOST_CHECK_NE( data.source.find( "53398e40b70ccfa411d292e9a7f752a01cbac1be07d4e7fb18906d6481e57dd8" ), std::string::npos );
    BOOST_CHECK_SMALL( data.terms.front( ).sineAmplitudes[ 2 ] - .07e-6, 1e-22 );
    BOOST_CHECK_SMALL( data.terms.front( ).cosineAmplitudes[ 3 ] + .34e-6, 1e-22 );
}
