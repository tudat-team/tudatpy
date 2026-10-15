#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include "tudat/io/readStationLoadingData.h"

using namespace tudat::input_output;
namespace
{
const std::string signature = "EPHEDISP  Format version of 2005.06.30";
std::string fixture( const std::string& name = "imls_merra2_extract.eph" )
{
    const std::string file = __FILE__;
    std::ifstream input( file.substr( 0, file.find_last_of( "/\\" ) + 1 ) + "station_loading_fixtures/" + name );
    return { std::istreambuf_iterator< char >( input ), std::istreambuf_iterator< char >( ) };
}
StationLoadingData parse( const std::string& text )
{
    std::istringstream stream( text );
    StationLoadingMetadata metadata;
    metadata.source = "loading parser fixture";
    return readEphedispStationLoadingData( stream, metadata );
}
std::string replace( std::string text, const std::string& old, const std::string& replacement )
{
    const auto position = text.find( old );
    BOOST_REQUIRE_NE( position, std::string::npos );
    text.replace( position, old.size( ), replacement );
    return text;
}
std::string record( const std::string& text, const std::string& prefix )
{
    const auto start = text.find( prefix );
    BOOST_REQUIRE_NE( start, std::string::npos );
    return text.substr( start, text.find( '\n', start ) - start );
}
std::string multiEpoch( )
{
    const auto original = fixture( );
    std::string text = signature + "\nP T 3 S          3 E      4 D          6\n";
    text += "T begin   58849     0.0  2020.01.01-00:00:00\nT end     58849 64800.0  2020.01.01-18:00:00\n";
    text += "T sample     0.25000000000\nA    3000.000000\n";
    for( const auto name : { "00NA", "weuc", "yela" } ) text += record( original, std::string( "S  " ) + name ) + "\n";
    for( int epoch = 1; epoch <= 4; ++epoch )
        for( const auto name : { "00NA", "weuc" } )
        {
            if( ( std::string( name ) == "00NA" && epoch == 4 ) || ( std::string( name ) == "weuc" && epoch == 1 ) ) continue;
            auto line = record( original, "D     1  58849     0.0  2020.01.01-00:00:00  " + std::string( name ) );
            std::ostringstream index;
            index << std::setw( 5 ) << epoch;
            line.replace( 2, 5, index.str( ) );

            // Specification says these epoch fields must be ignored, not checked against the index.
            line.replace( 9, 34, "xxxxx yyyyyyy  informational text!" );
            BOOST_REQUIRE_EQUAL( line.size( ), 80 );
            text += line + "\n";
        }
    return text + signature + "\n";
}
}  // namespace

BOOST_AUTO_TEST_CASE( pinned_real_provider_preserves_exact_authority_basis_and_unknown_inventory )
{
    const auto data = parse( fixture( ) );
    BOOST_REQUIRE_EQUAL( data.sites.size( ), 3 );
    BOOST_REQUIRE_EQUAL( data.epochs.size( ), 1 );
    BOOST_CHECK_EQUAL( data.epochs[ 0 ].modifiedJulianDay, 58849 );
    BOOST_CHECK_EQUAL( data.epochs[ 0 ].secondsOfDay, 0.0L );
    BOOST_CHECK_EQUAL( *data.nominalSampleIntervalSeconds, 0.0L );
    BOOST_CHECK( data.timeScale == StationLoadingTimeScale::tai );
    BOOST_CHECK( data.basis == StationLoadingBasis::geocentricRadialEastNorth );
    BOOST_CHECK( data.metadata.referenceOrigin == StationLoadingOrigin::unspecified );
    BOOST_CHECK( data.metadata.annual == StationLoadingSignalPresence::unknown );
    BOOST_CHECK( data.metadata.components.empty( ) );
    const auto& site = data.sites.at( "weuc" );
    BOOST_CHECK_EQUAL( site.validityRadiusMetres, 3000 );
    BOOST_CHECK_EQUAL( site.cartesianPosition[ 0 ], 4075560.075 );
    BOOST_CHECK_EQUAL( site.cartesianPosition[ 1 ], 931580.295 );
    BOOST_CHECK_EQUAL( site.cartesianPosition[ 2 ], 4801621.958 );
    BOOST_REQUIRE_EQUAL( site.samples.size( ), 1 );
    const auto& sample = site.samples.at( 0 );
    BOOST_CHECK_EQUAL( sample.displacement[ 0 ], -.00510 );
    BOOST_CHECK_EQUAL( sample.displacement[ 1 ], -.00202 );
    BOOST_CHECK_EQUAL( sample.displacement[ 2 ], .00170 );
    BOOST_CHECK_GT( sample.sourceLine, site.sourceLine );
    BOOST_CHECK_EQUAL( sample.sourceRecord.substr( 45, 8 ), "weuc    " );
    BOOST_CHECK_EQUAL( data.declaredEndEpoch->modifiedJulianDay, 58849 );
    BOOST_CHECK_EQUAL( data.informationalBeginEpoch, "2020.01.01-00:00:00" );
}

BOOST_AUTO_TEST_CASE( indexes_define_epochs_and_site_coverage_without_informational_fields )
{
    const auto text = multiEpoch( );
    const auto data = parse( text );
    BOOST_REQUIRE_EQUAL( data.epochs.size( ), 4 );
    BOOST_CHECK_EQUAL( data.epochs.back( ).secondsOfDay, 64800.0L );
    BOOST_CHECK_EQUAL( data.sites.at( "00NA" ).samples.begin( )->first, 0 );
    BOOST_CHECK_EQUAL( data.sites.at( "00NA" ).samples.rbegin( )->first, 2 );
    BOOST_CHECK_EQUAL( data.sites.at( "weuc" ).samples.begin( )->first, 1 );
    BOOST_CHECK_EQUAL( data.sites.at( "weuc" ).samples.rbegin( )->first, 3 );
    BOOST_CHECK( data.sites.at( "yela" ).samples.empty( ) );
    auto arbitrary = replace( fixture( ), "S  00NA    ", "S  A!b@9   " );
    arbitrary = replace( arbitrary, "  00NA    ", "  A!b@9   " );
    const auto siteLine = record( arbitrary, "S  A!b@9" );
    const std::string ignoredCoordinates = "xxxxxxxx xxxxxxxx xxxxxx";
    arbitrary = replace( arbitrary, siteLine, siteLine.substr( 0, 56 ) + ignoredCoordinates );
    const auto labels = parse( arbitrary );
    BOOST_CHECK_EQUAL( labels.sites.at( "A!b@9" ).cartesianPosition[ 0 ], -4073662.420 );
    BOOST_CHECK_EQUAL( labels.sites.at( "A!b@9" ).informationalCoordinates, ignoredCoordinates );
    auto fractional = replace( text, "58849     0.0", "58849     0.1" );
    fractional = replace( fractional, "58849 64800.0", "58849 64800.1" );
    BOOST_CHECK_EQUAL( parse( fractional ).epochs[ 1 ].secondsOfDay, 21600.1L );
    auto crlf = replace( fixture( ), signature + "\n", signature + "\r\n" );
    BOOST_CHECK_EQUAL( parse( crlf ).sites.size( ), 3 );
    auto cr = fixture( );
    for( auto& character : cr )
        if( character == '\n' ) character = '\r';
    BOOST_CHECK_EQUAL( parse( cr ).sites.size( ), 3 );
}

BOOST_AUTO_TEST_CASE( counted_sections_trailers_duplicates_and_internal_gaps_are_strict )
{
    const auto original = fixture( );
    const auto firstSample = record( original, "D     1" );
    const auto firstSite = record( original, "S  00NA" );
    BOOST_CHECK_THROW( parse( fixture( "imls_duplicate_identifier.eph" ) ), std::runtime_error );
    for( const auto& text : { replace( original, "S          3", "S          4" ),
                              replace( original, "D          3", "D          4" ),
                              replace( original, "E      1", "E      0" ),
                              replace( original, "A    3000.000000", "A      -1.000000" ),
                              replace( original, firstSample, firstSample.substr( 0, 79 ) ),
                              replace( original, firstSample, firstSample + "\n" + firstSample ),
                              replace( original, firstSite, firstSite + "\n" + firstSite ),
                              replace( original, "  00NA     -0.00104", "  XXXXX    -0.00104" ),
                              replace( original, "D     1", "D     0" ),
                              replace( original, "D     1", "D     2" ),
                              original.substr( 0, original.rfind( signature ) ),
                              original + "# after trailer\n",
                              "# before signature\n" + original,
                              replace( original, " -0.00104", "      nan" ),
                              replace( original, " -0.00104", "    1e999" ) } )
        BOOST_CHECK_THROW( parse( text ), std::runtime_error );
    const auto multi = multiEpoch( );
    const auto internal = record( multi, "D     2" );
    auto gap = replace( multi, internal + "\n", "" );
    gap = replace( gap, "D          6", "D          5" );
    BOOST_CHECK_THROW( parse( gap ), std::runtime_error );
    BOOST_CHECK_THROW( parse( replace( multi, "T sample     0.25000000000", "T sample     0.00000000000" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parse( replace( multi, "58849 64800.0", "58849 64799.0" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parse( replace( multi, "T sample     0.25000000000", "T sample     0.25000000000\nA    3000.000000" ) ),
                       std::runtime_error );
    BOOST_CHECK_THROW( readEphedispStationLoadingData( "/missing/provider.eph" ), std::runtime_error );
    try
    {
        parse( replace( original, "D     1", "D     0" ) );
        BOOST_FAIL( "invalid index should report its source and line" );
    }
    catch( const std::runtime_error& error )
    {
        BOOST_CHECK_NE( std::string( error.what( ) ).find( "loading parser fixture line " ), std::string::npos );
        BOOST_CHECK_NE( std::string( error.what( ) ).find( "sample index outside grid" ), std::string::npos );
    }
}

BOOST_AUTO_TEST_CASE( generic_preloaded_series_and_explicit_profile_remain_neutral )
{
    auto data = parse( multiEpoch( ) );
    data.nominalSampleIntervalSeconds.reset( );
    data.epochs[ 1 ].secondsOfDay = 12345.1234567890123L;
    data.basis = StationLoadingBasis::cartesian;
    BOOST_CHECK_NO_THROW( validateStationLoadingData( data ) );
    data.epochs[ 1 ] = data.epochs[ 0 ];
    BOOST_CHECK_THROW( validateStationLoadingData( data ), std::runtime_error );
    data = parse( fixture( ) );
    data.sites.at( "weuc" ).samples.at( 0 ).displacement[ 0 ] = std::numeric_limits< double >::infinity( );
    BOOST_CHECK_THROW( validateStationLoadingData( data ), std::runtime_error );
    const auto metadata = imlsMerra2AtmosphericLoadingMetadata( "explicit release identifier" );
    BOOST_CHECK( metadata.referenceOrigin == StationLoadingOrigin::wholeEarthCentreOfMass );
    BOOST_CHECK( metadata.annual == StationLoadingSignalPresence::present );
    BOOST_CHECK( metadata.diurnal == StationLoadingSignalPresence::absent );
    BOOST_CHECK( metadata.mean == StationLoadingSignalPresence::absent );
    BOOST_CHECK_EQUAL( metadata.removedHarmonics.size( ), 8 );
    BOOST_CHECK_EQUAL( metadata.release, "explicit release identifier" );
    BOOST_CHECK( metadata.referenceFrame.empty( ) );
    BOOST_CHECK_THROW( imlsMerra2AtmosphericLoadingMetadata( "" ), std::invalid_argument );
}
