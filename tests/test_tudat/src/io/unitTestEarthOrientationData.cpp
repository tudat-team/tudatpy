#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/readEarthOrientationData.h"
#include <cmath>
#include <fstream>
#include <sstream>

using namespace tudat::input_output;
namespace
{
constexpr double arcsecond = 3.141592653589793238462643383279502884 / 648000.0;
std::string fixturePath( const std::string& name )
{
    const std::string sourceFile = __FILE__;
    return sourceFile.substr( 0, sourceFile.find_last_of( "/\\" ) + 1 ) + "earth_orientation_fixtures/" + name;
}
std::string fixtureText( const std::string& name )
{
    std::ifstream stream( fixturePath( name ) );
    return { std::istreambuf_iterator< char >( stream ), std::istreambuf_iterator< char >( ) };
}
EarthOrientationData parse( const std::string& text, const EarthOrientationFormat format = EarthOrientationFormat::c04 )
{
    std::istringstream stream( text );
    return readEarthOrientationData( stream, format, "test source" );
}
std::string firstFinalsLine( )
{
    const auto text = fixtureText( "finals2000a.txt" );
    return text.substr( 0, text.find( '\n' ) );
}
}  // namespace

BOOST_AUTO_TEST_CASE( official_legacy_c04_units_metadata_and_raw_coverage )
{
    const auto data = readEarthOrientationData( fixturePath( "legacy_c04.txt" ) );
    BOOST_REQUIRE_EQUAL( data.records.size( ), 3 );
    BOOST_CHECK_EQUAL( data.records.begin( )->first, 37665.0 );
    BOOST_CHECK_EQUAL( data.records.rbegin( )->first, 37667.0 );
    BOOST_CHECK_EQUAL( data.source, fixturePath( "legacy_c04.txt" ) );
    BOOST_REQUIRE( !data.header.empty( ) );
    const auto& record = data.records.begin( )->second;
    BOOST_CHECK_SMALL( record.polarMotion[ 0 ] - ( -0.012700 * arcsecond ), 1e-20 );
    BOOST_CHECK_SMALL( record.polarMotion[ 1 ] - 0.213000 * arcsecond, 1e-20 );
    BOOST_CHECK_EQUAL( record.ut1MinusUtc, 0.0326338 );
    BOOST_REQUIRE( record.lengthOfDay );
    BOOST_CHECK_EQUAL( *record.lengthOfDay, 0.0017230 );
    BOOST_REQUIRE( record.celestialPoleOffsets );
    BOOST_CHECK_EQUAL( ( *record.celestialPoleOffsets )[ 0 ], 0.0 );  // Present zero differs from missing.
    BOOST_REQUIRE( record.celestialPoleOffsetUncertainty );
    BOOST_CHECK_SMALL( ( *record.celestialPoleOffsetUncertainty )[ 0 ] - 0.004774 * arcsecond, 1e-20 );
    BOOST_CHECK( !record.polarMotionRates );
    BOOST_CHECK( record.polarMotionStatus.empty( ) );
}

BOOST_AUTO_TEST_CASE( official_modern_c04_fractional_day_rates_and_uncertainties )
{
    const auto data = parse( fixtureText( "modern_c04.txt" ) );
    BOOST_REQUIRE_EQUAL( data.records.size( ), 3 );
    BOOST_CHECK_EQUAL( data.records.begin( )->first, 45700.5 );  // HH is already part of MJD.
    BOOST_CHECK_EQUAL( data.records.rbegin( )->first, 45702.5 );
    const auto& record = data.records.begin( )->second;
    BOOST_CHECK_SMALL( record.polarMotion[ 0 ] + 0.134064 * arcsecond, 1e-20 );
    BOOST_CHECK_EQUAL( record.ut1MinusUtc, 0.3968898 );
    BOOST_CHECK_EQUAL( *record.lengthOfDay, 0.0014139 );
    BOOST_CHECK_SMALL( ( *record.celestialPoleOffsets )[ 0 ] - 0.002255 * arcsecond, 1e-22 );
    BOOST_CHECK_SMALL( ( *record.celestialPoleOffsets )[ 1 ] + 0.002791 * arcsecond, 1e-22 );
    BOOST_REQUIRE( record.polarMotionRates );
    BOOST_CHECK_SMALL( ( *record.polarMotionRates )[ 0 ] + 0.002644 * arcsecond / 86400.0, 1e-24 );
    BOOST_REQUIRE( record.polarMotionRateUncertainty );
    BOOST_CHECK_SMALL( ( *record.polarMotionRateUncertainty )[ 1 ] - 0.000771 * arcsecond / 86400.0, 1e-24 );
    BOOST_CHECK_EQUAL( *record.ut1MinusUtcUncertainty, 0.0002417 );
    BOOST_CHECK_EQUAL( *record.lengthOfDayUncertainty, 0.0000730 );
    BOOST_CHECK( data.terrestrialFrame.find( "ITRF 2020" ) != std::string::npos );
    BOOST_CHECK( data.terrestrialFrame.find( "u2024" ) == std::string::npos );
}

BOOST_AUTO_TEST_CASE( official_finals_bulletin_a_units_flags_missing_fields_and_placeholders )
{
    const auto data = readEarthOrientationData( fixturePath( "finals2000a.txt" ), EarthOrientationFormat::finals2000a );
    BOOST_REQUIRE_EQUAL( data.records.size( ), 4 );  // Trailing date-only placeholder is not coverage.
    BOOST_CHECK_EQUAL( data.records.begin( )->first, 41684.0 );
    BOOST_CHECK_EQUAL( data.records.rbegin( )->first, 61396.0 );
    const auto& observed = data.records.begin( )->second;
    BOOST_CHECK_EQUAL( observed.polarMotionStatus, "I" );
    BOOST_CHECK_EQUAL( observed.ut1Status, "I" );
    BOOST_CHECK_EQUAL( observed.celestialPoleStatus, "P" );
    BOOST_CHECK_SMALL( observed.polarMotion[ 0 ] - 0.120733 * arcsecond, 1e-20 );
    BOOST_CHECK_EQUAL( observed.ut1MinusUtc, 0.8084178 );  // Not Bulletin B's 0.8075000.
    BOOST_CHECK_SMALL( ( *observed.celestialPoleOffsets )[ 0 ] + 0.766 * arcsecond / 1000.0, 1e-22 );
    BOOST_CHECK_EQUAL( *observed.lengthOfDay, 0.0 );  // Genuine observed zero remains present.
    const auto& next = std::next( data.records.begin( ) )->second;
    BOOST_CHECK_SMALL( *next.lengthOfDay - 0.0035563, 1e-18 );
    const auto& predicted = data.records.rbegin( )->second;
    BOOST_CHECK_EQUAL( predicted.polarMotionStatus, "P" );
    BOOST_CHECK_EQUAL( predicted.ut1Status, "P" );
    BOOST_CHECK( !predicted.lengthOfDay );
    BOOST_CHECK( !predicted.celestialPoleOffsets );
    BOOST_CHECK( !predicted.celestialPoleOffsetUncertainty );
    BOOST_CHECK( predicted.celestialPoleStatus.empty( ) );
}

BOOST_AUTO_TEST_CASE( c04_rejects_invalid_inputs_with_context )
{
    BOOST_CHECK_THROW( parse( "" ), std::runtime_error );
    BOOST_CHECK_THROW( parse( "# header only\n" ), std::runtime_error );
    BOOST_CHECK_THROW( parse( "Date MJD x y UT1 LOD dPsi dEps\n" ), std::runtime_error );
    BOOST_CHECK_THROW( parse( "# DPSI DEPS\n" ), std::runtime_error );
    BOOST_CHECK_THROW( parse( "# dpsi deps\n" ), std::runtime_error );
    BOOST_CHECK_THROW( parse( "1962 1 1 37665 0 0 0\n" ), std::runtime_error );
    const std::string text = fixtureText( "legacy_c04.txt" );
    const auto first = text.find( "1962   1   1" );
    const auto last = text.find( '\n', first );
    const std::string row = text.substr( first, last - first );
    BOOST_CHECK_THROW( parse( row + "\n" + row + "\n" ), std::runtime_error );
    const auto nextEnd = text.find( '\n', last + 1 );
    const auto nextRow = text.substr( last + 1, nextEnd - last - 1 );
    BOOST_CHECK_THROW( parse( nextRow + "\n" + row + "\n" ), std::runtime_error );
    auto mismatch = row;
    mismatch.replace( mismatch.find( "37665" ), 5, "37664" );
    BOOST_CHECK_THROW( parse( mismatch ), std::runtime_error );
    auto nonfinite = row;
    nonfinite.replace( nonfinite.find( "-0.012700" ), 9, "nan" );
    BOOST_CHECK_THROW( parse( nonfinite ), std::runtime_error );
    auto invalidDate = row;
    invalidDate.replace( invalidDate.find( "1962   1   1" ), 12, "1962   2  30" );
    BOOST_CHECK_THROW( parse( invalidDate ), std::runtime_error );
    BOOST_CHECK_THROW( parse( row + " 0\n" ), std::runtime_error );
    try
    {
        parse( "# good header\n" + mismatch );
        BOOST_FAIL( "Calendar inconsistency was accepted" );
    }
    catch( const std::runtime_error& error )
    {
        BOOST_CHECK( std::string( error.what( ) ).find( "test source:2:" ) != std::string::npos );
    }
    BOOST_CHECK_THROW( readEarthOrientationData( "/nonexistent/eop/file" ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( finals_rejects_truncation_partial_coordinates_and_invalid_flags )
{
    const auto line = firstFinalsLine( );
    BOOST_CHECK_THROW( parse( line.substr( 0, 90 ), EarthOrientationFormat::finals2000a ), std::runtime_error );
    auto partial = line;
    partial.replace( 18, 9, 9, ' ' );
    BOOST_CHECK_THROW( parse( partial, EarthOrientationFormat::finals2000a ), std::runtime_error );
    auto partialCpo = line;
    partialCpo.replace( 97, 9, 9, ' ' );
    BOOST_CHECK_THROW( parse( partialCpo, EarthOrientationFormat::finals2000a ), std::runtime_error );
    auto wrongFlag = line;
    wrongFlag[ 16 ] = 'Q';
    BOOST_CHECK_THROW( parse( wrongFlag, EarthOrientationFormat::finals2000a ), std::runtime_error );
    auto nan = line;
    nan.replace( 58, 10, "       nan" );
    BOOST_CHECK_THROW( parse( nan, EarthOrientationFormat::finals2000a ), std::runtime_error );
    BOOST_CHECK_THROW( parse( line + "\n" + line, EarthOrientationFormat::finals2000a ), std::runtime_error );
    auto badDate = line;
    badDate.replace( 2, 2, "13" );
    BOOST_CHECK_THROW( parse( badDate, EarthOrientationFormat::finals2000a ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( validates_calendar_hours_headers_and_optional_uncertainties )
{
    const std::string leapRow = "2020 2 29 12 58908.50 .1 .2 -.1 .001 -.002 .003 -.004 .0005 .01 .02 .03 .04 .05 .06 .07 .08";
    const auto leap = parse( "# YR MM DD HH MJD x y UT1 dX dY\n" + leapRow );
    BOOST_REQUIRE_EQUAL( leap.records.size( ), 1 );
    BOOST_CHECK_EQUAL( leap.records.begin( )->first, 58908.5 );
    auto invalidHour = leapRow;
    invalidHour.replace( invalidHour.find( "12" ), 2, "24" );
    BOOST_CHECK_THROW( parse( invalidHour ), std::runtime_error );
    auto hourCountedTwice = leapRow;
    hourCountedTwice.replace( hourCountedTwice.find( "58908.50" ), 8, "58909.00" );
    BOOST_CHECK_THROW( parse( hourCountedTwice ), std::runtime_error );
    BOOST_CHECK_THROW( parse( "# YR MM DD HH MJD x y UT1 dX dY\n1962 1 1 37665 0 0 0 0 0 0 0 0 0 0 0 0" ), std::runtime_error );
    auto negativeError = leapRow;
    negativeError += " ";
    negativeError.replace( negativeError.rfind( ".08" ), 3, "-.08" );
    BOOST_CHECK_THROW( parse( negativeError ), std::runtime_error );
    auto errorWithoutCpo = firstFinalsLine( );
    errorWithoutCpo.replace( 95, 1, " " );
    errorWithoutCpo.replace( 97, 9, 9, ' ' );
    errorWithoutCpo.replace( 116, 9, 9, ' ' );
    BOOST_CHECK_THROW( parse( errorWithoutCpo, EarthOrientationFormat::finals2000a ), std::runtime_error );
    auto negativeFinalsError = firstFinalsLine( );
    negativeFinalsError.replace( 27, 9, "-0.009786" );
    BOOST_CHECK_THROW( parse( negativeFinalsError, EarthOrientationFormat::finals2000a ), std::runtime_error );
    auto partialUncertainty = firstFinalsLine( );
    partialUncertainty.replace( 46, 9, 9, ' ' );
    BOOST_CHECK_THROW( parse( partialUncertainty, EarthOrientationFormat::finals2000a ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( validates_public_preloaded_data_contract )
{
    auto data = parse( fixtureText( "modern_c04.txt" ) );
    BOOST_CHECK_NO_THROW( validateEarthOrientationData( data ) );
    auto invalid = data;
    invalid.records.begin( )->second.polarMotion[ 0 ] = std::numeric_limits< double >::infinity( );
    BOOST_CHECK_THROW( validateEarthOrientationData( invalid ), std::runtime_error );
    invalid = data;
    invalid.records.begin( )->second.utcModifiedJulianDate += 1.0;
    BOOST_CHECK_THROW( validateEarthOrientationData( invalid ), std::runtime_error );
    invalid = data;
    invalid.records.begin( )->second.polarMotionRateUncertainty = std::array< double, 2 >{ 0.0, -1.0 };
    BOOST_CHECK_THROW( validateEarthOrientationData( invalid ), std::runtime_error );
    invalid = data;
    invalid.records.begin( )->second.celestialPoleOffsets = std::array< double, 2 >{ 0.0, std::numeric_limits< double >::quiet_NaN( ) };
    BOOST_CHECK_THROW( validateEarthOrientationData( invalid ), std::runtime_error );
    invalid = data;
    invalid.records.begin( )->second.lengthOfDay.reset( );
    BOOST_CHECK_THROW( validateEarthOrientationData( invalid ), std::runtime_error );
}
