#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "tudat/io/readOceanLoadingData.h"

namespace tudat
{
namespace unit_tests
{
namespace
{

// Provider-format coefficients from the corrected IERS HARDISP Onsala example.
// Global provider header syntax is retained, including variable comments and CRLF.
const std::string providerBlock =
        "$$ Ocean loading displacement\r\n"
        "$$ COLUMN ORDER: M2 S2 N2 K2 K1 O1 P1 Q1 MF MM SSA\r\n"
        "$$ AMPLITUDES (m)\r\n$$ PHASES (degrees)\r\n"
        "$$ Ocean tide model: TPXO.7.2, Ssa tide from TPXO.6.2\r\n"
        "$$ Gutenberg-Bullen Greens function is used.\r\n"
        "$$ CMC: NO\r\n$$ END HEADER\r\n\r\n"
        "  ONSALA\r\n"
        "$$ CSR4.0_f_PP ID: 2009-06-25 17:43:24\r\n"
        "$$ Computed by OLMPP by H G Scherneck, Onsala Space Observatory, 2009\r\n"
        "$$ An extra provider comment\r\n"
        "$$ Onsala, lon/lat: 11.9264 57.3958 0.00\r\n"
        " .00352 .00123 .00080 .00032 .00187 .00112 .00063 .00003 .00082 .00044 .00037\r\n"
        " .00144 .00035 .00035 .00008 .00053 .00049 .00018 .00009 .00012 .00005 .00006\r\n"
        " .00086 .00023 .00023 .00006 .00029 .00028 .00010 .00007 .00004 .00002 .00001\r\n"
        " -64.7 -52.0 -96.2 -55.2 -58.8 -151.4 -65.6 -138.1 8.4 5.2 2.1\r\n"
        " 85.5 114.5 56.5 113.6 99.4 19.1 94.1 -10.4 -167.4 -170.0 -177.7\r\n"
        " 109.5 147.0 92.7 148.8 50.5 -55.1 36.4 -170.4 -15.0 2.3 5.2\r\n$$\r\n";

std::string replaceText( std::string text, const std::string& oldText, const std::string& newText )
{
    const auto position = text.find( oldText );
    BOOST_REQUIRE( position != std::string::npos );
    text.replace( position, oldText.size( ), newText );
    return text;
}

input_output::OceanLoadingData parseText( const std::string& text )
{
    std::istringstream stream( text );
    return input_output::readOceanLoadingData( stream, "provider.blq" );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_ocean_loading_data )

BOOST_AUTO_TEST_CASE( testProviderHeadersAndCoefficientConventions )
{
    const auto data = parseText( providerBlock );
    BOOST_CHECK_EQUAL( data.stations.size( ), 1 );
    const auto& station = data.getStation( "ONSALA" );
    BOOST_CHECK_EQUAL( station.stationIdentifier, "ONSALA" );
    BOOST_CHECK_EQUAL( station.source, "provider.blq" );
    BOOST_CHECK_EQUAL( data.sourceHeaders.size( ), 8 );
    BOOST_CHECK_EQUAL( station.sourceHeaders.size( ), 5 );
    BOOST_CHECK_EQUAL( station.oceanTideModel, "CSR4.0_f_PP" );
    BOOST_CHECK_EQUAL( station.greensFunction, "Gutenberg-Bullen Greens function is used." );
    BOOST_CHECK( station.referenceOrigin == input_output::OceanLoadingReferenceOrigin::solidEarthCentre );
    BOOST_CHECK( station.hasProviderCoordinates );
    BOOST_CHECK_EQUAL( station.providerLongitudeLatitudeHeight( 0 ), 11.9264 );
    BOOST_CHECK_EQUAL( station.providerLongitudeLatitudeHeight( 1 ), 57.3958 );
    BOOST_CHECK_EQUAL( station.amplitudesMetres( 0, 0 ), 0.00352 );
    BOOST_CHECK_EQUAL( station.amplitudesMetres( 1, 1 ), 0.00035 );
    BOOST_CHECK_EQUAL( station.amplitudesMetres( 2, 10 ), 0.00001 );
    BOOST_CHECK_EQUAL( station.phaseLagsDegrees( 0, 0 ), -64.7 );
    BOOST_CHECK_EQUAL( station.phaseLagsDegrees( 1, 10 ), -177.7 );
    BOOST_CHECK_EQUAL( station.phaseLagsDegrees( 2, 10 ), 5.2 );
    BOOST_CHECK_EQUAL( input_output::getOceanLoadingConstituentNames( )[ 8 ], "Mf" );
    BOOST_CHECK_THROW( data.getStation( "MISSING" ), std::runtime_error );
    BOOST_CHECK_THROW( data.getStation( "onsala" ), std::runtime_error );

    const auto cm = parseText( replaceText( providerBlock, "CMC: NO", "CMC: YES" ) );
    BOOST_CHECK( cm.getStation( "ONSALA" ).referenceOrigin == input_output::OceanLoadingReferenceOrigin::wholeEarthCentreOfMass );
    const auto unspecified = parseText( replaceText( providerBlock, "$$ CMC: NO\r\n", "" ) );
    BOOST_CHECK( unspecified.getStation( "ONSALA" ).referenceOrigin == input_output::OceanLoadingReferenceOrigin::unspecified );
    const auto scientific = parseText( replaceText( providerBlock, ".00352", "3.52D-3" ) );
    BOOST_CHECK_EQUAL( scientific.getStation( "ONSALA" ).amplitudesMetres( 0, 0 ), 0.00352 );
}

BOOST_AUTO_TEST_CASE( testMultipleStationsAndOptionalHeaders )
{
    const std::string dataOnly = providerBlock.substr( providerBlock.find( " .00352" ) );
    const auto data = parseText( providerBlock + "\nSECOND\n$$ Complete FES2022b\n$$ FRAME: CM\n" + dataOnly );
    BOOST_CHECK_EQUAL( data.stations.size( ), 2 );
    BOOST_CHECK_EQUAL( data.getStation( "SECOND" ).oceanTideModel, "FES2022b" );
    BOOST_CHECK( data.getStation( "SECOND" ).referenceOrigin == input_output::OceanLoadingReferenceOrigin::wholeEarthCentreOfMass );
    const auto noComments = parseText( "NUMERIC_ID_123\n" + dataOnly );
    BOOST_CHECK_EQUAL( noComments.stations.size( ), 1 );
    BOOST_CHECK( noComments.getStation( "NUMERIC_ID_123" ).referenceOrigin == input_output::OceanLoadingReferenceOrigin::unspecified );
    const auto numericId = parseText( "12345\n" + dataOnly );
    BOOST_CHECK_EQUAL( numericId.getStation( "12345" ).stationIdentifier, "12345" );
}

BOOST_AUTO_TEST_CASE( testRejectMalformedBlqData )
{
    const std::string dataOnly = providerBlock.substr( providerBlock.find( " .00352" ) );
    BOOST_CHECK_THROW( parseText( "" ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( "$$ header only\n" ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( dataOnly ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( providerBlock + "ONSALA\n" + dataOnly ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( providerBlock.substr( 0, providerBlock.find( " 109.5" ) ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, ".00352 .00123", ".00352" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, ".00352", ".00352 0.0" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, ".00352", "garbage" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, ".00352", "-.00352" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, ".00352", "nan" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, "-64.7", "inf" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, "M2 S2", "S2 M2" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, "AMPLITUDES (m)", "AMPLITUDES (mm)" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, "PHASES (degrees)", "PHASES (radians)" ) ), std::runtime_error );
    BOOST_CHECK_THROW( input_output::readOceanLoadingData( "__missing_ocean_loading_fixture__.blq" ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( testDeclaredRowOrderAndProviderLatitude )
{
    const std::string rows =
            "$$ ROW ORDER:\n$$ AMPLITUDES (m)\n$$ RADIAL\n$$ TANGENTL EW\n$$ TANGENTL NS\n"
            "$$ PHASES (degrees)\n$$ RADIAL\n$$ TANGENTL EW\n$$ TANGENTL NS\n";
    BOOST_CHECK_EQUAL( parseText( rows + providerBlock ).stations.size( ), 1 );
    BOOST_CHECK_EQUAL( parseText( "$$ ROW ORDER: UP WEST SOUTH\n" + providerBlock ).stations.size( ), 1 );
    BOOST_CHECK_THROW( parseText( replaceText( rows, "TANGENTL EW\n$$ TANGENTL NS", "TANGENTL NS\n$$ TANGENTL EW" ) + providerBlock ),
                       std::runtime_error );
    BOOST_CHECK_THROW( parseText( "$$ ROW ORDER: UP EAST NORTH\n" + providerBlock ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( rows, "$$ RADIAL\n", "$$ NORTH\n" ) + providerBlock ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( rows, "AMPLITUDES (m)", "PHASES (degrees)" ) + providerBlock ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, "57.3958", "90.0001" ) ), std::runtime_error );
    BOOST_CHECK_THROW( parseText( replaceText( providerBlock, "57.3958", "-90.0001" ) ), std::runtime_error );
    auto data = parseText( providerBlock );
    data.stations.at( "ONSALA" ).providerLongitudeLatitudeHeight( 1 ) = 100.0;
    BOOST_CHECK_THROW( input_output::validateOceanLoadingData( data ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( testValidatePreloadedOceanLoadingData )
{
    auto data = parseText( providerBlock );
    data.stations.at( "ONSALA" ).phaseLagsDegrees( 1, 1 ) = std::numeric_limits< double >::quiet_NaN( );
    BOOST_CHECK_THROW( input_output::validateOceanLoadingData( data ), std::runtime_error );
    data = parseText( providerBlock );
    data.stations.at( "ONSALA" ).stationIdentifier = "DIFFERENT";
    BOOST_CHECK_THROW( input_output::validateOceanLoadingData( data ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
