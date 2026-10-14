#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <sstream>
#include "tudat/io/readAtmosphericTideData.h"
using namespace tudat::input_output;
namespace
{
std::string path( const std::string& name )
{
    return std::string( __FILE__ ).substr( 0, std::string( __FILE__ ).find_last_of( "/\\" ) ) +
            "/../../reference_files/atmospheric_tides/" + name;
}
std::string text( const std::string& name )
{
    std::ifstream input( path( name ) );
    BOOST_REQUIRE( input.good( ) );
    return std::string( std::istreambuf_iterator< char >( input ), {} );
}
std::string replace( std::string original, const std::string& from, const std::string& to )
{
    const auto i = original.find( from );
    BOOST_REQUIRE( i != std::string::npos );
    original.replace( i, from.size( ), to );
    return original;
}
void rejectGrid( const std::string& value )
{
    std::istringstream input( value );
    BOOST_CHECK_THROW( readAtmosphericTideGrid( input ), std::exception );
}
void rejectStation( const std::string& value )
{
    std::istringstream input( value );
    BOOST_CHECK_THROW( readAtmosphericTideStationData( input ), std::exception );
}
}  // namespace
BOOST_AUTO_TEST_CASE( realProviderGridNormalizesUnitsAndRetainsMetadata )
{
    AtmosphericTideMetadata metadata;
    metadata.referenceOrigin = AtmosphericLoadingReferenceOrigin::solidEarthCentre;
    const auto grid = readAtmosphericTideGrid( path( "rp03_ce_grid_extract.txt" ), metadata );
    BOOST_CHECK_EQUAL( grid.coefficientsMetres.size( ), 4 );
    BOOST_CHECK_EQUAL( grid.metadata.modelName, "Ray-Ponte 2003" );
    BOOST_CHECK_EQUAL( grid.metadata.originalCoefficientUnits, "mm" );
    BOOST_CHECK_EQUAL( grid.metadata.sourceHeaders.size( ), 3 );
    BOOST_CHECK( grid.metadata.referenceOrigin == AtmosphericLoadingReferenceOrigin::solidEarthCentre );
    BOOST_CHECK( grid.metadata.latitudeConvention == AtmosphericTideLatitudeConvention::unspecified );
    BOOST_CHECK_EQUAL( grid.longitudesDegrees.front( ), 47 );
    BOOST_CHECK_EQUAL( grid.latitudesDegrees.front( ), -19 );
    BOOST_CHECK_EQUAL( grid.latitudesDegrees.back( ), -20 );
    std::istringstream source( text( "rp03_ce_grid_extract.txt" ) );
    BOOST_CHECK( readAtmosphericTideGrid( source ).metadata.referenceOrigin == AtmosphericLoadingReferenceOrigin::unspecified );
    // Frozen provider literals in mm, independently converted to SI.
    BOOST_CHECK_SMALL( grid.coefficientsMetres[ 0 ]( 0, 0 ) + .0002344, 1e-18 );
    BOOST_CHECK_SMALL( grid.coefficientsMetres[ 0 ]( 1, 0 ) - .0000595, 1e-18 );
    BOOST_CHECK_SMALL( grid.coefficientsMetres[ 0 ]( 2, 3 ) - .0000597, 1e-18 );
    const auto originalHeader = text( "rp03_ce_grid_extract.txt" );
    std::istringstream indented( "  " + originalHeader );
    BOOST_CHECK_EQUAL( readAtmosphericTideGrid( indented ).metadata.sourceHeaders.front( ), "  " + grid.metadata.sourceHeaders.front( ) );
    std::string windows;
    for( char c : text( "rp03_ce_grid_extract.txt" ) )
    {
        if( c == '\n' )
        {
            windows += '\r';
        }
        windows += c;
    }
    std::istringstream crlf( windows );
    BOOST_CHECK_EQUAL( readAtmosphericTideGrid( crlf ).coefficientsMetres.size( ), 4 );
}
BOOST_AUTO_TEST_CASE( publishedManualGridNodeMatchesPinnedProvider )
{
    // Bernese manual Figure24.17 prints these exact radial/north values at lon0,lat90.
    // This published grid-node check is distinct from its incompatible station examples in Figure24.48.
    const auto grid = readAtmosphericTideGrid( path( "rp03_ce_polar_node_extract.txt" ) );
    const auto& coefficients = grid.coefficientsMetres.front( );
    BOOST_CHECK_SMALL( coefficients( 0, 0 ) - .0000183, 1e-18 );
    BOOST_CHECK_SMALL( coefficients( 0, 1 ) - .0000761, 1e-18 );
    BOOST_CHECK_SMALL( coefficients( 0, 2 ) + .0000787, 1e-18 );
    BOOST_CHECK_SMALL( coefficients( 0, 3 ) - .0000139, 1e-18 );
    BOOST_CHECK_SMALL( coefficients( 1, 0 ) - .0000077, 1e-18 );
}
BOOST_AUTO_TEST_CASE( stationIdentityCoordinatesAndSignedSiCoefficients )
{
    const auto data = readAtmosphericTideStationData( path( "bernese_manual_stations.atl" ) );
    BOOST_CHECK_EQUAL( data.stations.size( ), 2 );
    const auto& a = data.getStation( "ABPO" );
    BOOST_CHECK_EQUAL( a.providerStationIdentifierLine, "ABPO 33302M001" );
    BOOST_CHECK( a.hasProviderCoordinates );
    BOOST_CHECK_EQUAL( a.providerLongitudeLatitudeDegrees.x( ), 47.22921 );
    BOOST_CHECK_EQUAL( a.providerLongitudeLatitudeDegrees.y( ), -19.018307 );
    BOOST_CHECK_SMALL( a.coefficientsMetres( 0, 0 ) + .0003414, 1e-18 );
    BOOST_CHECK_SMALL( a.coefficientsMetres( 1, 2 ) - .0000310, 1e-18 );
    BOOST_CHECK_SMALL( data.getStation( "ADE1" ).coefficientsMetres( 2, 3 ) - .0001892, 1e-18 );
    BOOST_CHECK_THROW( data.getStation( "MISSING" ), std::runtime_error );
}
BOOST_AUTO_TEST_CASE( invalidGridAndStationStructures )
{
    const auto grid = text( "rp03_ce_grid_extract.txt" );
    rejectGrid( "" );
    rejectGrid( grid.substr( 0, grid.find_last_of( '\n', grid.size( ) - 2 ) ) );
    rejectGrid( grid + grid.substr( grid.find_last_of( '\n', grid.size( ) - 2 ) + 1 ) );
    rejectGrid( replace( grid, "47.000  -20.000", "47.000  -19.000" ) );
    rejectGrid( replace( grid, "47.000  -20.000", "47.000  -91.000" ) );
    rejectGrid( "0 0 nan 0 0 0 0 0 0 0 0 0 0 0\n" );
    rejectGrid( "$$ units metres\n" + grid );
    rejectGrid( "$$ Column order sins1 coss1 coss2 sins2\n" + grid );
    const auto station = text( "bernese_manual_stations.atl" );
    rejectStation( station.substr( 0, station.find_last_of( '\n', station.size( ) - 2 ) ) );
    rejectStation( station + "ABPO duplicate\n0 0 0 0\n0 0 0 0\n0 0 0 0\n" );
    rejectStation( replace( station, "-.3414", "nan" ) );
    rejectStation( replace( station, "-19.018307", "-91" ) );
    rejectStation( replace( station, "$$ TANGENTIAL NS\n$$ TANGENTIAL EW", "$$ TANGENTIAL EW\n$$ TANGENTIAL NS" ) );
    rejectStation( replace( station, "positive up, north and east", "positive up, north and west" ) );
    BOOST_CHECK_THROW( readAtmosphericTideGrid( "/missing/atmospheric/grid" ), std::runtime_error );
    AtmosphericTideMetadata expectedCounts;
    expectedCounts.expectedLongitudeCount = 361;
    std::istringstream incomplete( grid );
    BOOST_CHECK_THROW( readAtmosphericTideGrid( incomplete, expectedCounts ), std::runtime_error );
    AtmosphericTideMetadata wrong;
    wrong.modelName = "different model";
    std::istringstream declared( grid );
    BOOST_CHECK_THROW( readAtmosphericTideGrid( declared, wrong ), std::runtime_error );
}
BOOST_AUTO_TEST_CASE( diagnosticsAndStreamFailures )
{
    AtmosphericTideMetadata metadata;
    metadata.source = "provider-grid.dat";
    std::istringstream invalid( "$$ source\n0 0 broken 0 0 0 0 0 0 0 0 0 0 0\n" );
    try
    {
        readAtmosphericTideGrid( invalid, metadata );
        BOOST_FAIL( "Malformed number accepted" );
    }
    catch( const std::runtime_error& e )
    {
        BOOST_CHECK_NE( std::string( e.what( ) ).find( "provider-grid.dat:2:" ), std::string::npos );
    }
    std::istringstream failed( text( "rp03_ce_grid_extract.txt" ) );
    failed.setstate( std::ios::failbit );
    BOOST_CHECK_THROW( readAtmosphericTideGrid( failed ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( numericIdentifiersAndCoordinatesWithinCurrentBlock )
{
    std::istringstream input(
            "7839\n$$ station7839; coord.(long,lat)47 -19\n1 2 3 4\n5 6 7 8\n9 10 11 12\n7840\n0 0 0 0\n0 0 0 0\n0 0 0 0\n" );
    const auto data = readAtmosphericTideStationData( input );
    BOOST_CHECK_EQUAL( data.stations.size( ), 2 );
    BOOST_CHECK( data.getStation( "7839" ).hasProviderCoordinates );
    BOOST_CHECK_EQUAL( data.getStation( "7839" ).providerLongitudeLatitudeDegrees.y( ), -19 );
    BOOST_CHECK( !data.getStation( "7840" ).hasProviderCoordinates );
    BOOST_CHECK_SMALL( data.getStation( "7839" ).coefficientsMetres( 2, 3 ) - .012, 1e-18 );
    rejectStation( "1 2 3 4\n5 6 7 8\n9 10 11 12\n" );
    rejectStation( "7839\n$$ coord.(long,lat)47 -19\n$$ coord.(long,lat)48 -19\n0 0 0 0\n0 0 0 0\n0 0 0 0\n" );
    AtmosphericTideMetadata conflicting;
    conflicting.greensFunction = "another Green function";
    std::istringstream grid( text( "rp03_ce_grid_extract.txt" ) );
    BOOST_CHECK_THROW( readAtmosphericTideGrid( grid, conflicting ), std::runtime_error );
}
