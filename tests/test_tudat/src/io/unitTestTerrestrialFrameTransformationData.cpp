#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <sstream>
#include "tudat/io/readTerrestrialFrameTransformationData.h"
using namespace tudat::input_output;
namespace
{
std::string text( const std::string& name )
{
    const auto path = std::string( __FILE__ ).substr( 0, std::string( __FILE__ ).find_last_of( "/\\" ) ) +
            "/../../reference_files/itrf_transformations/" + name;
    std::ifstream file( path );
    BOOST_REQUIRE( file.good( ) );
    return std::string( std::istreambuf_iterator< char >( file ), {} );
}
std::string replace( std::string original, const std::string& from, const std::string& to )
{
    const auto i = original.find( from );
    BOOST_REQUIRE( i != std::string::npos );
    original.replace( i, from.size( ), to );
    return original;
}
void reject( const std::string& value )
{
    std::istringstream stream( value );
    BOOST_CHECK_THROW( readTerrestrialFrameTransformationData( stream ), std::runtime_error );
}
}  // namespace
BOOST_AUTO_TEST_CASE( officialTablesNormalizeUnitsAndRetainIdentityProvenance )
{
    const auto legacy = readOfficialItrf2014TransformationData( ), current = readOfficialItrf2020TransformationData( );
    BOOST_CHECK_EQUAL( legacy.records.size( ), 12 );
    BOOST_CHECK_EQUAL( current.records.size( ), 13 );
    BOOST_CHECK_EQUAL( legacy.referenceFrame.realization, "ITRF2014" );
    BOOST_CHECK_EQUAL( current.referenceFrame.realization, "ITRF2020" );
    BOOST_CHECK_EQUAL( current.provenance.sourceSha256, "0d8dc872e7ef994966063e614adb5b1b0fe20582359c6e77f1d83259d1ff1a68" );
    const auto& p = legacy.getParameters( "ITRF90" );
    BOOST_CHECK_EQUAL( p.referenceDecimalYear, 2010 );
    BOOST_CHECK_SMALL( p.translationMetres.y( ) - .0115, 1e-18 );
    BOOST_CHECK_SMALL( p.translationRateMetresPerSecond.y( ) * 31557600 + .0005, 1e-18 );
    BOOST_CHECK( p.rotationConvention == HelmertRotationConvention::positionVector );
    BOOST_CHECK_EQUAL( p.originalRotationUnits, "mas" );
    BOOST_CHECK_EQUAL( p.sourceFrame.release, "ITRF2014" );
    BOOST_CHECK( !p.originalParameterLine.empty( ) );
    BOOST_CHECK( p.rateSourceLine == p.parameterSourceLine + 1 );
    std::istringstream stream( text( "Transfo-ITRF2020_TRFs.txt" ) );
    const auto custom = readTerrestrialFrameTransformationData( stream, { "test-provider", "catalogue-release", "declared-source-hash" } );
    BOOST_CHECK_EQUAL( custom.referenceFrame.release, "catalogue-release" );
    BOOST_CHECK_EQUAL( custom.getParameters( "ITRF2014" ).provenance.source, "test-provider" );
    BOOST_CHECK_EQUAL( custom.getParameters( "ITRF2014" ).sourceFrame.release, "catalogue-release" );
    BOOST_CHECK_THROW( custom.getParameters( "unknown" ), std::runtime_error );
}
BOOST_AUTO_TEST_CASE( strictHeaderUnitsPairingAndFiniteChecks )
{
    const auto original = text( "Transfo-ITRF2014_ITRFs.txt" );
    reject( "" );
    reject( original.substr( 0, original.find( "________________________________________________________________" ) ) );
    reject( replace( original, "Tx       Ty       Tz", "Ty       Tx       Tz" ) );
    reject( replace( original, "mm       mm       mm", "cm       mm       mm" ) );
    reject( replace( original, "mm/y     mm/y", "mm/s     mm/y" ) );
    reject( replace( original, "  ITRF2008        1.6", "  ITRF2008        nan" ) );
    reject( replace( original, "       rates      0.0      0.0     -0.1      0.03      0.00      0.00      0.00\n", "" ) );
    reject( replace( original, "ITRF2005", "ITRF2008" ) );
    reject( replace( original, "-Rz   Ry", "Rz   Ry" ) );
    reject( replace( original, "0.26    2010.0", "0.26    2010.0junk" ) );
    std::istringstream dExponent( replace( original, "        1.6", "        1.6D0" ) );
    BOOST_CHECK_EQUAL( readTerrestrialFrameTransformationData( dExponent ).records.size( ), 12 );
    std::istringstream badStream( original );
    badStream.setstate( std::ios::badbit );
    BOOST_CHECK_THROW( readTerrestrialFrameTransformationData( badStream ), std::runtime_error );
    BOOST_CHECK_THROW( readTerrestrialFrameTransformationData( "/missing/itrf/table" ), std::runtime_error );
}
BOOST_AUTO_TEST_CASE( diagnosticsAndMutableRecordValidation )
{
    std::istringstream broken( replace( text( "Transfo-ITRF2014_ITRFs.txt" ), "  ITRF2008        1.6", "  ITRF2008        1.6e999" ) );
    try
    {
        readTerrestrialFrameTransformationData( broken, { "provider-table", "", "" } );
        BOOST_FAIL( "Expected finite parser rejection" );
    }
    catch( const std::runtime_error& e )
    {
        BOOST_CHECK( std::string( e.what( ) ).find( "provider-table:9:" ) != std::string::npos );
    }
    auto data = readOfficialItrf2014TransformationData( );
    data.records.begin( )->second.sourceFrame.realization = "ITRF2020";
    BOOST_CHECK_THROW( validateTerrestrialFrameTransformationData( data ), std::runtime_error );
    data = readOfficialItrf2014TransformationData( );
    data.records.begin( )->second.provenance.sourceRelease = "another-release";
    BOOST_CHECK_THROW( validateTerrestrialFrameTransformationData( data ), std::runtime_error );
    data = readOfficialItrf2014TransformationData( );
    data.records.begin( )->second.translationMetres.x( ) = NAN;
    BOOST_CHECK_THROW( validateTerrestrialFrameTransformationData( data ), std::runtime_error );
}
