#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/readIcgemGravityField.h"
#include <sstream>
#include <limits>
#include <algorithm>

using namespace tudat::input_output;
namespace
{
//! Locate the committed coefficient extracts independently of the executable's working directory.
std::string fixture( const std::string& name )
{
    const std::string source = __FILE__;
    return source.substr( 0, source.find_last_of( "/\\" ) ) + "/../../reference_files/icgem_gravity/" + name;
}

//! Supply a minimal controlled header with a selectable uncertainty-column declaration.
std::string header( const std::string& errors = "no" )
{
    return "product_type gravity_field\nmodelname controlled-static\nearth_gravity_constant 3.986004415D14\n"
           "radius 6.37813646D6\nmax_degree 2\nerrors " +
            errors + "\nend_of_head trailing comment\n";
}

//! Exercise the same stream parser used for provider files.
IcgemGravityFieldData parse( const std::string& text )
{
    std::istringstream stream( text );
    return readIcgemGravityField( stream, "controlled.gfc" );
}
}  // namespace

BOOST_AUTO_TEST_CASE( publishedGfzHeaderPreambleAndStaticSnapshotArePreserved )
{
    const auto raw = readIcgemGravityField( fixture( "gfz_rl061_ddk4_202104_degree2.gfc" ) );
    BOOST_CHECK_EQUAL( raw.maximumDegree, 60 );
    BOOST_CHECK_EQUAL( raw.modelName, "GSM-2_2021091-2021120_GRFO_GFZOP_BA01_0601" );
    BOOST_CHECK_EQUAL( raw.headerLines[ 15 ], "max_degree         60" );  // Comment section precedes begin_of_head.
    BOOST_CHECK( raw.normalization == IcgemNormalization::fullyNormalized );
    BOOST_CHECK( raw.tideSystem == tudat::basic_astrodynamics::TideSystem::tideFree );
    BOOST_CHECK_EQUAL( raw.coefficients.size( ), 6 );
    BOOST_CHECK_EQUAL( raw.coefficients.at( { 2, 0 } ).cosine, -4.841653458641e-4 );
    BOOST_CHECK_EQUAL( ( *raw.coefficients.at( { 2, 2 } ).formalErrors )[ 0 ], 1.0130e-12 );
    BOOST_CHECK_EQUAL( raw.coefficients.at( { 2, 0 } ).sourceLine, 59 );
}

BOOST_AUTO_TEST_CASE( rawPresenceProviderDefaultsCommentsAndDeclaredErrorColumns )
{
    const auto raw = parse( "modelname this_is_a_pre_header_comment\nbegin_of_head =====\n" + header( "calibrated_and_formal" ) +
                            "comment data section\ngfc 0 0 1D0 0 1D-8 2D-8 3D-8 4D-8 numeric 99 trailing\n"
                            "gfc 2 0 -4.84D-4 0 5D-10 6D-10 7D-10 8D-10\n" );
    BOOST_CHECK_EQUAL( raw.modelName, "controlled-static" );
    BOOST_CHECK_EQUAL( raw.gravitationalParameter, 3.986004415e14 );
    BOOST_CHECK_EQUAL( raw.referenceRadius, 6.37813646e6 );
    BOOST_CHECK_EQUAL( raw.maximumDegree, 2 );
    BOOST_CHECK( raw.normalization == IcgemNormalization::fullyNormalized );
    BOOST_CHECK( raw.tideSystem == tudat::basic_astrodynamics::TideSystem::unknown );
    BOOST_CHECK_EQUAL( raw.coefficients.size( ), 2 );
    BOOST_CHECK_EQUAL( raw.coefficients.count( { 1, 0 } ), 0 );
    BOOST_CHECK_EQUAL( raw.coefficients.at( { 0, 0 } ).sine, 0.0 );
    BOOST_CHECK_EQUAL( ( *raw.coefficients.at( { 0, 0 } ).calibratedErrors )[ 1 ], 2e-8 );
    BOOST_CHECK_EQUAL( ( *raw.coefficients.at( { 0, 0 } ).formalErrors )[ 0 ], 3e-8 );
    BOOST_CHECK_EQUAL( raw.headerLines.front( ), "modelname this_is_a_pre_header_comment" );
    BOOST_CHECK( raw.coefficients.at( { 0, 0 } ).originalRecord.find( "numeric 99 trailing" ) != std::string::npos );
}

BOOST_AUTO_TEST_CASE( readerRejectsConflictingPhysicalRecords )
{
    BOOST_CHECK_THROW( parse( header( ) + "gfc 2 3 1 0\n" ), std::runtime_error );
    BOOST_CHECK_THROW( parse( header( ) + "gfc 2 0 nan 0\n" ), std::runtime_error );
    BOOST_CHECK_THROW( parse( header( "formal" ) + "gfc 2 0 1 0 0.1\n" ), std::runtime_error );
    auto raw = parse( header( ) + "gfc 0 0 1 0\n" );
    raw.coefficients.at( { 0, 0 } ).cosine = std::numeric_limits< double >::infinity( );
    BOOST_CHECK_THROW( validateIcgemGravityFieldData( raw ), std::runtime_error );
}

//! Preserve independent provider intervals, including subday tokens that must not be silently rolled over.
BOOST_AUTO_TEST_CASE( publishedTemporalC20KeepsIndependentIntervalsAndSparsePresence )
{
    const auto data = readIcgemGravityField( fixture( "eigen6s4_v2_c20_1950_2006.gfc" ) );

    BOOST_CHECK_EQUAL( data.modelName, "EIGEN-6S4v2" );
    BOOST_CHECK_EQUAL( data.maximumDegree, 300 );
    BOOST_CHECK_EQUAL( data.format, "icgem2.0" );
    BOOST_CHECK_EQUAL( data.coefficients.size( ), 1 );
    BOOST_CHECK_EQUAL( data.coefficients.count( { 2, 0 } ), 0 );
    BOOST_REQUIRE_EQUAL( data.temporalCoefficients.size( ), 60 );

    const auto annual = std::find_if( data.temporalCoefficients.begin( ), data.temporalCoefficients.end( ), []( const auto& record ) {
        return record.type == "acos" && record.referenceDate == "19500101.0000" && record.periodYears == 1.0;
    } );

    BOOST_REQUIRE( annual != data.temporalCoefficients.end( ) );
    BOOST_REQUIRE( annual->coefficient.calibratedErrors );
    BOOST_CHECK_EQUAL( annual->validityEndDate, "20030101.0000" );
    BOOST_CHECK_EQUAL( annual->coefficient.cosine, 3.87007395388e-11 );
    BOOST_CHECK_EQUAL( ( *annual->coefficient.calibratedErrors )[ 0 ], 1.1170e-12 );
    BOOST_CHECK_EQUAL( annual->coefficient.sourceLine, 110 );
    BOOST_CHECK( std::any_of( data.temporalCoefficients.begin( ), data.temporalCoefficients.end( ), []( const auto& record ) {
        return record.referenceDate == "20041226.0060";
    } ) );

    // The long annual interval requires the union of adjacent reference pieces, including this removed one.
    auto incomplete = data;

    incomplete.temporalCoefficients.erase(
            std::remove_if( incomplete.temporalCoefficients.begin( ),
                            incomplete.temporalCoefficients.end( ),
                            []( const auto& record ) { return record.type == "gfct" && record.referenceDate == "19850109.1751"; } ),
            incomplete.temporalCoefficients.end( ) );
    BOOST_CHECK_THROW( validateIcgemGravityFieldData( incomplete ), std::runtime_error );
}

//! Associate 1.0 terms with gfct without manufacturing printed dates; retain the old dot-only convention.
BOOST_AUTO_TEST_CASE( temporalV1FieldsAndOldTrendRestrictions )
{
    const auto modern = parse( header( "calibrated_and_formal" ) +
                               "gfct 2 0 -4.84D-4 0 1D-12 0 2D-12 0 20050101\n"
                               "trnd 2 0 3D-12 0 4D-13 0 5D-13 0\n"
                               "acos 2 0 6D-12 0 7D-13 0 8D-13 0 1.0 trailing comment\n"
                               "asin 2 0 9D-12 0 1D-13 0 2D-13 0 0.5\n" );

    BOOST_REQUIRE_EQUAL( modern.temporalCoefficients.size( ), 4 );
    BOOST_CHECK_EQUAL( modern.temporalCoefficients.front( ).referenceDate, "20050101" );
    BOOST_CHECK( modern.temporalCoefficients.at( 2 ).referenceDate.empty( ) );
    BOOST_REQUIRE( modern.temporalCoefficients.at( 2 ).periodYears );
    BOOST_REQUIRE( modern.temporalCoefficients.at( 2 ).coefficient.formalErrors );
    BOOST_CHECK_EQUAL( *modern.temporalCoefficients.at( 2 ).periodYears, 1.0 );
    BOOST_CHECK_EQUAL( ( *modern.temporalCoefficients.at( 2 ).coefficient.formalErrors )[ 0 ], 8e-13 );

    const auto old = parse( header( ) + "gfct 2 0 -4.84e-4 0 20050101\ndot 2 0 1e-12 0\n" );

    BOOST_CHECK_EQUAL( old.temporalCoefficients.at( 1 ).type, "dot" );
    BOOST_CHECK_THROW( parse( header( ) + "gfct 2 0 -4.84e-4 0 20050101\ndot 2 0 1e-12 0\nasin 2 0 1e-12 0 1.0\n" ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( publishedCostgExtractRetainsDuplicateC00AsARejection )
{
    const auto file = fixture( "costg_eigen6_c4_duplicate_c00.gfc" );
    try
    {
        readIcgemGravityField( file );
        BOOST_FAIL( "Published duplicate C00 must reject" );
    }
    catch( const std::runtime_error& error )
    {
        BOOST_CHECK( std::string( error.what( ) ).find( file + ":16: Duplicate gfc degree/order" ) != std::string::npos );
    }
}
