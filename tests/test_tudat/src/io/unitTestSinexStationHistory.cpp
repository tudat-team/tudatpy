#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/readSinexStationHistory.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

using namespace tudat::input_output;
namespace
{
std::string fixture( const std::string& name )
{
    return std::string( __FILE__ ).substr( 0, std::string( __FILE__ ).find_last_of( "/\\" ) ) + "/../../reference_files/station_history/" +
            name;
}
std::string contents( const std::string& name )
{
    std::ifstream file( fixture( name ) );
    return { std::istreambuf_iterator< char >( file ), std::istreambuf_iterator< char >( ) };
}
std::vector< std::string > lines( const std::string& text )
{
    std::istringstream stream( text );
    std::vector< std::string > result;
    std::string line;
    while( std::getline( stream, line ) ) result.push_back( line );
    return result;
}
SinexStationHistory parse( const std::string& text )
{
    std::istringstream stream( text );
    return readSinexStationHistory( stream, "controlled.snx" );
}
std::string firstSite( )
{
    return lines( contents( "ab07_solutions.snx" ) ).at( 2 );
}
std::string firstParameter( )
{
    for( const auto& line : lines( contents( "ab07_solutions.snx" ) ) )
        if( line.size( ) == 80 && line.substr( 7, 4 ) == "STAX" ) return line;
    throw std::runtime_error( "Fixture has no STAX" );
}
std::string document( const std::vector< std::string >& records, bool addSite = true )
{
    auto header = lines( contents( "ab07_solutions.snx" ) ).front( );
    const auto count = std::to_string( records.size( ) );
    header.replace( 60, 5, std::string( 5 - count.size( ), ' ' ) + count );
    std::string text = header + "\n";
    if( addSite ) text += "+SITE/ID\n" + firstSite( ) + "\n-SITE/ID\n";
    text += "+SOLUTION/ESTIMATE\n";
    for( const auto& record : records ) text += record + "\n";
    return text + "-SOLUTION/ESTIMATE\n%ENDSNX\n";
}
void errorContains( const std::string& text, const std::string& diagnostic )
{
    try
    {
        parse( text );
        BOOST_FAIL( "Expected parser failure" );
    }
    catch( const std::runtime_error& error )
    {
        BOOST_CHECK( std::string( error.what( ) ).find( "controlled.snx:" ) != std::string::npos );
        BOOST_CHECK( std::string( error.what( ) ).find( diagnostic ) != std::string::npos );
    }
}
}  // namespace

BOOST_AUTO_TEST_CASE( real_fixed_width_ssc_preserves_six_solutions_components_and_units )
{
    const auto history = readSinexStationHistory( fixture( "ab07_solutions.snx" ) );
    BOOST_CHECK_EQUAL( history.sites.size( ), 1 );
    BOOST_CHECK_EQUAL( history.parameters.size( ), 36 );
    BOOST_CHECK_EQUAL( history.observationSpans.size( ), 6 );
    BOOST_CHECK_EQUAL( history.discontinuities.size( ), 0 );
    const auto keys = getSinexCartesianSolutionKeys( history );
    BOOST_REQUIRE_EQUAL( keys.size( ), 6 );
    for( const auto& key : keys )
    {
        const auto solution = getSinexCartesianSolution( history, key, true );
        BOOST_CHECK_EQUAL( solution.site.domesId, "40423M002" );
        BOOST_CHECK_EQUAL( solution.key.site.siteCode, "AB07" );
        BOOST_CHECK_EQUAL( solution.key.site.pointCode, "A" );
        BOOST_CHECK_EQUAL( solution.positionReferenceEpoch.civilUtc->year, 2015 );
        BOOST_CHECK_EQUAL( solution.positionReferenceEpoch.civilUtc->dayOfYear, 1 );
        for( int i = 0; i < 3; ++i )
        {
            const auto index = *solution.velocityParameterIndices[ i ];
            BOOST_CHECK_EQUAL( history.parameters.at( index ).originalUnit, "m/y" );
            BOOST_CHECK_CLOSE_FRACTION(
                    *solution.velocityMetresPerSecond[ i ], history.parameters.at( index ).originalValue / 31557600.0, 1e-15 );
        }
    }
    const auto first = getSinexCartesianSolution( history, keys.front( ) );
    BOOST_CHECK_CLOSE_FRACTION( first.positionMetres[ 0 ], -3425750.19097981, 1e-15 );
    BOOST_CHECK_CLOSE_FRACTION( first.positionMetres[ 1 ], -1214685.85537427, 1e-15 );
    BOOST_CHECK_CLOSE_FRACTION( first.positionMetres[ 2 ], 5223662.74510764, 1e-15 );
    BOOST_CHECK_EQUAL( history.parameters.begin( )->second.provenance.source, fixture( "ab07_solutions.snx" ) );
}

BOOST_AUTO_TEST_CASE( raw_observation_spans_and_position_velocity_validity_are_distinct )
{
    const auto observed = readSinexStationHistory( fixture( "ab07_solutions.snx" ) );
    const auto validity = readSinexStationHistory( fixture( "ab07_discontinuities.snx" ) );
    BOOST_CHECK_EQUAL( validity.sites.size( ), 0 );
    BOOST_CHECK_EQUAL( validity.discontinuities.size( ), 10 );
    BOOST_CHECK_EQUAL( std::count_if( validity.discontinuities.begin( ),
                                      validity.discontinuities.end( ),
                                      []( const auto& record ) { return record.kind == "P"; } ),
                       6 );
    BOOST_CHECK_EQUAL( std::count_if( validity.discontinuities.begin( ),
                                      validity.discontinuities.end( ),
                                      []( const auto& record ) { return record.kind == "V"; } ),
                       4 );
    BOOST_CHECK( validity.discontinuities.front( ).start.meaning == SinexEpochMeaning::openStart );
    BOOST_CHECK( !validity.discontinuities.front( ).start.civilUtc );
    BOOST_CHECK( validity.discontinuities[ 5 ].end.meaning == SinexEpochMeaning::openEnd );
    BOOST_CHECK( !( *validity.discontinuities[ 2 ].end.civilUtc == *observed.observationSpans[ 2 ].end.civilUtc ) );
}

BOOST_AUTO_TEST_CASE( parameter_only_psd_and_seasonal_records_preserve_unresolved_provider_units )
{
    const auto psd = readSinexStationHistory( fixture( "ab07_psd.snx" ) );
    const auto seasonal = readSinexStationHistory( fixture( "albh_seasonal.snx" ) );
    BOOST_CHECK_EQUAL( psd.parameters.size( ), 4 );
    BOOST_CHECK_EQUAL( seasonal.parameters.size( ), 72 );
    BOOST_CHECK_EQUAL( getSinexCartesianSolutionKeys( psd ).size( ), 0 );
    BOOST_CHECK_EQUAL( getSinexCartesianSolutionKeys( seasonal ).size( ), 0 );
    for( const auto& row : psd.parameters )
    {
        const auto& parameter = row.second;
        BOOST_CHECK_EQUAL( parameter.key.solutionId, "----" );
        BOOST_CHECK_EQUAL( parameter.originalUnit, "m" );
        if( parameter.type.rfind( "TEXP", 0 ) == 0 )
        {
            BOOST_CHECK( parameter.normalization == SinexNormalization::unresolvedExtension );
            BOOST_CHECK( !parameter.valueSi );
            BOOST_CHECK( !parameter.standardDeviationSi );
            BOOST_CHECK( parameter.siUnit.empty( ) );
        }
    }
    auto repeated = psd;
    auto parameter = psd.parameters.begin( )->second;
    parameter.index = 99;
    repeated.parameters.emplace( 99, parameter );
    BOOST_CHECK_NO_THROW( validateSinexStationHistory( repeated ) );  // Multiple physical PSD terms can share type/event.
}

BOOST_AUTO_TEST_CASE( fixed_columns_and_reordered_blocks_records_are_independent_of_file_order )
{
    const auto source = lines( contents( "ab07_solutions.snx" ) );
    std::vector< std::string > parameters;
    for( const auto& line : source )
        if( line.size( ) == 80 ) parameters.push_back( line );
    std::reverse( parameters.begin( ), parameters.end( ) );
    auto text = document( parameters, false );
    text.insert( text.find( "%ENDSNX" ), "+SITE/ID\n" + firstSite( ) + "\n-SITE/ID\n" );
    const auto shuffled = parse( text );
    const auto original = parse( contents( "ab07_solutions.snx" ) );
    for( const auto& key : getSinexCartesianSolutionKeys( shuffled ) )
    {
        const auto a = getSinexCartesianSolution( original, key, true ), b = getSinexCartesianSolution( shuffled, key, true );
        for( int i = 0; i < 3; ++i )
        {
            BOOST_CHECK_EQUAL( a.positionMetres[ i ], b.positionMetres[ i ] );
            BOOST_CHECK_EQUAL( *a.velocityMetresPerSecond[ i ], *b.velocityMetresPerSecond[ i ] );
        }
    }
}

BOOST_AUTO_TEST_CASE( epoch_pivot_and_field_specific_zero_token_semantics )
{
    auto parameter = firstParameter( );
    parameter.replace( 27, 12, "50:365:86399" );
    BOOST_CHECK_EQUAL( parse( document( { parameter } ) ).parameters.begin( )->second.referenceEpoch.civilUtc->year, 2050 );
    parameter.replace( 27, 12, "51:001:00000" );
    BOOST_CHECK_EQUAL( parse( document( { parameter } ) ).parameters.begin( )->second.referenceEpoch.civilUtc->year, 1951 );
    parameter.replace( 27, 12, "00:000:00000" );
    const auto raw = parse( document( { parameter } ) );
    BOOST_CHECK( raw.parameters.begin( )->second.referenceEpoch.meaning == SinexEpochMeaning::unspecified );
    auto text = document( { firstParameter( ) } );
    text.insert( text.find( "%ENDSNX" ), "+SOLUTION/EPOCHS\n AB07  A    1 C 00:000:00000 00:000:00000 00:000:00000\n-SOLUTION/EPOCHS\n" );
    const auto spans = parse( text ).observationSpans;
    BOOST_CHECK( spans.front( ).start.meaning == SinexEpochMeaning::headerStartDefault );
    BOOST_CHECK( spans.front( ).end.meaning == SinexEpochMeaning::headerEndDefault );
    BOOST_CHECK( !spans.front( ).start.civilUtc );
}

BOOST_AUTO_TEST_CASE( optional_cartesian_components_remain_missing_and_complete_selection_is_explicit )
{
    auto history = parse( document( { firstParameter( ) } ) );
    BOOST_CHECK_EQUAL( history.parameters.size( ), 1 );
    BOOST_CHECK_THROW( getSinexCartesianSolution( history, { { "AB07", "A" }, "1" } ), std::runtime_error );
    const auto complete = readSinexStationHistory( fixture( "ab07_solutions.snx" ) );
    history.parameters.clear( );
    for( const auto& row : complete.parameters )
        if( row.first <= 3 ) history.parameters.insert( row );
    const auto selected = getSinexCartesianSolution( history, { { "AB07", "A" }, "1" } );
    for( int i = 0; i < 3; ++i ) BOOST_CHECK( !selected.velocityMetresPerSecond[ i ] );
    BOOST_CHECK_THROW( getSinexCartesianSolution( history, { { "AB07", "A" }, "1" }, true ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( exponent_variants_nonfinite_units_and_duplicate_components_are_checked )
{
    auto parameter = firstParameter( );
    parameter.replace( parameter.find( 'E', 47 ), 1, "D" );
    BOOST_CHECK_EQUAL( *parse( document( { parameter } ) ).parameters.begin( )->second.valueSi, -3425750.19097981 );
    auto duplicate = parameter;
    duplicate.replace( 1, 5, "   99" );
    errorContains( document( { parameter, duplicate } ), "Duplicate/conflicting Cartesian" );
    errorContains( document( { parameter, parameter } ), "Duplicate parameter index" );
    auto bad = parameter;
    bad.replace( 40, 4, "m/y " );
    errorContains( document( { bad } ), "Invalid unit" );
    bad = parameter;
    bad.replace( 47, 21, "                  nan" );
    errorContains( document( { bad } ), "nonfinite" );
    bad = parameter;
    bad.replace( 69, 11, "-1.0000E-03" );
    errorContains( document( { bad } ), "Negative parameter" );
}

BOOST_AUTO_TEST_CASE( malformed_calendar_record_and_block_diagnostics_include_source_line )
{
    auto bad = firstParameter( );
    bad.replace( 27, 12, "19:366:00000" );
    errorContains( document( { bad } ), "civil UTC" );
    bad.replace( 27, 12, "20:366:00000" );
    BOOST_CHECK_NO_THROW( parse( document( { bad } ) ) );
    bad = firstParameter( );
    bad.pop_back( );
    errorContains( document( { bad } ), "Truncated" );
    errorContains( "", "Empty/truncated" );
    errorContains( "+SITE/ID\n", "requires a SINEX header" );
    auto text = document( { firstParameter( ) } );
    text.erase( text.find( "%ENDSNX" ) );
    errorContains( text, "Empty/truncated" );
    text = document( { firstParameter( ) } );
    text.insert( text.find( "-SITE/ID" ), "+SOLUTION/EPOCHS\n" );
    errorContains( text, "Nested" );
}

BOOST_AUTO_TEST_CASE( preloaded_mutable_records_are_validated_before_cartesian_consumption )
{
    auto history = readSinexStationHistory( fixture( "ab07_solutions.snx" ) );
    history.parameters.begin( )->second.originalValue = std::numeric_limits< double >::infinity( );
    BOOST_CHECK_THROW( validateSinexStationHistory( history ), std::runtime_error );
    history = readSinexStationHistory( fixture( "ab07_solutions.snx" ) );
    history.parameters.begin( )->second.valueSi = 0.0;
    BOOST_CHECK_THROW( getSinexCartesianSolution( history, { { "AB07", "A" }, "1" } ), std::runtime_error );
    history = readSinexStationHistory( fixture( "ab07_solutions.snx" ) );
    history.parameters.begin( )->second.referenceEpoch.civilUtc->year = 2016;
    BOOST_CHECK_THROW( validateSinexStationHistory( history ), std::runtime_error );
    history = readSinexStationHistory( fixture( "ab07_solutions.snx" ) );
    auto duplicate = history.parameters.begin( )->second;
    duplicate.index = 99;
    duplicate.referenceEpoch.original.clear( );
    history.parameters.emplace( 99, duplicate );
    BOOST_CHECK_THROW( validateSinexStationHistory( history ), std::runtime_error );
    history.parameters.begin( )->second.referenceEpoch.original.clear( );
    history.parameters.at( 99 ).referenceEpoch.civilUtc->year = 2016;
    BOOST_CHECK_NO_THROW( validateSinexStationHistory( history ) );  // Distinct civil epochs stay distinct without raw strings.
    history = readSinexStationHistory( fixture( "ab07_discontinuities.snx" ) );
    history.discontinuities.front( ).start.meaning = SinexEpochMeaning::headerStartDefault;
    BOOST_CHECK_THROW( validateSinexStationHistory( history ), std::runtime_error );
    history = readSinexStationHistory( fixture( "ab07_solutions.snx" ) );
    history.headerStart.civilUtc.reset( );
    history.headerStart.original = "00:000:00000";
    history.headerStart.meaning = SinexEpochMeaning::openStart;
    BOOST_CHECK_THROW( validateSinexStationHistory( history ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( eccentricity_component_basis_header_defaults_and_opaque_blocks_are_retained )
{
    auto text = document( {} );
    text.insert( text.find( "%ENDSNX" ),
                 "+SITE/ECCENTRICITY\n AB07  A    1 P 00:000:00000 00:000:00000 UNE   0.1000   0.2000   0.3000\n-SITE/ECCENTRICITY\n"
                 "+SITE/RECEIVER\n* equipment retained\n AB07  A    1 P 15:001:00000 00:000:00000 TEST RECEIVER        12345 FW 1.0     "
                 "\n-SITE/RECEIVER\n"
                 "+SOLUTION/MATRIX_ESTIMATE L COVA\n     1     1  0.123456789012345E-05\n-SOLUTION/MATRIX_ESTIMATE L COVA\n" );
    const auto history = parse( text );
    BOOST_REQUIRE_EQUAL( history.eccentricities.size( ), 1 );
    const auto& eccentricity = history.eccentricities.front( );
    BOOST_CHECK( eccentricity.basis == SinexEccentricityBasis::une );
    BOOST_CHECK_EQUAL( eccentricity.componentsMetres[ 0 ], 0.1 );
    BOOST_CHECK_EQUAL( eccentricity.componentsMetres[ 1 ], 0.2 );
    BOOST_CHECK_EQUAL( eccentricity.componentsMetres[ 2 ], 0.3 );
    BOOST_CHECK( eccentricity.start.meaning == SinexEpochMeaning::headerStartDefault );
    BOOST_CHECK( eccentricity.end.meaning == SinexEpochMeaning::headerEndDefault );
    BOOST_REQUIRE_EQUAL( history.equipmentHistories.size( ), 1 );
    BOOST_CHECK( history.equipmentHistories.front( ).kind == SinexEquipmentKind::receiver );
    BOOST_CHECK_EQUAL( history.equipmentHistories.front( ).equipmentType, "TEST RECEIVER" );
    BOOST_CHECK_EQUAL( history.equipmentHistories.front( ).serialNumber, "12345" );
    BOOST_CHECK_EQUAL( history.equipmentHistories.front( ).firmwareVersion, "FW 1.0" );
    BOOST_CHECK( history.equipmentHistories.front( ).end.meaning == SinexEpochMeaning::headerEndDefault );
    BOOST_REQUIRE_EQUAL( history.opaqueBlocks.size( ), 1 );
    BOOST_CHECK_EQUAL( history.opaqueBlocks[ 0 ].name, "SOLUTION/MATRIX_ESTIMATE L COVA" );
    BOOST_CHECK( !history.opaqueBlocks[ 0 ].records.front( ).text.empty( ) );
}

BOOST_AUTO_TEST_CASE( real_combined_product_duplicate_site_codes_and_repeated_psd_terms_are_preserved )
{
    const auto history = readSinexStationHistory( fixture( "duplicate_sites.snx" ) );
    const auto gala = history.sites.equal_range( { "GALA", "A" } );
    BOOST_REQUIRE_EQUAL( std::distance( gala.first, gala.second ), 2 );
    BOOST_CHECK_EQUAL( gala.first->second.domesId, "42005M001" );
    BOOST_CHECK_EQUAL( std::next( gala.first )->second.domesId, "42004S001" );
    BOOST_CHECK_THROW( getSinexCartesianSolution( history, { { "GALA", "A" }, "1" } ), std::runtime_error );
    const auto yuma = history.sites.equal_range( { "7894", "A" } );
    BOOST_REQUIRE_EQUAL( std::distance( yuma.first, yuma.second ), 2 );
    BOOST_CHECK_EQUAL( yuma.first->second.domesId, std::next( yuma.first )->second.domesId );
    const auto repeated = readSinexStationHistory( fixture( "repeated_psd.snx" ) );
    std::size_t northExponentialAmplitudes = 0;
    for( const auto& row : repeated.parameters )
        if( row.second.key.site.siteCode == "AB48" && row.second.type == "AEXP_N" ) ++northExponentialAmplitudes;
    BOOST_CHECK_EQUAL( northExponentialAmplitudes, 2 );
    BOOST_CHECK_NO_THROW( validateSinexStationHistory( repeated ) );
    // Identical SITE/ID records remain separate raw provenance, but do not make a Cartesian selection ambiguous.
    auto text = contents( "ab07_solutions.snx" );
    text.insert( text.find( "-SITE/ID" ), firstSite( ) + "\n" );
    const auto duplicate = parse( text );
    BOOST_CHECK_EQUAL( duplicate.sites.size( ), 2 );
    BOOST_CHECK_NO_THROW( getSinexCartesianSolution( duplicate, { { "AB07", "A" }, "1" } ) );
}

BOOST_AUTO_TEST_CASE( equipment_histories_validate_standard_columns_and_mutable_interval_roles )
{
    auto text = document( {} );
    std::string antenna( 68, ' ' );
    antenna.replace( 1, 4, "AB07" );
    antenna.replace( 6, 2, " A" );
    antenna.replace( 9, 4, "   1" );
    antenna.replace( 14, 1, "P" );
    antenna.replace( 16, 12, "00:000:00000" );
    antenna.replace( 29, 12, "00:000:00000" );
    antenna.replace( 42, 20, "TEST ANTENNA    NONE" );
    antenna.replace( 63, 5, "12345" );
    text.insert( text.find( "%ENDSNX" ), "+SITE/ANTENNA\n" + antenna + "\n-SITE/ANTENNA\n" );
    auto history = parse( text );
    BOOST_REQUIRE_EQUAL( history.equipmentHistories.size( ), 1 );
    BOOST_CHECK( history.equipmentHistories.front( ).kind == SinexEquipmentKind::antenna );
    BOOST_CHECK_EQUAL( history.equipmentHistories.front( ).equipmentType, "TEST ANTENNA    NONE" );
    BOOST_CHECK( history.equipmentHistories.front( ).firmwareVersion.empty( ) );
    BOOST_CHECK( history.equipmentHistories.front( ).start.meaning == SinexEpochMeaning::headerStartDefault );
    history.equipmentHistories.front( ).end.meaning = SinexEpochMeaning::openEnd;
    BOOST_CHECK_THROW( validateSinexStationHistory( history ), std::runtime_error );
    text.erase( text.find( "12345" ) + 4, 1 );
    errorContains( text, "Truncated" );
}

BOOST_AUTO_TEST_CASE( published_igs_equipment_history_preserves_identity_changes_and_source_fields )
{
    const auto history = readSinexStationHistory( fixture( "abmf_equipment.snx" ) );
    BOOST_CHECK_EQUAL( history.sites.size( ), 1 );
    BOOST_CHECK_EQUAL( history.parameters.size( ), 0 );
    BOOST_CHECK( history.headerStart.meaning == SinexEpochMeaning::unspecified );
    BOOST_REQUIRE_EQUAL( history.equipmentHistories.size( ), 28 );
    const auto& first = history.equipmentHistories.front( );
    BOOST_CHECK_EQUAL( first.key.site.siteCode, "abmf" );
    BOOST_CHECK_EQUAL( first.key.site.pointCode, "A" );
    BOOST_CHECK_EQUAL( first.key.solutionId, "----" );
    BOOST_CHECK_EQUAL( first.equipmentType, "ASHTECH UZ-12" );
    BOOST_CHECK_EQUAL( first.serialNumber, "UC220" );
    BOOST_CHECK_EQUAL( first.firmwareVersion, "CN00" );
    BOOST_CHECK_EQUAL( first.start.civilUtc->year, 2008 );
    BOOST_CHECK_EQUAL( first.start.civilUtc->dayOfYear, 197 );
    BOOST_CHECK_EQUAL( first.provenance.line, 6 );
    BOOST_CHECK_EQUAL( first.provenance.text.size( ), 80 );
    const auto& lastReceiver = history.equipmentHistories[ 24 ];
    BOOST_CHECK( lastReceiver.kind == SinexEquipmentKind::receiver );
    BOOST_CHECK_EQUAL( lastReceiver.firmwareVersion, "5.6.0" );
    BOOST_CHECK( lastReceiver.end.meaning == SinexEpochMeaning::headerEndDefault );
    const auto& antenna = history.equipmentHistories[ 25 ];
    BOOST_CHECK( antenna.kind == SinexEquipmentKind::antenna );
    BOOST_CHECK_EQUAL( antenna.equipmentType, "AERAT2775_43    SPKE" );
    BOOST_CHECK_EQUAL( antenna.serialNumber, "5546" );
    BOOST_CHECK( antenna.firmwareVersion.empty( ) );
    // The IGS extension after the standard 68 columns is retained as raw provenance.
    BOOST_CHECK_EQUAL( antenna.provenance.text.substr( 68 ), "    0" );
    BOOST_REQUIRE_EQUAL( history.eccentricities.size( ), 3 );
    BOOST_CHECK_EQUAL( history.eccentricities.front( ).componentsMetres[ 0 ], 0.05 );
    BOOST_CHECK( history.eccentricities.back( ).end.meaning == SinexEpochMeaning::headerEndDefault );

    auto text = contents( "abmf_equipment.snx" );
    const auto receiver = text.find( first.provenance.text );
    text.replace( receiver + 29, 12, "08:196:00000" );
    errorContains( text, "reversed station interval" );
    text = contents( "abmf_equipment.snx" );
    text.erase( text.find( first.provenance.text ) + 79, 1 );
    errorContains( text, "Truncated" );
}

BOOST_AUTO_TEST_CASE( declared_estimate_counts_and_header_intervals_reject_incomplete_products )
{
    auto text = document( { firstParameter( ) } );
    text.replace( 60, 5, "    2" );
    errorContains( text, "declared parameter count" );
    text.replace( 60, 5, "    0" );
    errorContains( text, "declared parameter count" );
    text = document( {} );
    text.replace( 32, 12, "20:001:00000" );
    text.replace( 45, 12, "19:365:86399" );
    errorContains( text, "reversed station interval" );
}

BOOST_AUTO_TEST_CASE( published_rounded_approximate_seconds_carry_exactly_with_raw_fields_retained )
{
    const auto history = readSinexStationHistory( fixture( "cedu_equipment.snx" ) );
    BOOST_REQUIRE_EQUAL( history.sites.size( ), 1 );
    const auto& site = history.sites.begin( )->second;
    BOOST_CHECK_EQUAL( site.provenance.text.substr( 63, 4 ), "60.0" );
    BOOST_REQUIRE( site.approximateLatitudeRadians );
    BOOST_CHECK_CLOSE_FRACTION( *site.approximateLatitudeRadians, -( 31.0 + 52.0 / 60.0 ) * std::acos( -1.0 ) / 180.0, 2e-16 );
    auto text = contents( "cedu_equipment.snx" );
    text.replace( text.find( site.provenance.text ) + 63, 4, "60.1" );
    errorContains( text, "Invalid approximate DMS coordinate" );
}

BOOST_AUTO_TEST_CASE( malformed_published_igs_site_fields_are_rejected_without_column_repair )
{
    // Pinned IGS template line 389 has a 21-column rather than 22-column
    // description. Preserve the source row and reject its shifted coordinates.
    const auto text = contents( "igs_malformed_site.snx" );
    errorContains( text, "Truncated SITE/ID coordinates" );
}
