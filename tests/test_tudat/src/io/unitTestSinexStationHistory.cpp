#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <iterator>
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readSinexStationHistory.h"

using namespace tudat::input_output;
namespace
{

//! Locate the published SINEX extracts unpacked with the other test data.
std::string fixture( const std::string& name )
{
    return tudat::paths::getTudatTestDataPath( ) + "pr1092/reference_files/station_history/" + name;
}
}  // namespace

//! Preserve the published Cartesian solutions, independent P/V intervals and provider-specific motion units.
BOOST_AUTO_TEST_CASE( published_coordinate_and_motion_products )
{
    const auto history = readSinexStationHistory( fixture( "ab07_solutions.snx" ) );

    const auto keys = getSinexStationCoordinateSolutionKeys( history );

    // Six solutions retain the published Cartesian coordinates and convert m/year velocities to SI.
    BOOST_REQUIRE_EQUAL( keys.size( ), 6 );
    for( const auto& key : keys )
    {
        const auto solution = getSinexStationCoordinateSolution( history, key, true );

        BOOST_CHECK_EQUAL( solution.site.domesId, "40423M002" );
        for( int i = 0; i < 3; ++i )
        {
            const auto& parameter = history.parameters.at( *solution.velocityParameterIndices[ i ] );

            BOOST_CHECK_EQUAL( parameter.originalUnit, "m/y" );
            BOOST_CHECK_CLOSE_FRACTION( *solution.velocityMetresPerSecond[ i ], parameter.originalValue / 31557600.0, 1e-15 );
        }
    }
    const auto first = getSinexStationCoordinateSolution( history, keys.front( ) );

    BOOST_CHECK_CLOSE_FRACTION( first.positionMetres[ 0 ], -3425750.19097981, 1e-15 );
    BOOST_CHECK_CLOSE_FRACTION( first.positionMetres[ 1 ], -1214685.85537427, 1e-15 );
    BOOST_CHECK_CLOSE_FRACTION( first.positionMetres[ 2 ], 5223662.74510764, 1e-15 );

    const auto validity = readSinexStationHistory( fixture( "ab07_discontinuities.snx" ) );

    // Observation spans do not replace the separate six position and four velocity validity intervals.
    BOOST_CHECK_EQUAL( validity.discontinuities.size( ), 10 );
    BOOST_CHECK_EQUAL(
            std::count_if(
                    validity.discontinuities.begin( ), validity.discontinuities.end( ), []( const auto& row ) { return row.kind == "P"; } ),
            6 );
    BOOST_CHECK_EQUAL(
            std::count_if(
                    validity.discontinuities.begin( ), validity.discontinuities.end( ), []( const auto& row ) { return row.kind == "V"; } ),
            4 );
    BOOST_CHECK( !( *validity.discontinuities[ 2 ].end.civilUtc == *history.observationSpans[ 2 ].end.civilUtc ) );

    const auto psd = readSinexStationHistory( fixture( "ab07_psd.snx" ) );

    const auto seasonal = readSinexStationHistory( fixture( "albh_seasonal.snx" ) );

    // PSD relaxation times retain their printed units until the consuming source profile resolves them.
    BOOST_CHECK_EQUAL( psd.parameters.size( ), 4 );
    BOOST_CHECK_EQUAL( seasonal.parameters.size( ), 72 );
    for( const auto& row : psd.parameters )
    {
        if( row.second.type.rfind( "TEXP", 0 ) == 0 )
        {
            BOOST_CHECK_EQUAL( row.second.originalUnit, "m" );
            BOOST_CHECK( row.second.normalization == SinexNormalization::unresolvedExtension );
            BOOST_CHECK( !row.second.valueSi );
        }
    }
}

//! Preserve distinct monuments sharing a site code and multiple physical PSD terms sharing a type.
BOOST_AUTO_TEST_CASE( published_site_and_event_identities )
{
    const auto history = readSinexStationHistory( fixture( "duplicate_sites.snx" ) );

    const auto gala = history.sites.equal_range( { "GALA", "A" } );

    // Different DOMES identities cannot be silently selected by their shared four-character site code.
    BOOST_REQUIRE_EQUAL( std::distance( gala.first, gala.second ), 2 );
    BOOST_CHECK_EQUAL( gala.first->second.domesId, "42005M001" );
    BOOST_CHECK_EQUAL( std::next( gala.first )->second.domesId, "42004S001" );
    BOOST_CHECK_THROW( getSinexStationCoordinateSolution( history, { { "GALA", "A" }, "1" } ), std::runtime_error );

    const auto repeated = readSinexStationHistory( fixture( "repeated_psd.snx" ) );

    // Original parameter indices keep both exponential amplitudes rather than overwriting one.
    BOOST_CHECK_EQUAL(
            std::count_if( repeated.parameters.begin( ),
                           repeated.parameters.end( ),
                           []( const auto& row ) { return row.second.key.site.siteCode == "AB48" && row.second.type == "AEXP_N"; } ),
            2 );
}

//! Retain actual IGS equipment identifiers, reference vectors and standard-column extensions.
BOOST_AUTO_TEST_CASE( published_igs_equipment_history )
{
    const auto history = readSinexStationHistory( fixture( "abmf_equipment.snx" ) );

    // All published receiver and antenna changes remain available for station construction.
    BOOST_REQUIRE_EQUAL( history.equipmentHistories.size( ), 28 );

    const auto& receiver = history.equipmentHistories.front( );

    const auto& antenna = history.equipmentHistories[ 25 ];

    // Lowercase station codes and embedded equipment spaces are provider identities, not normalization hints.
    BOOST_CHECK_EQUAL( receiver.key.site.siteCode, "abmf" );
    BOOST_CHECK_EQUAL( receiver.equipmentType, "ASHTECH UZ-12" );
    BOOST_CHECK_EQUAL( receiver.serialNumber, "UC220" );
    BOOST_CHECK_EQUAL( receiver.firmwareVersion, "CN00" );
    BOOST_CHECK_EQUAL( antenna.equipmentType, "AERAT2775_43    SPKE" );
    BOOST_CHECK_EQUAL( antenna.serialNumber, "5546" );

    // Equipment changes preserve the published eccentricity and trailing IGS source fields.
    BOOST_CHECK_EQUAL( antenna.provenance.text.substr( 68 ), "    0" );
    BOOST_REQUIRE_EQUAL( history.eccentricities.size( ), 3 );
    BOOST_CHECK_EQUAL( history.eccentricities.front( ).componentsMetres[ 0 ], 0.05 );
    BOOST_CHECK( history.eccentricities.back( ).end.meaning == SinexEpochMeaning::headerEndDefault );
}
