/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/readAod1bGravityField.h"

namespace
{

//! Resolve the labelled provider extract independently of the test's working directory.
std::string fixture( )
{
    const std::string source = __FILE__;

    return source.substr( 0, source.find_last_of( "/\\" ) ) +
            "/../../reference_files/aod1b_gravity/AOD1B_2026-01-01_degree2_extract_07.asc";
}

}  // namespace

//! Preserve the published GPS dates, header constants and the physical ATM+OCN=GLO relation.
BOOST_AUTO_TEST_CASE( publishedRl07RecordsAndComponentSum )
{
    const auto data = tudat::input_output::readAod1bGravityField( fixture( ) );

    BOOST_CHECK_EQUAL( data.release, 7 );
    BOOST_CHECK_EQUAL( data.softwareVersion, "atm_ocean_dealise.07_v1" );
    BOOST_CHECK_EQUAL( data.gravitationalParameter, 3.986004418e14 );
    BOOST_CHECK_EQUAL( data.referenceRadius, 6378136.6 );
    BOOST_CHECK( data.fullyNormalized );
    BOOST_REQUIRE_EQUAL( data.samples.size( ), 12 );
    BOOST_CHECK_EQUAL( data.headerLines.back( ), "END OF HEADER" );
    BOOST_CHECK_EQUAL( data.samples.front( ).gpsDate, "2026-01-01 00:00:00" );
    BOOST_CHECK_EQUAL( data.samples.front( ).gpsSecondsSinceJ2000, 820497600.0 );
    BOOST_CHECK_EQUAL( data.samples.front( ).coefficients.at( { 2, 0 } ).values[ 0 ], 0.805907484e-10 );
    BOOST_CHECK_EQUAL( data.samples.front( ).coefficients.at( { 2, 0 } ).sourceLine, 35 );

    // Three independently printed nine-digit fields allow at most 1.5e-18 absolute rounding here.
    for( std::size_t epoch = 0; epoch < 3; ++epoch )
    {
        const auto& atmosphere = data.samples.at( epoch * 4 ).coefficients;

        const auto& ocean = data.samples.at( epoch * 4 + 1 ).coefficients;

        const auto& global = data.samples.at( epoch * 4 + 2 ).coefficients;

        for( const auto& entry : global )
        {
            for( std::size_t column = 0; column < 2; ++column )
            {
                BOOST_CHECK_SMALL( entry.second.values[ column ] - atmosphere.at( entry.first ).values[ column ] -
                                           ocean.at( entry.first ).values[ column ],
                                   1.5e-18 );
            }
        }
    }
}

//! Reject incompatible releases and missing supplied coefficients instead of filling them silently.
BOOST_AUTO_TEST_CASE( mutablePreloadedRecordsRetainPresenceAndReleaseRequirements )
{
    auto data = tudat::input_output::readAod1bGravityField( fixture( ) );

    data.release = 6;
    BOOST_CHECK_THROW( tudat::input_output::validateAod1bGravityFieldData( data ), std::runtime_error );
    data.release = 7;
    data.samples.front( ).coefficients.erase( { 2, 1 } );
    BOOST_CHECK_THROW( tudat::input_output::validateAod1bGravityFieldData( data ), std::runtime_error );
}
