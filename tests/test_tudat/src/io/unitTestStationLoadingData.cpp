#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readStationLoadingData.h"

using namespace tudat::input_output;

//! Preserve the published IMLS station coordinates, signed REN displacements and declared TAI epoch.
BOOST_AUTO_TEST_CASE( published_station_loading_coordinates_and_displacements )
{
    const auto data = readEphedispStationLoadingData( tudat::paths::getTudatTestDataPath( ) +
                                                      "pr1092/src/io/station_loading_fixtures/imls_merra2_extract.eph" );

    // The epoch is the declared TAI epoch (MJD 58849).
    BOOST_REQUIRE_EQUAL( data.epochs.size( ), 1 );
    BOOST_CHECK( data.epochs[ 0 ] == tudat::Time( ( 58849 - 51544 ) * 24 - 12, 0.0L ) );

    const auto& site = data.sites.at( "weuc" );

    // The station and its sample retain metres and the provider's radial/east/north signs.
    BOOST_CHECK_EQUAL( site.cartesianPosition[ 0 ], 4075560.075 );
    BOOST_CHECK_EQUAL( site.cartesianPosition[ 1 ], 931580.295 );
    BOOST_CHECK_EQUAL( site.cartesianPosition[ 2 ], 4801621.958 );
    BOOST_REQUIRE_EQUAL( site.displacements.size( ), 1 );
    BOOST_CHECK_EQUAL( site.displacements.at( 0 )[ 0 ], -.00510 );
    BOOST_CHECK_EQUAL( site.displacements.at( 0 )[ 1 ], -.00202 );
    BOOST_CHECK_EQUAL( site.displacements.at( 0 )[ 2 ], .00170 );
}
