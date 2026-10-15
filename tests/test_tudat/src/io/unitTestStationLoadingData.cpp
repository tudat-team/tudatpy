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

    // Epochs and the displacement basis come from the product, without inventing a seasonal inventory.
    BOOST_REQUIRE_EQUAL( data.epochs.size( ), 1 );
    BOOST_CHECK_EQUAL( data.epochs[ 0 ].modifiedJulianDay, 58849 );
    BOOST_CHECK_EQUAL( data.epochs[ 0 ].secondsOfDay, 0.0L );
    BOOST_CHECK( data.timeScale == StationLoadingTimeScale::tai );
    BOOST_CHECK( data.basis == StationLoadingBasis::geocentricRadialEastNorth );
    BOOST_CHECK( data.metadata.annual == StationLoadingSignalPresence::unknown );

    const auto& site = data.sites.at( "weuc" );

    // The station and its sample retain metres and the provider's radial/east/north signs.
    BOOST_CHECK_EQUAL( site.cartesianPosition[ 0 ], 4075560.075 );
    BOOST_CHECK_EQUAL( site.cartesianPosition[ 1 ], 931580.295 );
    BOOST_CHECK_EQUAL( site.cartesianPosition[ 2 ], 4801621.958 );
    BOOST_REQUIRE_EQUAL( site.samples.size( ), 1 );
    BOOST_CHECK_EQUAL( site.samples.at( 0 ).displacement[ 0 ], -.00510 );
    BOOST_CHECK_EQUAL( site.samples.at( 0 ).displacement[ 1 ], -.00202 );
    BOOST_CHECK_EQUAL( site.samples.at( 0 ).displacement[ 2 ], .00170 );
}
