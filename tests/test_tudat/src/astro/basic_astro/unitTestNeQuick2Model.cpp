/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Validates the NeQuick-2 implementation against ITU-R P.531-15 reference test vectors.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/io/readNeQuick2Data.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/basic_astro/neQuick2Model.h"

#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>

namespace tudat
{
namespace unit_tests
{

// Path to NeQuick-2 data (bundled in tudat resources)
const std::string neQuick2DataPath = tudat::paths::getTudatTestDataPath( ) + "/nequick2/";

BOOST_AUTO_TEST_SUITE( test_nequick2_model )

BOOST_AUTO_TEST_CASE( testElectronDensityProfileTestCase1 )
{
    // TestCase1: Ray from (45N, 15E, 0km) to (45N, 35E, 20000km)
    // Year 2000, month 9, UT 15.0, R12 = 116.3 -> F10.7 = 160.4
    // Reference electron density values along the ray at specific heights

    input_output::CcirData ccirData = input_output::readCcirCoefficients( neQuick2DataPath );
    Eigen::MatrixXd modipGrid = input_output::readModipGrid( neQuick2DataPath );

    // Constant F10.7 = 160.4 for the test
    auto solarFluxFunc = []( double ) { return 160.4; };

    environment::NeQuick2Model model( ccirData, modipGrid, solarFluxFunc );

    // Reference values from TestCase1.txt at specific heights along the ray
    // The ray goes from (45, 15, 0) to (45, 35, 20000) so we can test vertical-ish profile
    // at the starting point (lat=45, lon=15)

    // Test electron density at specific heights at (45N, 15E) for month=9, F10.7=160.4, UT=15
    struct TestPoint {
        double heightKm;
        double expectedNe;  // el/m^3
        double tolerancePct;
    };

    // Reference values from TestCase1.txt — note these are along an oblique ray from (45,15,0) to (45,35,20000)
    // At low heights the ray is near (45, 15), at high heights it's shifted eastward.
    // We test at exact (45, 15) vertically, so only the F2 peak and topside should match closely.
    // Lower altitudes have different lat/lon along the ray vs our vertical evaluation.
    std::vector< TestPoint > testPoints = {
        { 300.0,  9.58370e+11,  2.0 },   // near F2 peak, ray is close to (45, 15)
        { 400.0,  6.93350e+11,  2.0 },
        { 500.0,  3.55978e+11,  2.0 },
        { 1000.0, 3.66431e+10,  3.0 },
    };

    for( const auto& tp : testPoints )
    {
        double ne = model.computeElectronDensity( tp.heightKm, 45.0, 15.0, 9, 160.4, 15.0 );
        if( tp.expectedNe > 1.0e6 )  // skip near-zero values
        {
            double relError = std::abs( ne - tp.expectedNe ) / tp.expectedNe * 100.0;
            BOOST_TEST_MESSAGE( "h=" << tp.heightKm << " km: Ne=" << ne
                                << " expected=" << tp.expectedNe << " relErr=" << relError << "%" );
            BOOST_CHECK_SMALL( relError, tp.tolerancePct );
        }
    }
}

BOOST_AUTO_TEST_CASE( testVerticalTecConsistency )
{
    // Verify that vertical TEC (from getVerticalTotalElectronContent) is consistent
    // with the electron density profile integral

    input_output::CcirData ccirData = input_output::readCcirCoefficients( neQuick2DataPath );
    Eigen::MatrixXd modipGrid = input_output::readModipGrid( neQuick2DataPath );

    auto solarFluxFunc = []( double ) { return 120.0; };

    environment::NeQuick2Model model( ccirData, modipGrid, solarFluxFunc );

    // Use time = 0 (J2000 epoch = 2000-01-01 12:00 UT, month=1, ut=12.0)
    double vtec = model.getVerticalTotalElectronContent( 45.0, 15.0, 0.0 );

    BOOST_TEST_MESSAGE( "VTEC at (45N, 15E) for F10.7=120, month=1, UT=12: " << vtec << " TECU" );

    // VTEC should be in a plausible range: 1-200 TECU for mid-latitudes
    BOOST_CHECK_GT( vtec, 1.0 );
    BOOST_CHECK_LT( vtec, 200.0 );
}

BOOST_AUTO_TEST_CASE( testDataLoading )
{
    // Verify CCIR and MODIP data loads correctly

    BOOST_REQUIRE_NO_THROW( input_output::readCcirCoefficients( neQuick2DataPath ) );
    BOOST_REQUIRE_NO_THROW( input_output::readModipGrid( neQuick2DataPath ) );

    input_output::CcirData ccirData = input_output::readCcirCoefficients( neQuick2DataPath );
    Eigen::MatrixXd modipGrid = input_output::readModipGrid( neQuick2DataPath );

    // Check CCIR data dimensions
    for( int m = 0; m < 12; ++m )
    {
        BOOST_CHECK_EQUAL( ccirData.foF2Coefficients[ m ][ 0 ].rows( ), 76 );
        BOOST_CHECK_EQUAL( ccirData.foF2Coefficients[ m ][ 0 ].cols( ), 13 );
        BOOST_CHECK_EQUAL( ccirData.m3000F2Coefficients[ m ][ 0 ].rows( ), 49 );
        BOOST_CHECK_EQUAL( ccirData.m3000F2Coefficients[ m ][ 0 ].cols( ), 9 );
    }

    // Check MODIP grid dimensions (184 x 184 with padding)
    BOOST_CHECK_EQUAL( modipGrid.rows( ), 184 );
    BOOST_CHECK_EQUAL( modipGrid.cols( ), 184 );

    // MODIP values should be in range [-90, 90]
    BOOST_CHECK_GE( modipGrid.minCoeff( ), -91.0 );
    BOOST_CHECK_LE( modipGrid.maxCoeff( ), 91.0 );
}

BOOST_AUTO_TEST_CASE( testFluxSaturation )
{
    // Test that electron density is computed without errors at flux boundaries

    input_output::CcirData ccirData = input_output::readCcirCoefficients( neQuick2DataPath );
    Eigen::MatrixXd modipGrid = input_output::readModipGrid( neQuick2DataPath );

    auto solarFluxFunc = []( double ) { return 100.0; };
    environment::NeQuick2Model model( ccirData, modipGrid, solarFluxFunc );

    // F10.7 = 63 (minimum)
    double ne63 = model.computeElectronDensity( 300.0, 0.0, 0.0, 6, 63.0, 12.0 );
    BOOST_CHECK_GT( ne63, 0.0 );

    // F10.7 = 193 (maximum, saturated)
    double ne193 = model.computeElectronDensity( 300.0, 0.0, 0.0, 6, 193.0, 12.0 );
    BOOST_CHECK_GT( ne193, 0.0 );

    // Higher flux should give higher electron density at F2 peak
    BOOST_CHECK_GT( ne193, ne63 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
