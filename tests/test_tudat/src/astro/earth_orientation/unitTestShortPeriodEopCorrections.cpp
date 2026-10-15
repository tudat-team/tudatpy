/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/io/basicInputOutput.h"

#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"
#include "tudat/astro/earth_orientation/precessionNutationCalculator.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "short_period_reference/fixtures.h"

namespace tudat
{

namespace unit_tests
{

using namespace basic_astrodynamics;
using namespace sofa_interface;
using namespace earth_orientation;
using namespace input_output;
using namespace unit_conversions;

//! NOTE on tolerances: 1 mm position error on Earth's surface due to single angle error in Earth orientation corresponds to
//! 0.15 nrad angle uncertainty. This corresponds to 32 microarcseconds. For Earth rotation, this corresponds to 2.2 microseconds
//! in UT1. The applied tolerances (1 microarc seconds and 50 ns) correspond to about 20-30 microns difference at Earth's surface.
BOOST_AUTO_TEST_SUITE( test_short_period_eop_corrections )

//! Test short-periodic ut1-utc variations by comparing to test output of UTLIBR.f file
BOOST_AUTO_TEST_CASE( testShortPeriodLibrationalPolarMotion )
{
    // Create polar motion correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionCalculator(
            convertArcSecondsToRadians< double >( 1.0E-6 ),
            0.0,
            { getEarthOrientationDataFilesPath( ) + "/polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
            { getEarthOrientationDataFilesPath( ) + "/polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" } );

    //  Define test time
    double testMjd = 54335.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute polar motion correction
    double microAsToRadians = mathematical_constants::PI / ( 180.0 * 1.0E6 * 3600.0 );
    Eigen::Vector2d polarMotionCorrections = polarMotionCalculator.getCorrections( testEphemerisTime ) / microAsToRadians;

    // Compare with the corrected 2011 PMSDNUT2 reference. The 0.01 microarcsecond
    // tolerance covers its single-MJD time/argument approximation; see independent fixtures below.
    BOOST_CHECK_SMALL( polarMotionCorrections( 0 ) - 24.83144238273364834, 0.01 );
    BOOST_CHECK_SMALL( polarMotionCorrections( 1 ) + 14.09240692041837661, 0.01 );
}

//! Test short-periodic polar motion libration variations by comparing to test output of ORTHO_EOP.f file
BOOST_AUTO_TEST_CASE( testShortPeriodOceanTidePolarMotion )
{
    // Create polar motion correction ocean tide corrections
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionCalculator(
            convertArcSecondsToRadians< double >( 1.0E-6 ),
            0.0,
            {
                    getEarthOrientationDataFilesPath( ) + "/polarMotionOceanTidesAmplitudes.txt",
            },
            { getEarthOrientationDataFilesPath( ) + "/polarMotionOceanTidesFundamentalArgumentMultipliers.txt" } );

    //  Define test time
    double testMjd = 47100.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute polar motion correction
    double microAsToRadians = mathematical_constants::PI / ( 180.0 * 1.0E6 * 3600.0 );
    Eigen::Vector2d polarMotionCorrections = polarMotionCalculator.getCorrections( testEphemerisTime ) / microAsToRadians;

    // Compare against IERS reference code. Difference between IERS code and this code occurs since the reference uses a
    // different algorothm for calculation (ortho-weights vs. Delaunay arguments). Difference (1 micro arc seconds) is well below
    // observable threshold
    BOOST_CHECK_SMALL( polarMotionCorrections.x( ) + 162.8386373279636530, 1.0 );
    BOOST_CHECK_SMALL( polarMotionCorrections.y( ) - 117.7907525842668974, 1.0 );
}

//! Test short-periodic polar motion variations by checking if multiple contributions are properly combined
BOOST_AUTO_TEST_CASE( testShortPeriodPolarMotionCombinedCorrections )
{
    // Create polar motion correction ocean tide corrections
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionOceanTideCorrectionCalculator(
            convertArcSecondsToRadians< double >( 1.0E-6 ),
            0.0,
            { getEarthOrientationDataFilesPath( ) + "/polarMotionOceanTidesAmplitudes.txt" },
            { getEarthOrientationDataFilesPath( ) + "/polarMotionOceanTidesFundamentalArgumentMultipliers.txt" } );

    // Create polar motion correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionLibrationCorrectionCalculator(
            convertArcSecondsToRadians< double >( 1.0E-6 ),
            0.0,
            { getEarthOrientationDataFilesPath( ) + "/polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
            { getEarthOrientationDataFilesPath( ) + "/polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" } );

    // Create polar motion correction combined corrections
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionTotalCorrectionCalculator(
            convertArcSecondsToRadians< double >( 1.0E-6 ),
            0.0,
            { getEarthOrientationDataFilesPath( ) + "/polarMotionOceanTidesAmplitudes.txt",
              getEarthOrientationDataFilesPath( ) + "/polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
            { getEarthOrientationDataFilesPath( ) + "/polarMotionOceanTidesFundamentalArgumentMultipliers.txt",
              getEarthOrientationDataFilesPath( ) + "/polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" } );

    //  Define test time
    double testMjd = 47100.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute corrections and check validity of combination
    Eigen::Vector2d polarMotionCorrectionOceanTides = polarMotionOceanTideCorrectionCalculator.getCorrections( testEphemerisTime );
    Eigen::Vector2d polarMotionCorrectionLibration = polarMotionLibrationCorrectionCalculator.getCorrections( testEphemerisTime );
    Eigen::Vector2d polarMotionCorrectionTotal = polarMotionTotalCorrectionCalculator.getCorrections( testEphemerisTime );

    BOOST_CHECK_SMALL(
            std::fabs( polarMotionCorrectionTotal( 0 ) - ( polarMotionCorrectionLibration( 0 ) + polarMotionCorrectionOceanTides( 0 ) ) ),
            1.0E-24 );
    BOOST_CHECK_SMALL(
            std::fabs( polarMotionCorrectionTotal( 1 ) - ( polarMotionCorrectionLibration( 1 ) + polarMotionCorrectionOceanTides( 1 ) ) ),
            1.0E-24 );
}

//! Test short-periodic ut1-utc variations by comparing to test output of UTLIBR.f file
BOOST_AUTO_TEST_CASE( testShortPeriodLibrationalUt1 )
{
    // Create UT1 correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator< double > ut1CorrectionCalculator(
            1.0E-6,
            0.0,
            { getEarthOrientationDataFilesPath( ) + "/utcLibrationAmplitudes.txt" },
            { getEarthOrientationDataFilesPath( ) + "/utcLibrationFundamentalArgumentMultipliers.txt" } );

    //  Define test time
    double testMjd = 44239.1;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute UT1 correction
    double ut1Correction = ut1CorrectionCalculator.getCorrections( testEphemerisTime );

    // Compare against IERS reference code. Difference between IERS code and this code occurs since the reference uses a slightly
    // different implementation for computation. Difference (10 ns) is well below observable threshold for Earth rotation
    BOOST_CHECK_SMALL( ut1Correction - 2.441143834386761746E-6, 1.0E-8 );

    //  Define second test time
    testMjd = 55227.4;
    testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute UT1 correction
    ut1Correction = ut1CorrectionCalculator.getCorrections( testEphemerisTime );

    // Compare against IERS reference code. Difference between IERS code and this code occurs since the reference uses a slightly
    // different implementation for computation. Difference (10 ns) is well below observable threshold for Earth rotation
    BOOST_CHECK_SMALL( ut1Correction + 2.655705844335680244E-6, 5.0E-8 );
}

//! Test short-periodic ut1-utc libration variations by comparing to test output of ORTHO_EOP.f file
BOOST_AUTO_TEST_CASE( testShortPeriodOceanTideUt1 )
{
    // Create UT1 correction ocean tide corrections
    ShortPeriodEarthOrientationCorrectionCalculator< double > ut1CorrectionCalculator =
            ShortPeriodEarthOrientationCorrectionCalculator< double >(
                    1.0E-6,
                    0.0,
                    {
                            getEarthOrientationDataFilesPath( ) + "/utcOceanTidesAmplitudes.txt",
                    },
                    { getEarthOrientationDataFilesPath( ) + "/utcOceanTidesFundamentalArgumentMultipliers.txt" } );

    //  Define test time
    double testMjd = 47100.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );
    double ut1Correction = ut1CorrectionCalculator.getCorrections( testEphemerisTime );

    // Compare against IERS reference code. Difference between IERS code and this code occurs since the reference uses a
    // different algorothm for calculation (ortho-weights vs. Delaunay arguments). Difference (50 ns) is well below observable
    // threshold for Earth rotation
    BOOST_CHECK_SMALL( ut1Correction + 23.39092370609808214E-6, 5.0E-8 );
}

//! Test short-periodic ut1-utc variations by checking if multiple contributions are properly combined
BOOST_AUTO_TEST_CASE( testShortPeriodUt1CombinedCorrections )
{
    // Create UT1 correction ocean tide corrections
    ShortPeriodEarthOrientationCorrectionCalculator< double > ut1OceanTideCorrectionCalculator =
            ShortPeriodEarthOrientationCorrectionCalculator< double >(
                    1.0E-6,
                    0.0,
                    { getEarthOrientationDataFilesPath( ) + "/utcOceanTidesAmplitudes.txt" },
                    { getEarthOrientationDataFilesPath( ) + "/utcOceanTidesFundamentalArgumentMultipliers.txt" } );

    // Create UT1 correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator< double > ut1LibrationCorrectionCalculator(
            1.0E-6,
            0.0,
            { getEarthOrientationDataFilesPath( ) + "/utcLibrationAmplitudes.txt" },
            { getEarthOrientationDataFilesPath( ) + "/utcLibrationFundamentalArgumentMultipliers.txt" } );

    // Create UT1 correction combined corrections
    ShortPeriodEarthOrientationCorrectionCalculator< double > ut1TotalCorrectionCalculator(
            1.0E-6,
            0.0,
            { getEarthOrientationDataFilesPath( ) + "/utcLibrationAmplitudes.txt",
              getEarthOrientationDataFilesPath( ) + "/utcOceanTidesAmplitudes.txt" },
            { getEarthOrientationDataFilesPath( ) + "/utcLibrationFundamentalArgumentMultipliers.txt",
              getEarthOrientationDataFilesPath( ) + "/utcOceanTidesFundamentalArgumentMultipliers.txt" } );

    //  Define test time
    double testMjd = 47100.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute corrections and check validity of combination
    double ut1CorrectionOceanTides = ut1OceanTideCorrectionCalculator.getCorrections( testEphemerisTime );
    double ut1CorrectionLibration = ut1LibrationCorrectionCalculator.getCorrections( testEphemerisTime );
    double ut1CorrectionTotal = ut1TotalCorrectionCalculator.getCorrections( testEphemerisTime );

    BOOST_CHECK_SMALL( std::fabs( ut1CorrectionTotal - ( ut1CorrectionLibration + ut1CorrectionOceanTides ) ), 1.0E-20 );
}

BOOST_AUTO_TEST_CASE( test_independent_short_period_references_across_nodal_cycle )
{
    const std::string data = getEarthOrientationDataFilesPath( ) + "/";
    const std::string sourceFile = __FILE__;
    const std::string fixtures = sourceFile.substr( 0, sourceFile.find_last_of( "/\\" ) + 1 ) + "short_period_reference/";
    const double uas = convertArcSecondsToRadians< double >( 1e-6 );
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarLibration(
            uas,
            0.0,
            { data + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
            { data + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" } );
    ShortPeriodEarthOrientationCorrectionCalculator< double > ut1Libration(
            1e-6, 0.0, { data + "utcLibrationAmplitudes.txt" }, { data + "utcLibrationFundamentalArgumentMultipliers.txt" } );
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > oceanPolar(
            uas,
            0.0,
            { data + "polarMotionOceanTidesAmplitudes.txt" },
            { data + "polarMotionOceanTidesFundamentalArgumentMultipliers.txt" } );
    ShortPeriodEarthOrientationCorrectionCalculator< double > oceanUt1(
            1e-6, 0.0, { data + "utcOceanTidesAmplitudes.txt" }, { data + "utcOceanTidesFundamentalArgumentMultipliers.txt" } );
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > referenceUt1Lod(
            1e-6, 0.0, { fixtures + "libration_ut1_lod_amplitudes.txt" }, { fixtures + "libration_ut1_lod_multipliers.txt" } );
    const auto totalPolar = getDefaultPolarMotionCorrectionCalculator( );
    const auto totalUt1 = getDefaultUT1CorrectionCalculator( );
    for( const auto& fixture : short_period_reference::fixtures )
    {
        BOOST_TEST_CONTEXT( "UTC MJD " << fixture.mjdUtc )
        {
            Eigen::Vector6d arguments;
            for( int i = 0; i < 6; ++i )
            {
                arguments( i ) = fixture.arguments[ i ];
            }
            const Eigen::Vector2d polarExact = polarLibration.getCorrectionsFromFundamentalArgument( arguments ) / uas;
            const double ut1Exact = ut1Libration.getCorrectionsFromFundamentalArgument( arguments );

            // Published Table 5.1b T2 sine =0.05 us, while UTLIBR.F uses0.04.
            // This exact, independently generated phase accounts for that sole
            // UT1 coefficient difference; all other residuals remain tightly tested.
            const double t2Phase =
                    2.0 * arguments( 0 ) - arguments( 2 ) - 2.0 * arguments( 3 ) + 2.0 * arguments( 4 ) - 2.0 * arguments( 5 );
            const double acceptedUt1 = fixture.ut1LibrationUs * 1e-6 + 1e-8 * std::sin( t2Phase );
            BOOST_CHECK_SMALL( ut1Exact - acceptedUt1, 1e-12 );
            const Eigen::Vector2d ut1LodExact = referenceUt1Lod.getCorrectionsFromFundamentalArgument( arguments );
            BOOST_CHECK_SMALL( ut1LodExact.x( ) - fixture.ut1LibrationUs * 1e-6, 1e-12 );
            BOOST_CHECK_SMALL( ut1LodExact.y( ) - fixture.lodLibrationUs * 1e-6, 1e-11 );
            const Eigen::Vector2d oceanExact = oceanPolar.getCorrectionsFromFundamentalArgument( arguments ) / uas;
            const double oceanUt1Exact = oceanUt1.getCorrectionsFromFundamentalArgument( arguments );
            const Eigen::Vector2d polarDefault = polarLibration.getCorrections( fixture.ttSeconds ) / uas;
            const Eigen::Vector2d oceanDefault = oceanPolar.getCorrections( fixture.ttSeconds ) / uas;
            const Eigen::Vector2d summedDefault = totalPolar->getCorrections( fixture.ttSeconds ) / uas;
            for( int i = 0; i < 2; ++i )
            {
                BOOST_CHECK_SMALL( polarExact( i ) - fixture.polarLibration[ i ], 1e-5 );
                BOOST_CHECK_SMALL( oceanExact( i ) - fixture.ocean[ i ] * 1e6, 1e-3 );

                // Physical production path uses TT lunisolar arguments and UTC≈UT1
                // sidereal time. Original single-MJD routines approximate these scales.
                BOOST_CHECK_SMALL( polarDefault( i ) - fixture.polarLibration[ i ], 0.01 );
                BOOST_CHECK_SMALL( oceanDefault( i ) - fixture.ocean[ i ] * 1e6, 0.2 );
                BOOST_CHECK_SMALL( summedDefault( i ) - ( fixture.ocean[ i ] * 1e6 + fixture.polarLibration[ i ] ), 0.21 );

                // Alternative ortho-weight formulation has a quantified model/table
                // difference. Tight matching-table checks above are the primary oracle.
                BOOST_CHECK_SMALL( oceanDefault( i ) - fixture.ortho[ i ], 2.0 );
            }
            BOOST_CHECK_SMALL( oceanUt1Exact - fixture.ocean[ 2 ], 1e-10 );
            BOOST_CHECK_SMALL( oceanUt1.getCorrections( fixture.ttSeconds ) - fixture.ortho[ 2 ] * 1e-6, 1.2e-7 );
            BOOST_CHECK_SMALL( ut1Libration.getCorrections( fixture.ttSeconds ) - acceptedUt1, 2e-9 );
            BOOST_CHECK_SMALL( oceanUt1.getCorrections( fixture.ttSeconds ) - fixture.ocean[ 2 ], 3e-8 );
            BOOST_CHECK_SMALL( totalUt1->getCorrections( fixture.ttSeconds ) - ( fixture.ocean[ 2 ] + acceptedUt1 ), 3.2e-8 );
        }
    }
}

BOOST_AUTO_TEST_CASE( test_positive_amplitude_threshold_preserves_scalar_and_vector_terms )
{
    const std::string data = getEarthOrientationDataFilesPath( ) + "/";
    const auto ut1 = readAmplitudesAndFundamentalArgumentMultipliers(
            data + "utcLibrationAmplitudes.txt", data + "utcLibrationFundamentalArgumentMultipliers.txt", 1.0 );
    BOOST_REQUIRE_EQUAL( ut1.first.rows( ), 1 );
    BOOST_REQUIRE_EQUAL( ut1.first.cols( ), 2 );
    BOOST_REQUIRE_EQUAL( ut1.second.cols( ), 6 );
    BOOST_CHECK_EQUAL( ut1.first( 0, 0 ), 1.75 );
    BOOST_CHECK_EQUAL( ut1.second( 0, 3 ), -2.0 );
    ShortPeriodEarthOrientationCorrectionCalculator< double > scalar(
            1e-6, 1.0, { data + "utcLibrationAmplitudes.txt" }, { data + "utcLibrationFundamentalArgumentMultipliers.txt" } );
    Eigen::Vector6d zeroArguments = Eigen::Vector6d::Zero( );
    BOOST_CHECK_SMALL( scalar.getCorrectionsFromFundamentalArgument( zeroArguments ) + 1.01e-6, 1e-20 );
    const auto polar = readAmplitudesAndFundamentalArgumentMultipliers(
            data + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt",
            data + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt",
            20.0 );
    BOOST_REQUIRE_EQUAL( polar.first.rows( ), 1 );
    BOOST_REQUIRE_EQUAL( polar.first.cols( ), 4 );
    BOOST_REQUIRE_EQUAL( polar.second.cols( ), 6 );
    BOOST_CHECK_EQUAL( polar.first( 0, 0 ), 14.3 );
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > vector(
            1.0,
            20.0,
            { data + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
            { data + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" } );
    BOOST_CHECK_SMALL( vector.getCorrectionsFromFundamentalArgument( zeroArguments ).x( ) + 8.2, 1e-14 );
    BOOST_CHECK_SMALL( vector.getCorrectionsFromFundamentalArgument( zeroArguments ).y( ) - 14.3, 1e-14 );
    ShortPeriodEarthOrientationCorrectionCalculator< double > noScalar(
            1.0, 1e6, { data + "utcLibrationAmplitudes.txt" }, { data + "utcLibrationFundamentalArgumentMultipliers.txt" } );
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > noVector(
            1.0,
            1e6,
            { data + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
            { data + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" } );
    BOOST_CHECK_EQUAL( noScalar.getCorrectionsFromFundamentalArgument( zeroArguments ), 0.0 );
    BOOST_CHECK_EQUAL( noVector.getCorrectionsFromFundamentalArgument( zeroArguments ).norm( ), 0.0 );
    const auto empty = filterRawDataForAmplitudes( Eigen::MatrixXd( 0, 2 ), Eigen::MatrixXd( 0, 6 ), 1.0 );
    BOOST_CHECK_EQUAL( empty.first.rows( ), 0 );
    BOOST_CHECK_EQUAL( empty.second.cols( ), 6 );
}

BOOST_AUTO_TEST_CASE( test_daily_pole_and_cpo_are_applied_once_without_long_period_or_fcn_additions )
{
    const std::string data = getEarthOrientationDataFilesPath( ) + "/";
    const auto libration = readAmplitudesAndFundamentalArgumentMultipliers(
            data + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt",
            data + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt",
            0.0 );
    BOOST_REQUIRE_EQUAL( libration.second.rows( ), 10 );  // PMSDNUT2 IBAND=1, no long periods/secular drift.
    for( int i = 0; i < libration.second.rows( ); ++i )
    {
        BOOST_CHECK_EQUAL( libration.second( i, 0 ), 1.0 );
    }
    const auto ocean = readAmplitudesAndFundamentalArgumentMultipliers(
            data + "polarMotionOceanTidesAmplitudes.txt", data + "polarMotionOceanTidesFundamentalArgumentMultipliers.txt", 0.0 );
    BOOST_REQUIRE_EQUAL( ocean.second.rows( ), 71 );
    for( int i = 0; i < ocean.second.rows( ); ++i )
    {
        BOOST_CHECK( ocean.second( i, 0 ) == 1.0 || ocean.second( i, 0 ) == 2.0 );
    }
    const auto shortPeriods = getDefaultPolarMotionCorrectionCalculator( );
    const double uas = convertArcSecondsToRadians< double >( 1e-6 );
    for( const auto& fixture : short_period_reference::fixtures )
    {
        const double utc = ( fixture.mjdUtc - 51544.5 ) * 86400.0;
        const Time tai = sofa_interface::convertUTCtoTAI( Time( utc ) );
        const Eigen::Vector2d daily( 0.12 * 1e6 * uas, -0.34 * 1e6 * uas );

        // Synthetic observed CPO includes a 430-day oscillation representing an
        // FCN signal already in observations. It must enter X/Y exactly once.
        const Eigen::Vector2d cpo( 200.0 * uas * std::sin( 2.0 * mathematical_constants::PI * fixture.mjdUtc / 430.0 ), -80.0 * uas );

        // Daily records with constant pole and CPO, on the UTC days around the fixture epoch.
        input_output::EarthOrientationData data;
        for( int day = -1; day <= 2; ++day )
        {
            input_output::EarthOrientationRecord record;
            record.utcModifiedJulianDate = std::floor( fixture.mjdUtc ) + day;
            record.polarMotion = { daily.x( ), daily.y( ) };
            record.ut1MinusUtc = 0.0;
            record.celestialPoleOffsets = std::array< double, 2 >{ cpo.x( ), cpo.y( ) };
            data.records.emplace( record.utcModifiedJulianDate, record );
        }
        const auto eopData = std::make_shared< EarthOrientationDataEvaluator >( data );

        PolarMotionCalculator polarModel( eopData, shortPeriods );
        const auto actual = polarModel.getPositionOfCipInItrsFromTai(
                tai, sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst( fixture.ttSeconds ) );
        for( int i = 0; i < 2; ++i )
        {
            BOOST_CHECK_SMALL( actual( i ) - daily( i ) - uas * ( fixture.polarLibration[ i ] + fixture.ocean[ i ] * 1e6 ), 0.21 * uas );
        }

        PrecessionNutationCalculator precession( iau_2006, eopData );
        double x, y, s;
        iauXys06a( JULIAN_DAY_ON_J2000, fixture.ttSeconds / 86400.0, &x, &y, &s );
        const auto corrected = precession.getPositionOfCipInGcrsFromTai( fixture.ttSeconds, tai );
        BOOST_CHECK_SMALL( corrected.x( ) - x - cpo.x( ), 1e-17 );
        BOOST_CHECK_SMALL( corrected.y( ) - y - cpo.y( ), 1e-17 );
        BOOST_CHECK_SMALL( corrected.z( ) - s, 1e-17 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
