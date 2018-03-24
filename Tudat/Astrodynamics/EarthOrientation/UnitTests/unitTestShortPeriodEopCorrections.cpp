/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/EarthOrientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"

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
BOOST_AUTO_TEST_CASE( testShortPeriodLibrationalPolarMotion)
{
    // Create polar motion correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionCalculator(
               convertArcSecondsToRadians< double >( 1.0E-6 ), 0.0,
    { getEarthOrientationDataFilesPath( ) + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
    { getEarthOrientationDataFilesPath( ) + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" } );

    //  Define test time
    double testMjd = 54335.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute polar motion correction
    double microAsToRadians =  mathematical_constants::PI / ( 180.0 * 1.0E6 * 3600.0 );
    Eigen::Vector2d polarMotionCorrections = polarMotionCalculator.getCorrections( testEphemerisTime ) / microAsToRadians;

    // Compare against IERS reference code. Difference between IERS code and this code occurs since the reference uses a slightly
    // different implementation for computation. Difference (1 micro arc seconds) is well below observable threshold
    BOOST_CHECK_SMALL( polarMotionCorrections( 0 ) - 24.65518398386097942, 1.0 );
    BOOST_CHECK_SMALL( polarMotionCorrections( 1 ) + 14.11070254891893327, 1.0 );

}

//! Test short-periodic polar motion libration variations by comparing to test output of ORTHO_EOP.f file
BOOST_AUTO_TEST_CASE( testShortPeriodOceanTidePolarMotion)
{
    // Create polar motion correction ocean tide corrections
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionCalculator(
               convertArcSecondsToRadians< double >( 1.0E-6 ), 0.0,
    { getEarthOrientationDataFilesPath( ) + "polarMotionOceanTidesAmplitudes.txt", },
    { getEarthOrientationDataFilesPath( ) + "polarMotionOceanTidesFundamentalArgumentMultipliers.txt" } );

    //  Define test time
    double testMjd = 47100.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute polar motion correction
    double microAsToRadians =  mathematical_constants::PI / ( 180.0 * 1.0E6 * 3600.0 );
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
               convertArcSecondsToRadians< double >( 1.0E-6 ), 0.0,
    { getEarthOrientationDataFilesPath( ) + "polarMotionOceanTidesAmplitudes.txt" },
    { getEarthOrientationDataFilesPath( ) + "polarMotionOceanTidesFundamentalArgumentMultipliers.txt" } );

    // Create polar motion correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionLibrationCorrectionCalculator(
               convertArcSecondsToRadians< double >( 1.0E-6 ), 0.0,
    { getEarthOrientationDataFilesPath( ) + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
    { getEarthOrientationDataFilesPath( ) + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" } );

    // Create polar motion correction combined corrections
    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionTotalCorrectionCalculator(
               convertArcSecondsToRadians< double >( 1.0E-6 ), 0.0,
    { getEarthOrientationDataFilesPath( ) + "polarMotionOceanTidesAmplitudes.txt",
      getEarthOrientationDataFilesPath( ) + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
    { getEarthOrientationDataFilesPath( ) + "polarMotionOceanTidesFundamentalArgumentMultipliers.txt",
      getEarthOrientationDataFilesPath( ) + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" } );

    //  Define test time
    double testMjd = 47100.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute corrections and check validity of combination
    Eigen::Vector2d polarMotionCorrectionOceanTides = polarMotionOceanTideCorrectionCalculator.getCorrections(
                testEphemerisTime );
    Eigen::Vector2d polarMotionCorrectionLibration = polarMotionLibrationCorrectionCalculator.getCorrections(
                testEphemerisTime );
    Eigen::Vector2d polarMotionCorrectionTotal = polarMotionTotalCorrectionCalculator.getCorrections( testEphemerisTime );


    BOOST_CHECK_SMALL( std::fabs( polarMotionCorrectionTotal( 0 ) -
                                  ( polarMotionCorrectionLibration( 0 ) + polarMotionCorrectionOceanTides( 0 ) ) ), 1.0E-24 );
    BOOST_CHECK_SMALL( std::fabs( polarMotionCorrectionTotal( 1 ) -
                                  ( polarMotionCorrectionLibration( 1 ) + polarMotionCorrectionOceanTides( 1 ) ) ), 1.0E-24 );
}

//! Test short-periodic ut1-utc variations by comparing to test output of UTLIBR.f file
BOOST_AUTO_TEST_CASE( testShortPeriodLibrationalUt1)
{
    // Create UT1 correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator < double > ut1CorrectionCalculator(
                1.0E-6, 0.0,
    { getEarthOrientationDataFilesPath( ) + "utcLibrationAmplitudes.txt" },
    { getEarthOrientationDataFilesPath( ) + "utcLibrationFundamentalArgumentMultipliers.txt" } );

    //  Define test time
    double testMjd = 44239.1 ;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay, JULIAN_DAY_ON_J2000 );
    double testEphemerisTime = convertUTCtoTT( testUtc );

    // Compute UT1 correction
    double ut1Correction = ut1CorrectionCalculator.getCorrections( testEphemerisTime );

    // Compare against IERS reference code. Difference between IERS code and this code occurs since the reference uses a slightly
    // different implementation for computation. Difference (10 ns) is well below observable threshold for Earth rotation
    BOOST_CHECK_SMALL( ut1Correction - 2.441143834386761746E-6, 1.0E-8 );

    //  Define second test time
    testMjd = 55227.4 ;
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
    ShortPeriodEarthOrientationCorrectionCalculator < double > ut1CorrectionCalculator =
            ShortPeriodEarthOrientationCorrectionCalculator < double >(
                1.0E-6, 0.0,
    { getEarthOrientationDataFilesPath( ) + "utcOceanTidesAmplitudes.txt", },
    { getEarthOrientationDataFilesPath( ) + "utcOceanTidesFundamentalArgumentMultipliers.txt" } );

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
    ShortPeriodEarthOrientationCorrectionCalculator < double > ut1OceanTideCorrectionCalculator =
            ShortPeriodEarthOrientationCorrectionCalculator < double >(
                1.0E-6, 0.0,
    { getEarthOrientationDataFilesPath( ) + "utcOceanTidesAmplitudes.txt" },
    { getEarthOrientationDataFilesPath( ) + "utcOceanTidesFundamentalArgumentMultipliers.txt" } );

    // Create UT1 correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator < double > ut1LibrationCorrectionCalculator(
                1.0E-6, 0.0,
    { getEarthOrientationDataFilesPath( ) + "utcLibrationAmplitudes.txt" },
    { getEarthOrientationDataFilesPath( ) + "utcLibrationFundamentalArgumentMultipliers.txt" } );

    // Create UT1 correction combined corrections
    ShortPeriodEarthOrientationCorrectionCalculator < double > ut1TotalCorrectionCalculator(
                1.0E-6, 0.0,
    { getEarthOrientationDataFilesPath( ) + "utcLibrationAmplitudes.txt",
                getEarthOrientationDataFilesPath( ) + "utcOceanTidesAmplitudes.txt" },
    { getEarthOrientationDataFilesPath( ) + "utcLibrationFundamentalArgumentMultipliers.txt",
                getEarthOrientationDataFilesPath( ) + "utcOceanTidesFundamentalArgumentMultipliers.txt" } );

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

BOOST_AUTO_TEST_SUITE_END( )

}

}



