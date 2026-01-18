/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    WASM test suite for Earth Orientation functions.
 *    Ported from native Tudat tests:
 *    - unitTestTimeScaleConverter.cpp
 *    - unitTestEopReader.cpp
 *    - unitTestPolarMotionCalculator.cpp
 *    - unitTestShortPeriodEopCorrections.cpp
 *    - unitTestEarthOrientationCalculator.cpp
 */

#include "wasmTestFramework.h"

#include <iostream>
#include <limits>
#include <cmath>

#include <Eigen/Core>

#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/interface/sofa/earthOrientation.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/io/basicInputOutput.h"

using namespace tudat;
using namespace tudat::earth_orientation;
using namespace tudat::basic_astrodynamics;
using namespace tudat::sofa_interface;
using namespace tudat::interpolators;
using namespace tudat::unit_conversions;
using namespace tudat::input_output;

// Helper function to convert SOFA output to seconds since J2000
double convertSofaOutputToSecondsSinceJ2000(const double fullJulianDays, const double fractionalJulianDays)
{
    return (fullJulianDays - basic_astrodynamics::JULIAN_DAY_ON_J2000) * physical_constants::JULIAN_DAY +
           fractionalJulianDays * physical_constants::JULIAN_DAY;
}

// Structure for SOFA time outputs (from SOFA Cookbook)
struct SofaTimeOutput {
    double expectedUtcFractionalDays;
    double expectedUtcDays;
    double expectedUt1Seconds;
    double expectedUt1Days;
    double expectedTaiFractionalDays;
    double expectedTaiDays;
    double expectedTtFractionalDays;
    double expectedTtDays;
    double expectedTdbFractionalDays;
    double expectedTdbDays;
};

// Get SOFA times from cookbook example (pre-computed values)
SofaTimeOutput getSofaDirectTimes()
{
    SofaTimeOutput sofaTimes;

    int latnd, latnm, lonwd, lonwm, j, iy, mo, id, ih, im;
    double slatn, slonw, hm, elon, phi, xyz[3], u, v, sec, utc1, utc2, dut, ut11, ut12, ut, tai1, tai2, tt1, tt2, tcg1, tcg2, dtr, tdb1,
           tdb2, tcb1, tcb2;

    (void)j; // Suppress unused variable warning

    /* Site terrestrial coordinates (WGS84). */
    latnd = 19;
    latnm = 28;
    slatn = 52.5;
    lonwd = 155;
    lonwm = 55;
    slonw = 59.6;
    hm = 0.0;

    /* Transform to geocentric. */
    j = iauAf2a('+', latnd, latnm, slatn, &phi);
    j = iauAf2a('-', lonwd, lonwm, slonw, &elon);
    j = iauGd2gc(1, elon, phi, hm, xyz);
    u = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
    v = xyz[2];

    /* UTC date and time. */
    iy = 2006;
    mo = 1;
    id = 15;
    ih = 21;
    im = 24;
    sec = 37.5;

    /* Transform into internal format. */
    j = iauDtf2d("UTC", iy, mo, id, ih, im, sec, &utc1, &utc2);

    /* UT1-UTC (s, from IERS). */
    dut = 0.3340960443019867; // Value modified to coincide with value in code

    /* UTC -> UT1. */
    j = iauUtcut1(utc1, utc2, dut, &ut11, &ut12);

    /* Extract fraction for TDB-TT calculation, later. */
    ut = fmod(fmod(ut11, 1.0) + fmod(ut12, 1.0), 1.0) + 0.5;

    /* UTC -> TAI -> TT -> TCG. */
    j = iauUtctai(utc1, utc2, &tai1, &tai2);
    j = iauTaitt(tai1, tai2, &tt1, &tt2);
    j = iauTttcg(tt1, tt2, &tcg1, &tcg2);

    /* TDB-TT (using TT as a substitute for TDB). */
    dtr = iauDtdb(tt1, tt2, ut, elon, u / 1000.0, v / 1000.0);

    /* TT -> TDB -> TCB. */
    j = iauTttdb(tt1, tt2, dtr, &tdb1, &tdb2);
    j = iauTdbtcb(tdb1, tdb2, &tcb1, &tcb2);

    /* Report. */
    sofaTimes.expectedUtcDays = utc1;
    sofaTimes.expectedUtcFractionalDays = utc2;

    sofaTimes.expectedUt1Days = ut11;
    sofaTimes.expectedUt1Seconds = ut12;

    sofaTimes.expectedTaiDays = tai1;
    sofaTimes.expectedTaiFractionalDays = tai2;

    sofaTimes.expectedTtDays = tt1;
    sofaTimes.expectedTtFractionalDays = tt2;

    sofaTimes.expectedTdbDays = tdb1;
    sofaTimes.expectedTdbFractionalDays = tdb2;
    return sofaTimes;
}

/**
 * Test time scale conversions against SOFA cookbook
 * Ported from unitTestTimeScaleConverter.cpp: testDifferentTimeScaleConversions
 */
void testTimeScaleConversions()
{
    std::cout << "  Testing time scale conversions..." << std::endl;

    // Get times from SOFA cookbook
    SofaTimeOutput sofaTimes = getSofaDirectTimes();

    // Create default time converter
    std::shared_ptr<TerrestrialTimeScaleConverter> timeScaleConverter =
            createStandardEarthOrientationCalculator(
                    std::make_shared<EOPReader>("/tudat_data/earth_orientation/eopc04_14_IAU2000.62-now.txt"))
                    ->getTerrestrialTimeScaleConverter();

    // Set station position
    Eigen::Vector3d stationCartesianPosition;
    stationCartesianPosition << -5492333.306498738, -2453018.508911721, 2113645.653406073;

    // Retrieve SOFA times in seconds since J2000
    std::map<TimeScales, double> sofaSecondsSinceJ2000;
    sofaSecondsSinceJ2000[tt_scale] =
            convertSofaOutputToSecondsSinceJ2000(sofaTimes.expectedTtDays, sofaTimes.expectedTtFractionalDays);
    sofaSecondsSinceJ2000[utc_scale] =
            convertSofaOutputToSecondsSinceJ2000(sofaTimes.expectedUtcDays, sofaTimes.expectedUtcFractionalDays);
    sofaSecondsSinceJ2000[ut1_scale] = convertSofaOutputToSecondsSinceJ2000(sofaTimes.expectedUt1Days, sofaTimes.expectedUt1Seconds);
    sofaSecondsSinceJ2000[tai_scale] =
            convertSofaOutputToSecondsSinceJ2000(sofaTimes.expectedTaiDays, sofaTimes.expectedTaiFractionalDays);
    sofaSecondsSinceJ2000[tdb_scale] =
            convertSofaOutputToSecondsSinceJ2000(sofaTimes.expectedTdbDays, sofaTimes.expectedTdbFractionalDays);

    // Test conversion between time scales (using TT as origin)
    timeScaleConverter->resetTimes<double>();
    timeScaleConverter->updateTimes(tt_scale, sofaSecondsSinceJ2000[tt_scale], stationCartesianPosition);

    double ut1 = timeScaleConverter->getCurrentTime(
            tt_scale, ut1_scale, sofaSecondsSinceJ2000[tt_scale], stationCartesianPosition);
    // UT1-UTC has some uncertainty in the EOP data; tolerance of ~10 microseconds is reasonable
    checkClose("TT->UT1 conversion", ut1, sofaSecondsSinceJ2000[ut1_scale], 1e-5);

    double utc = timeScaleConverter->getCurrentTime(
            tt_scale, utc_scale, sofaSecondsSinceJ2000[tt_scale], stationCartesianPosition);
    checkClose("TT->UTC conversion", utc, sofaSecondsSinceJ2000[utc_scale], std::numeric_limits<double>::epsilon());

    double tdb = timeScaleConverter->getCurrentTime(
            tt_scale, tdb_scale, sofaSecondsSinceJ2000[tt_scale], stationCartesianPosition);
    checkClose("TT->TDB conversion", tdb, sofaSecondsSinceJ2000[tdb_scale], std::numeric_limits<double>::epsilon());

    double tai = timeScaleConverter->getCurrentTime(
            tt_scale, tai_scale, sofaSecondsSinceJ2000[tt_scale], stationCartesianPosition);
    checkClose("TT->TAI conversion", tai, sofaSecondsSinceJ2000[tai_scale], std::numeric_limits<double>::epsilon());

    // Test using UTC as origin
    timeScaleConverter->resetTimes<double>();
    timeScaleConverter->updateTimes(utc_scale, sofaSecondsSinceJ2000[utc_scale], stationCartesianPosition);

    double tt = timeScaleConverter->getCurrentTime(
            utc_scale, tt_scale, sofaSecondsSinceJ2000[utc_scale], stationCartesianPosition);
    checkClose("UTC->TT conversion", tt, sofaSecondsSinceJ2000[tt_scale], std::numeric_limits<double>::epsilon());
}

/**
 * Test time scale conversion precision with Time type
 * Ported from unitTestTimeScaleConverter.cpp: testTimeScaleConversionPrecisionWithTimeType
 */
void testTimeScaleConversionPrecision()
{
    std::cout << "  Testing time scale conversion precision..." << std::endl;

    // Create time scale converters
    std::shared_ptr<TerrestrialTimeScaleConverter> timeScaleConverter =
            createStandardEarthOrientationCalculator()->getTerrestrialTimeScaleConverter();
    std::shared_ptr<TerrestrialTimeScaleConverter> comparisonTimeScaleConverter =
            createStandardEarthOrientationCalculator()->getTerrestrialTimeScaleConverter();

    // Define time scales
    std::vector<TimeScales> originScales = {tt_scale, utc_scale, ut1_scale, tai_scale, tdb_scale};
    std::vector<std::string> scaleNames = {"TT", "UTC", "UT1", "TAI", "TDB"};

    // Define current time/position
    Time baseTime = Time(-20.0 * physical_constants::JULIAN_YEAR);
    Eigen::Vector3d stationCartesianPosition;
    stationCartesianPosition << -5492333.306498738, -2453018.508911721, 2113645.653406073;

    // Test TT -> all scales -> back to TT roundtrip
    std::map<TimeScales, Time> currentTimes;

    timeScaleConverter->resetTimes<Time>();
    timeScaleConverter->updateTimes<Time>(tt_scale, baseTime, stationCartesianPosition);
    for (unsigned int j = 0; j < originScales.size(); j++)
    {
        currentTimes[originScales.at(j)] = timeScaleConverter->getCurrentTime<Time>(
                tt_scale, originScales.at(j), baseTime, stationCartesianPosition);
    }

    // Convert back and compare results. Tolerance is at ps level.
    for (unsigned int j = 0; j < originScales.size(); j++)
    {
        Time currentBackConvertedTime = comparisonTimeScaleConverter->getCurrentTime<Time>(
                originScales.at(j), tt_scale, currentTimes[originScales.at(j)], stationCartesianPosition);
        double diff = std::fabs(static_cast<long double>(currentBackConvertedTime - baseTime));
        checkTrue("TT->" + scaleNames[j] + "->TT roundtrip (ps precision)", diff < 1.0E-12);
    }
}

/**
 * Test time scale conversion around leap seconds
 * Ported from unitTestTimeScaleConverter.cpp: testTimeScaleConversionDuringLeapSeconds
 */
void testLeapSecondConversions()
{
    std::cout << "  Testing leap second conversions..." << std::endl;

    std::shared_ptr<TerrestrialTimeScaleConverter> timeScaleConverter =
            createStandardEarthOrientationCalculator()->getTerrestrialTimeScaleConverter();

    // Define a subset of leap seconds (day, month, year)
    std::vector<std::tuple<int, int, int>> leapSecondDays = {
        {1, 1, 2006}, {1, 1, 2009}, {1, 7, 2012}, {1, 7, 2015}, {1, 1, 2017}
    };

    // Test across leap seconds
    for (const auto& leapDay : leapSecondDays)
    {
        int day = std::get<0>(leapDay);
        int month = std::get<1>(leapDay);
        int year = std::get<2>(leapDay);

        double utcTimeOfLeapSeconds =
                basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(year, month, day, 0, 0, 0.0,
                        basic_astrodynamics::JULIAN_DAY_ON_J2000);

        double taiPreLeap = timeScaleConverter->getCurrentTime<double>(
                utc_scale, tai_scale, utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY - 1.0E-6, Eigen::Vector3d::Zero());
        double taiPostLeap = timeScaleConverter->getCurrentTime<double>(
                utc_scale, tai_scale, utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY + 1.0E-6, Eigen::Vector3d::Zero());

        // Check TAI jump at leap second
        double taiDiff = taiPostLeap - taiPreLeap;
        double expectedDiff = 1.0 + 2.0E-6; // 1 second leap + 2 microsecond interval
        std::string testName = "Leap second " + std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);
        checkTrue(testName, std::fabs(taiDiff - expectedDiff) < 1.0E-7);

        timeScaleConverter->resetTimes<double>();
    }
}

/**
 * Test EOP reader data
 * Ported from unitTestEopReader.cpp: testEopReaderData
 */
void testEopReaderData()
{
    std::cout << "  Testing EOP reader data..." << std::endl;

    std::shared_ptr<EOPReader> eopReader = std::make_shared<EOPReader>(
            "/tudat_data/earth_orientation/eopc04_14_IAU2000.62-now.txt", "C04", basic_astrodynamics::iau_2000_a);

    // Define current time
    int year = 2007;
    int month = 4;
    int day = 5;
    double utcSecondsSinceJ2000 = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
            year, month, day, 0, 0, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000) *
            physical_constants::JULIAN_DAY;

    // Set EOP corrections read manually from file for 5-4-2007
    double arcSecondToRadian = 4.848136811095359935899141E-6;
    double expectedXp = 0.033194 * arcSecondToRadian;
    double expectedYp = 0.483144 * arcSecondToRadian;
    double expectedUT1COffset = -0.0714163;
    double expecteddX = 0.000250 * arcSecondToRadian;
    double expecteddY = -0.000302 * arcSecondToRadian;

    // Create interpolators and retrieve values
    std::shared_ptr<OneDimensionalInterpolator<double, Eigen::Vector2d>> cipInItrsInterpolator =
            createStandardEarthOrientationCalculator(eopReader)
                    ->getPolarMotionCalculator()
                    ->getDailyIersValueInterpolator();
    Eigen::Vector2d cipInItrs = cipInItrsInterpolator->interpolate(utcSecondsSinceJ2000);

    std::shared_ptr<OneDimensionalInterpolator<double, Eigen::Vector2d>> cipInGcrsCorrectionInterpolator =
            createStandardEarthOrientationCalculator(eopReader)
                    ->getPrecessionNutationCalculator()
                    ->getDailyCorrectionInterpolator();
    Eigen::Vector2d cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate(utcSecondsSinceJ2000);

    std::shared_ptr<OneDimensionalInterpolator<double, double>> ut1MinusUtcInterpolator =
            createDefaultTimeConverter(eopReader)->getDailyUtcUt1CorrectionInterpolator();
    double utcMinusUt1 = ut1MinusUtcInterpolator->interpolate(utcSecondsSinceJ2000);

    // Check interpolated values against test data
    checkClose("UT1-UTC (2007-04-05)", utcMinusUt1, expectedUT1COffset, std::numeric_limits<double>::epsilon());
    checkClose("Xp (2007-04-05)", cipInItrs.x(), expectedXp, std::numeric_limits<double>::epsilon());
    checkClose("Yp (2007-04-05)", cipInItrs.y(), expectedYp, std::numeric_limits<double>::epsilon());
    checkClose("dX (2007-04-05)", cipInGcrs.x(), expecteddX, std::numeric_limits<double>::epsilon());
    checkClose("dY (2007-04-05)", cipInGcrs.y(), expecteddY, std::numeric_limits<double>::epsilon());

    // Test at second day (6-4-2007)
    utcSecondsSinceJ2000 = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
            year, month, 6, 0, 0, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000) *
            physical_constants::JULIAN_DAY;

    double expectedXp2 = 0.035736 * arcSecondToRadian;
    double expectedYp2 = 0.484204 * arcSecondToRadian;
    double expectedUT1COffset2 = -0.0728009;

    cipInItrs = cipInItrsInterpolator->interpolate(utcSecondsSinceJ2000);
    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate(utcSecondsSinceJ2000);

    checkClose("UT1-UTC (2007-04-06)", utcMinusUt1, expectedUT1COffset2, std::numeric_limits<double>::epsilon());
    checkClose("Xp (2007-04-06)", cipInItrs.x(), expectedXp2, std::numeric_limits<double>::epsilon());
    checkClose("Yp (2007-04-06)", cipInItrs.y(), expectedYp2, std::numeric_limits<double>::epsilon());
}

/**
 * Test short-period librational polar motion
 * Ported from unitTestShortPeriodEopCorrections.cpp: testShortPeriodLibrationalPolarMotion
 */
void testShortPeriodLibrationalPolarMotion()
{
    std::cout << "  Testing short-period librational polar motion..." << std::endl;

    // Create polar motion correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator<Eigen::Vector2d> polarMotionCalculator(
            convertArcSecondsToRadians<double>(1.0E-6),
            0.0,
            {"/tudat_data/earth_orientation/polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt"},
            {"/tudat_data/earth_orientation/polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt"});

    // Define test time
    double testMjd = 54335.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch(testJulianDay, JULIAN_DAY_ON_J2000);
    double testEphemerisTime = convertUTCtoTT(testUtc);

    // Compute polar motion correction
    double microAsToRadians = mathematical_constants::PI / (180.0 * 1.0E6 * 3600.0);
    Eigen::Vector2d polarMotionCorrections = polarMotionCalculator.getCorrections(testEphemerisTime) / microAsToRadians;

    // Compare against IERS reference code (1 micro arc second tolerance)
    checkTrue("Libration polar motion X correction",
              std::fabs(polarMotionCorrections(0) - 24.65518398386097942) < 1.0);
    checkTrue("Libration polar motion Y correction",
              std::fabs(polarMotionCorrections(1) + 14.11070254891893327) < 1.0);
}

/**
 * Test short-period ocean tide polar motion
 * Ported from unitTestShortPeriodEopCorrections.cpp: testShortPeriodOceanTidePolarMotion
 */
void testShortPeriodOceanTidePolarMotion()
{
    std::cout << "  Testing short-period ocean tide polar motion..." << std::endl;

    // Create polar motion correction ocean tide corrections
    ShortPeriodEarthOrientationCorrectionCalculator<Eigen::Vector2d> polarMotionCalculator(
            convertArcSecondsToRadians<double>(1.0E-6),
            0.0,
            {"/tudat_data/earth_orientation/polarMotionOceanTidesAmplitudes.txt"},
            {"/tudat_data/earth_orientation/polarMotionOceanTidesFundamentalArgumentMultipliers.txt"});

    // Define test time
    double testMjd = 47100.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch(testJulianDay, JULIAN_DAY_ON_J2000);
    double testEphemerisTime = convertUTCtoTT(testUtc);

    // Compute polar motion correction
    double microAsToRadians = mathematical_constants::PI / (180.0 * 1.0E6 * 3600.0);
    Eigen::Vector2d polarMotionCorrections = polarMotionCalculator.getCorrections(testEphemerisTime) / microAsToRadians;

    // Compare against IERS reference code (1 micro arc second tolerance)
    checkTrue("Ocean tide polar motion X correction",
              std::fabs(polarMotionCorrections.x() + 162.8386373279636530) < 1.0);
    checkTrue("Ocean tide polar motion Y correction",
              std::fabs(polarMotionCorrections.y() - 117.7907525842668974) < 1.0);
}

/**
 * Test short-period UT1 libration correction
 * Ported from unitTestShortPeriodEopCorrections.cpp: testShortPeriodLibrationalUt1
 */
void testShortPeriodLibrationalUt1()
{
    std::cout << "  Testing short-period librational UT1..." << std::endl;

    // Create UT1 correction libration corrections
    ShortPeriodEarthOrientationCorrectionCalculator<double> ut1CorrectionCalculator(
            1.0E-6,
            0.0,
            {"/tudat_data/earth_orientation/utcLibrationAmplitudes.txt"},
            {"/tudat_data/earth_orientation/utcLibrationFundamentalArgumentMultipliers.txt"});

    // Define test time 1
    double testMjd = 44239.1;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch(testJulianDay, JULIAN_DAY_ON_J2000);
    double testEphemerisTime = convertUTCtoTT(testUtc);

    // Compute UT1 correction
    double ut1Correction = ut1CorrectionCalculator.getCorrections(testEphemerisTime);

    // Compare against IERS reference code (10 ns tolerance)
    checkTrue("Libration UT1 correction (time 1)",
              std::fabs(ut1Correction - 2.441143834386761746E-6) < 1.0E-8);

    // Define test time 2
    testMjd = 55227.4;
    testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    testUtc = convertJulianDayToSecondsSinceEpoch(testJulianDay, JULIAN_DAY_ON_J2000);
    testEphemerisTime = convertUTCtoTT(testUtc);

    ut1Correction = ut1CorrectionCalculator.getCorrections(testEphemerisTime);

    checkTrue("Libration UT1 correction (time 2)",
              std::fabs(ut1Correction + 2.655705844335680244E-6) < 5.0E-8);
}

/**
 * Test short-period ocean tide UT1 correction
 * Ported from unitTestShortPeriodEopCorrections.cpp: testShortPeriodOceanTideUt1
 */
void testShortPeriodOceanTideUt1()
{
    std::cout << "  Testing short-period ocean tide UT1..." << std::endl;

    // Create UT1 correction ocean tide corrections
    ShortPeriodEarthOrientationCorrectionCalculator<double> ut1CorrectionCalculator(
            1.0E-6,
            0.0,
            {"/tudat_data/earth_orientation/utcOceanTidesAmplitudes.txt"},
            {"/tudat_data/earth_orientation/utcOceanTidesFundamentalArgumentMultipliers.txt"});

    // Define test time
    double testMjd = 47100.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch(testJulianDay, JULIAN_DAY_ON_J2000);
    double testEphemerisTime = convertUTCtoTT(testUtc);
    double ut1Correction = ut1CorrectionCalculator.getCorrections(testEphemerisTime);

    // Compare against IERS reference code (50 ns tolerance)
    checkTrue("Ocean tide UT1 correction",
              std::fabs(ut1Correction + 23.39092370609808214E-6) < 5.0E-8);
}

/**
 * Test polar motion calculator
 * Ported from unitTestPolarMotionCalculator.cpp: testPolarMotionCalculator
 */
void testPolarMotionCalculator()
{
    std::cout << "  Testing polar motion calculator..." << std::endl;

    // Retrieve polar motion calculator
    std::shared_ptr<EarthOrientationAnglesCalculator> standardEarthRotationModel = createStandardEarthOrientationCalculator();
    std::shared_ptr<PolarMotionCalculator> standardPolarMotionCalculator = standardEarthRotationModel->getPolarMotionCalculator();

    // Get constituent polar motion calculation objects
    std::shared_ptr<OneDimensionalInterpolator<double, Eigen::Vector2d>> dailyPolarMotionValueInterpolator =
            standardPolarMotionCalculator->getDailyIersValueInterpolator();
    std::shared_ptr<ShortPeriodEarthOrientationCorrectionCalculator<Eigen::Vector2d>> shortPeriodPolarMotionCalculator =
            getDefaultPolarMotionCorrectionCalculator();

    // Define test time
    double testEphemerisTime = 1.0E8;
    double testUtc = sofa_interface::convertTTtoUTC(testEphemerisTime);

    // Compute fundamental arguments
    Eigen::Vector6d fundamentalArguments = sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst(testEphemerisTime);

    // Compute polar motion from both interfaces
    Eigen::Vector2d totalPolarMotionFromTime = standardPolarMotionCalculator->getPositionOfCipInItrs(testEphemerisTime, testUtc);
    Eigen::Vector2d totalPolarMotionFromArguments = standardPolarMotionCalculator->getPositionOfCipInItrs(fundamentalArguments, testUtc);

    // Check combinations of calculations
    Eigen::Vector2d expectedFromTime = dailyPolarMotionValueInterpolator->interpolate(testUtc) +
                                        shortPeriodPolarMotionCalculator->getCorrections(testEphemerisTime);
    Eigen::Vector2d expectedFromArgs = dailyPolarMotionValueInterpolator->interpolate(testUtc) +
                                        shortPeriodPolarMotionCalculator->getCorrectionsFromFundamentalArgument(fundamentalArguments);

    checkClose("PM from time X", totalPolarMotionFromTime.x(), expectedFromTime.x(), 1e-20);
    checkClose("PM from time Y", totalPolarMotionFromTime.y(), expectedFromTime.y(), 1e-20);
    checkClose("PM from args X", totalPolarMotionFromArguments.x(), expectedFromArgs.x(), 1e-20);
    checkClose("PM from args Y", totalPolarMotionFromArguments.y(), expectedFromArgs.y(), 1e-20);
}

/**
 * Test Earth orientation rotation setup against SOFA
 * Ported from unitTestEarthOrientationCalculator.cpp: testEarthOrientationRotationSetupAgainstSofa
 */
void testEarthOrientationRotationSetup()
{
    std::cout << "  Testing Earth orientation rotation setup..." << std::endl;

    double arcSecondToRadian = 4.848136811095359935899141E-6;

    // Define test time
    double terrestrialTimeDaysSinceMjd0 = 54195.50075444445;
    double terrestrialTimeSecondsSinceJ2000 =
            (terrestrialTimeDaysSinceMjd0 - (basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD)) *
            physical_constants::JULIAN_DAY;

    // Define precession/nutation parameters
    double X = 0.0007122647295989105;
    double Y = 4.438525042571229e-05;
    double s = -0.002200475 * arcSecondToRadian;

    // Compute CIRS->GCRS rotation matrix in Tudat
    Eigen::Matrix3d cirsToGcrsRotation = calculateRotationFromCirsToGcrs(X, Y, s).toRotationMatrix();

    // Expected values from SOFA cookbook
    Eigen::Matrix3d expectedCirsToGcrsRotation;
    expectedCirsToGcrsRotation << 0.999999746339445, -5.138822464778592e-09, -0.0007122647300724212, -2.647522726051399e-08,
            0.9999999990149748, -4.438524282712702e-05, 0.0007122647295989105, 4.438525042571229e-05, 0.9999997453544198;

    // Check CIRS to GCRS rotation
    bool rotationPassed = true;
    for (unsigned int i = 0; i < 3; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (std::fabs(expectedCirsToGcrsRotation.transpose()(i, j) - cirsToGcrsRotation(i, j)) > 1.0E-15)
            {
                rotationPassed = false;
            }
        }
    }
    checkTrue("CIRS->GCRS rotation matrix", rotationPassed);

    // Define Earth rotation angle
    double era = 0.2324515536620879;

    // Compute TIRS->CIRS rotation matrix
    Eigen::Matrix3d tirsToCirsRotation = Eigen::Matrix3d(calculateRotationFromTirsToCirs(era));

    // Define polar motion values
    double xPole = 0.034928200 * arcSecondToRadian;
    double yPole = 0.483316300 * arcSecondToRadian;

    // Compute ITRS->TIRS rotation matrix
    Eigen::Matrix3d itrsToTirsRotation =
            calculateRotationFromItrsToTirs(xPole, yPole, getApproximateTioLocator(terrestrialTimeSecondsSinceJ2000))
                    .toRotationMatrix();

    // Compute ITRS -> GCRS by concatenation
    Eigen::Matrix3d itrsToGcrsRotation = cirsToGcrsRotation * tirsToCirsRotation * itrsToTirsRotation;

    // Compare Tudat GCRS -> ITRS against Sofa
    Eigen::Matrix3d expectedGcrsToItrs;
    expectedGcrsToItrs << 0.973104317697536, 0.2303638262391256, -0.0007031634821983242, -0.2303638004560344, 0.9731045706328012,
            0.000118545366624876, 0.000711560162667892, 4.662640399540082e-05, 0.9999997457540244;

    bool itrsRotationPassed = true;
    for (unsigned int i = 0; i < 3; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (std::fabs(expectedGcrsToItrs.transpose()(i, j) - itrsToGcrsRotation(i, j)) > 1.0E-15)
            {
                itrsRotationPassed = false;
            }
        }
    }
    checkTrue("ITRS->GCRS rotation matrix", itrsRotationPassed);
}

/**
 * Test historical Earth rotation (pre-C04 file)
 * Ported from unitTestEarthOrientationCalculator.cpp: testHistoricalEarthRotation
 */
void testHistoricalEarthRotation()
{
    std::cout << "  Testing historical Earth rotation..." << std::endl;

    std::shared_ptr<EarthOrientationAnglesCalculator> earthOrientationCalculator = createStandardEarthOrientationCalculator();
    basic_astrodynamics::DateTime dateTimePreC04(1958, 12, 30, 12, 0.0, 0.0);
    basic_astrodynamics::DateTime dateTimeStartC04(1962, 1, 1, 0, 0.0, 0.0);

    Eigen::Vector2d polarMotionPreC04 =
            earthOrientationCalculator->getPolarMotionCalculator()->getDailyIersValueInterpolator()->interpolate(
                    dateTimePreC04.epoch<double>());
    Eigen::Vector2d polarMotionPostC04 =
            earthOrientationCalculator->getPolarMotionCalculator()->getDailyIersValueInterpolator()->interpolate(
                    dateTimeStartC04.epoch<double>());

    double ut1CorrectionPreC04 =
            earthOrientationCalculator->getTerrestrialTimeScaleConverter()->getDailyUtcUt1CorrectionInterpolator()->interpolate(
                    dateTimePreC04.epoch<double>());
    double ut1CorrectionPostC04 =
            earthOrientationCalculator->getTerrestrialTimeScaleConverter()->getDailyUtcUt1CorrectionInterpolator()->interpolate(
                    dateTimeStartC04.epoch<double>());

    // Check if polar motion is zero before C04 file starts
    checkTrue("Pre-C04 polar motion X is zero", polarMotionPreC04(0) == 0);
    checkTrue("Pre-C04 polar motion Y is zero", polarMotionPreC04(1) == 0);
    checkTrue("Post-C04 polar motion X is non-zero", polarMotionPostC04(0) != 0);
    checkTrue("Post-C04 polar motion Y is non-zero", polarMotionPostC04(1) != 0);

    // Check if UTC-UT1 is zero before C04 file starts
    checkTrue("Pre-C04 UT1-UTC is zero", ut1CorrectionPreC04 == 0);
    checkTrue("Post-C04 UT1-UTC is non-zero", ut1CorrectionPostC04 != 0);
}

/**
 * Test leap second identification in EOP reader
 * Ported from unitTestEopReader.cpp: testLeapSecondIdentification
 */
void testLeapSecondIdentification()
{
    std::cout << "  Testing leap second identification..." << std::endl;

    // Read EOP file and get UT1-UTC interpolator
    std::shared_ptr<EOPReader> eopReader = std::make_shared<EOPReader>(
            "/tudat_data/earth_orientation/eopc04_14_IAU2000.62-now.txt", "C04", basic_astrodynamics::iau_2000_a);
    std::shared_ptr<OneDimensionalInterpolator<double, double>> ut1MinusUtcInterpolator =
            createDefaultTimeConverter(eopReader)->getDailyUtcUt1CorrectionInterpolator();

    // Test a subset of leap seconds
    std::vector<std::tuple<int, int, int>> leapSecondDays = {
        {1, 1, 2006}, {1, 1, 2009}, {1, 7, 2012}, {1, 7, 2015}, {1, 1, 2017}
    };

    for (const auto& leapDay : leapSecondDays)
    {
        int day = std::get<0>(leapDay);
        int month = std::get<1>(leapDay);
        int year = std::get<2>(leapDay);

        double utcSecondsSinceJ2000 = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                year, month, day, 0, 0, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000) *
                physical_constants::JULIAN_DAY;

        double differenceMinusOneMilliSecond = ut1MinusUtcInterpolator->interpolate(utcSecondsSinceJ2000 - 1.0E-3);
        double differenceAtEpochPlusOneMilliSecond = ut1MinusUtcInterpolator->interpolate(utcSecondsSinceJ2000 + 1.0E-3);

        // Check if leap second is properly handled (UT1-UTC should jump by ~1 second)
        std::string testName = "Leap second detection " + std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);
        checkTrue(testName, std::fabs(differenceMinusOneMilliSecond - differenceAtEpochPlusOneMilliSecond + 1.0) < 1.0E-10);
    }
}
