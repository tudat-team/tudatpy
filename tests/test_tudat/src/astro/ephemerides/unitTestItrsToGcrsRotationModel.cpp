/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/earth_orientation/sofaEarthOrientationCookbookExamples.h"
#include "tudat/interface/sofa/earthOrientation.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace earth_orientation;
using namespace basic_astrodynamics;
using namespace simulation_setup;

namespace
{

// SOFA Tools for Earth Attitude, C edition, revision 1.64, section 5.5:
// 2007-04-05 12:00:00 UTC, xp=0.0349282 arcsec, yp=0.4833163 arcsec,
// UT1-UTC=-0.072073685 s, dX=0.1750 mas, dY=-0.2259 mas.
// The calculator receives these values as constant daily EOP data, without
// short-period correction files.
std::shared_ptr< EarthOrientationAnglesCalculator > createCookbookAnglesCalculator( )
{
    const Eigen::Vector2d pole = unit_conversions::convertArcSecondsToRadians( Eigen::Vector2d( 0.0349282, 0.4833163 ) );
    const Eigen::Vector2d correction = unit_conversions::convertArcSecondsToRadians( Eigen::Vector2d( 0.0001750, -0.0002259 ) );

    // Daily records with these constant values, from 1996 to 2009, which covers all test epochs.
    input_output::EarthOrientationData data;
    for( double mjd = 50300.0; mjd <= 55100.0; mjd += 1.0 )
    {
        input_output::EarthOrientationRecord record;
        record.utcModifiedJulianDate = mjd;
        record.polarMotion = { pole.x( ), pole.y( ) };
        record.ut1MinusUtc = -0.072073685;
        record.celestialPoleOffsets = std::array< double, 2 >{ correction.x( ), correction.y( ) };
        data.records.emplace( mjd, record );
    }
    const auto eopData = std::make_shared< EarthOrientationDataEvaluator >( data );
    const auto zeroArguments = []( const double ) -> Eigen::Vector6d { return Eigen::Vector6d::Zero( ); };
    const auto zeroPolarCorrections = std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
            1.0, 0.0, std::vector< std::string >( ), std::vector< std::string >( ), zeroArguments );
    const auto zeroUt1Corrections = std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< double > >(
            1.0, 0.0, std::vector< std::string >( ), std::vector< std::string >( ), zeroArguments );
    return std::make_shared< EarthOrientationAnglesCalculator >(
            std::make_shared< PolarMotionCalculator >( eopData, zeroPolarCorrections ),
            std::make_shared< PrecessionNutationCalculator >( iau_2006, eopData ),
            std::make_shared< TerrestrialTimeScaleConverter >( eopData, zeroUt1Corrections ) );
}

std::map< std::string, Eigen::Matrix3d > getIndependentBaseFrameRotations( )
{
    // Time-independent rb fixture from the IAU SOFA 2023-10-11 validation
    // program t_sofa_c.c, t_bp06. Do not obtain the oracle from getFrameBias.
    const Eigen::Matrix3d gcrsToEme2000 = ( Eigen::Matrix3d( ) << 0.9999999999999942497,
                                            -0.7078368960971557145E-7,
                                            0.8056213977613185606E-7,
                                            0.7078368694637674333E-7,
                                            0.9999999999999969484,
                                            0.3305943742989134124E-7,
                                            -0.8056214211620056792E-7,
                                            -0.3305943172740586950E-7,
                                            0.9999999999999962084 )
                                                  .finished( );

    // NAIF's built-in ECLIPJ2000 frame: rotation about x by 84381.448 arcsec.
    // The fixed values below are independent of the production SPICE wrapper.
    const Eigen::Matrix3d j2000ToEclipj2000 = ( Eigen::Matrix3d( ) << 1.0,
                                                0.0,
                                                0.0,
                                                0.0,
                                                0.9174820620691818,
                                                0.3977771559319137,
                                                0.0,
                                                -0.3977771559319137,
                                                0.9174820620691818 )
                                                      .finished( );

    // J2000 (and its synonym EME2000) includes the frame bias; ECLIPJ2000 is defined with respect to J2000.
    return { { "GCRS", Eigen::Matrix3d::Identity( ) },
             { "ICRF", Eigen::Matrix3d::Identity( ) },
             { "J2000", gcrsToEme2000 },
             { "EME2000", gcrsToEme2000 },
             { "ECLIPJ2000", j2000ToEclipj2000 * gcrsToEme2000 } };
}

void checkRotationMatrix( const Eigen::Matrix3d& actual, const Eigen::Matrix3d& expected, const double tolerance )
{
    for( int row = 0; row < 3; ++row )
    {
        for( int column = 0; column < 3; ++column )
        {
            BOOST_CHECK_SMALL( actual( row, column ) - expected( row, column ), tolerance );
        }
    }
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_itrs_to_gcrs_rotation )

//! Test ITRS <-> GCRS rotation by compariong against Spice
BOOST_AUTO_TEST_CASE( test_ItrsToGcrsRotationAgainstSpice )
{
    spice_interface::loadSpiceKernelInTudat( paths::getSpiceKernelPath( ) + "/naif0012.tls" );
    spice_interface::loadSpiceKernelInTudat( paths::getSpiceKernelPath( ) + "/earth_latest_high_prec.bpc" );
    spice_interface::loadSpiceKernelInTudat( paths::getSpiceKernelPath( ) + "/earth_fixed.tf" );

    // Create rotation model
    // SPICE J2000 has ICRF axes, which are the GCRS axes of this model.
    std::shared_ptr< GcrsToItrsRotationModel > earthRotationModel =
            std::make_shared< GcrsToItrsRotationModel >( earth_orientation::createStandardEarthOrientationCalculator( ) );

    // Compare spice vs. Tudat for list of evaluation times
    std::vector< double > testTimes;
    testTimes.push_back( 1.0E8 );
    testTimes.push_back( 1.0E7 );
    testTimes.push_back( 0.0 );
    for( unsigned test = 0; test < testTimes.size( ); test++ )
    {
        Eigen::Matrix3d sofaRotation = earthRotationModel->getRotationToBaseFrame( testTimes.at( test ) ).toRotationMatrix( );
        Eigen::Matrix3d sofaRotationDerivative = earthRotationModel->getDerivativeOfRotationToBaseFrame( testTimes.at( test ) );
        Eigen::Matrix3d spiceRotation =
                spice_interface::computeRotationQuaternionBetweenFrames( "ITRF93", "J2000", testTimes.at( test ) ).toRotationMatrix( );
        Eigen::Matrix3d spiceRotationDerivative =
                spice_interface::computeRotationMatrixDerivativeBetweenFrames( "ITRF93", "J2000", testTimes.at( test ) );

        // Check whether Spice and Tudat give same result. Note that Spice model is not accurate up to IERS standards. Comparison
        // is done at 10 cm position difference on Earth surface (per component).
        double tolerance = 0.1 / 6378.0E3;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( sofaRotation( i, j ) - spiceRotation( i, j ), tolerance );
                BOOST_CHECK_SMALL( sofaRotationDerivative( i, j ) - spiceRotationDerivative( i, j ), 1.0E-11 );
            }
        }
    }
}

//! Test ITRS <-> GCRS rotation by compariong against Sofa
BOOST_AUTO_TEST_CASE( test_ItrsToGcrsRotationAgainstSofaCookbook )
{
    // Get UTC time for evaluation
    int year = 2007;
    int month = 4;
    int day = 5;
    int hour = 12;
    int minutes = 0;
    double seconds = 0.0;
    double sofaCookbookTime = convertCalendarDateToJulianDaysSinceEpoch( year, month, day, hour, minutes, seconds, JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Create Earth rotation model
    std::shared_ptr< GcrsToItrsRotationModel > earthRotationModelFromUtc =
            std::make_shared< GcrsToItrsRotationModel >( earth_orientation::createStandardEarthOrientationCalculator( ), utc_scale );

    // Test Tudat vs. Sofa implementations, with default Sofa EOP corrections (as defined in cookbook
    {
        Eigen::Matrix3d sofaCookbookResult = getSofaEarthOrientationExamples( 3 ).transpose( );
        Eigen::Matrix3d tudatResult = earthRotationModelFromUtc->getRotationToBaseFrame( sofaCookbookTime ).toRotationMatrix( );

        // Check sofa against Tudat result, small difference due to slightly different values of EOP corrections.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( sofaCookbookResult( i, j ) - tudatResult( i, j ), 2.0E-9 );
            }
        }
    }

    // Test Tudat vs. Sofa implementations with identical EOP corrections.
    {
        // Set current time in UTC and TT
        double interpolationUtc = sofaCookbookTime;
        double interpolationTt = earthRotationModelFromUtc->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( )->getCurrentTime(
                utc_scale, tt_scale, interpolationUtc );

        // Get EOP corrections
        const auto eop =
                earthRotationModelFromUtc->getAnglesCalculator( )->getPrecessionNutationCalculator( )->getEarthOrientationDataEvaluator( );
        const Eigen::Vector2d cpo = eop->getCelestialPoleOffsets( sofa_interface::convertUTCtoTAI( Time( interpolationUtc ) ) );
        const double Xcorrection = cpo.x( ), Ycorrection = cpo.y( );
        const Eigen::Vector2d polarMotion =
                earthRotationModelFromUtc->getAnglesCalculator( )->getPolarMotionCalculator( )->getPositionOfCipInItrsAtTt(
                        interpolationTt );
        const double xPolarMotion = polarMotion.x( ), yPolarMotion = polarMotion.y( );
        double ut1Correction = earthRotationModelFromUtc->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( )->getUt1Correction(
                utc_scale, Time( sofaCookbookTime ) );

        // Compute Sofa rotation matrix
        Eigen::Matrix3d sofaCookbookResult =
                getSofaEarthOrientationExamples( 3,
                                                 unit_conversions::convertRadiansToArcSeconds( Xcorrection ) * 1000.0,
                                                 unit_conversions::convertRadiansToArcSeconds( Ycorrection ) * 1000.0,
                                                 unit_conversions::convertRadiansToArcSeconds( xPolarMotion ),
                                                 unit_conversions::convertRadiansToArcSeconds( yPolarMotion ),
                                                 ut1Correction )
                        .transpose( );

        // Compute Tudat rotation matrix
        Eigen::Matrix3d tudatResult = earthRotationModelFromUtc->getRotationToBaseFrame( sofaCookbookTime ).toRotationMatrix( );

        // Check sofa against Tudat result, small difference due to rounding errors, in particular in Earth rotation angle
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                if( i < 2 && j < 2 )
                {
                    BOOST_CHECK_SMALL( sofaCookbookResult( i, j ) - tudatResult( i, j ), 1.0E-11 );
                }
                else
                {
                    BOOST_CHECK_SMALL( sofaCookbookResult( i, j ) - tudatResult( i, j ), 1.0E-14 );
                }
            }
        }

        // Compute Tudat rotation matrix with high-precision time input
        long double sofaCookbookExtendedTime =
                convertCalendarDateToJulianDaysSinceEpoch< long double >( year, month, day, hour, minutes, seconds, JULIAN_DAY_ON_J2000 ) *
                physical_constants::JULIAN_DAY_LONG;
        Eigen::Matrix3d tudatResultPrecise =
                earthRotationModelFromUtc->getRotationToBaseFrameFromExtendedTime( Time( sofaCookbookExtendedTime ) ).toRotationMatrix( );

        // Check sofa against Tudat result
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( sofaCookbookResult( i, j ) - tudatResultPrecise( i, j ), 1.0E-15 );
            }
        }
    }
}

//! Test each base frame against an independent absolute orientation reference.
BOOST_AUTO_TEST_CASE( test_InertialFramesAgainstFrozenSofaCookbook )
{
    const auto calculator = createCookbookAnglesCalculator( );
    const double utc =
            convertCalendarDateToJulianDaysSinceEpoch( 2007, 4, 5, 12, 0, 0.0, JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;
    const Time extendedUtc( convertCalendarDateToJulianDaysSinceEpoch< long double >( 2007, 4, 5, 12, 0, 0.0L, JULIAN_DAY_ON_J2000 ) *
                            physical_constants::JULIAN_DAY_LONG );

    // Published celestial-to-terrestrial matrix, SOFA cookbook section 5.5.
    // Transposition maps the independent reference from ITRS to GCRS axes.
    const Eigen::Matrix3d gcrsToItrs = ( Eigen::Matrix3d( ) << 0.973104317697535,
                                         0.230363826239128,
                                         -0.000703163482198,
                                         -0.230363800456037,
                                         0.973104570632801,
                                         0.000118545366625,
                                         0.000711560162668,
                                         0.000046626403995,
                                         0.999999745754024 )
                                               .finished( );

    for( const auto& frame : getIndependentBaseFrameRotations( ) )
    {
        BOOST_TEST_CONTEXT( "base frame = " << frame.first )
        {
            GcrsToItrsRotationModel model( calculator, utc_scale, frame.first );
            const Eigen::Matrix3d expected = frame.second * gcrsToItrs.transpose( );
            BOOST_CHECK_EQUAL( model.getBaseFrameOrientation( ), frame.first );
            checkRotationMatrix( model.getRotationToBaseFrame( utc ).toRotationMatrix( ), expected, 2.0E-11 );
            checkRotationMatrix( model.getRotationToTargetFrame( utc ).toRotationMatrix( ), expected.transpose( ), 2.0E-11 );
            checkRotationMatrix( model.getRotationToBaseFrameFromExtendedTime( extendedUtc ).toRotationMatrix( ), expected, 2.0E-14 );
            checkRotationMatrix(
                    model.getRotationToTargetFrameFromExtendedTime( extendedUtc ).toRotationMatrix( ), expected.transpose( ), 2.0E-14 );
        }
    }
    BOOST_CHECK_THROW( GcrsToItrsRotationModel( calculator, utc_scale, "UNKNOWN" ), std::runtime_error );
}

//! Test fixed frame mappings for rotations, inverse rotations, extended time and derivatives.
BOOST_AUTO_TEST_CASE( test_InertialFrameRotationAndDerivativeMapping )
{
    const auto calculator = createCookbookAnglesCalculator( );
    GcrsToItrsRotationModel gcrsModel( calculator, utc_scale, "GCRS" );
    const std::vector< double > testTimes = { -1.0E8, 0.0, 1.0E8, 3.0E8 };
    for( const double time : testTimes )
    {
        const Eigen::Matrix3d gcrsRotation = gcrsModel.getRotationToBaseFrame( time ).toRotationMatrix( );
        const Eigen::Matrix3d gcrsDerivative = gcrsModel.getDerivativeOfRotationToBaseFrame( time );
        const Time extendedTime( static_cast< long double >( time ) + 0.123456789L );
        const Eigen::Matrix3d extendedGcrsRotation = gcrsModel.getRotationToBaseFrameFromExtendedTime( extendedTime ).toRotationMatrix( );
        for( const auto& frame : getIndependentBaseFrameRotations( ) )
        {
            BOOST_TEST_CONTEXT( "base frame = " << frame.first << ", UTC seconds since J2000 = " << time )
            {
                GcrsToItrsRotationModel model( calculator, utc_scale, frame.first );
                const Eigen::Matrix3d expectedRotation = frame.second * gcrsRotation;
                const Eigen::Matrix3d expectedDerivative = frame.second * gcrsDerivative;
                const Eigen::Matrix3d expectedExtendedRotation = frame.second * extendedGcrsRotation;
                checkRotationMatrix( model.getRotationToBaseFrame( time ).toRotationMatrix( ), expectedRotation, 2.0E-15 );
                checkRotationMatrix( model.getRotationToTargetFrame( time ).toRotationMatrix( ), expectedRotation.transpose( ), 2.0E-15 );
                checkRotationMatrix( model.getDerivativeOfRotationToBaseFrame( time ), expectedDerivative, 2.0E-19 );
                checkRotationMatrix( model.getDerivativeOfRotationToTargetFrame( time ), expectedDerivative.transpose( ), 2.0E-19 );
                checkRotationMatrix( model.getRotationToBaseFrameFromExtendedTime( extendedTime ).toRotationMatrix( ),
                                     expectedExtendedRotation,
                                     2.0E-15 );
                checkRotationMatrix( model.getRotationToTargetFrameFromExtendedTime( extendedTime ).toRotationMatrix( ),
                                     expectedExtendedRotation.transpose( ),
                                     2.0E-15 );
            }
        }
    }
}

//! Test rotations between IERS intermediate frames.
BOOST_AUTO_TEST_CASE( test_RotationBetweenIntermediateFrames )
{
    double testTime = 1.0E8;
    std::shared_ptr< GcrsToItrsRotationModel > earthRotationModel =
            std::make_shared< GcrsToItrsRotationModel >( createStandardEarthOrientationCalculator( ) );

    // The model evaluates the TIO locator at TT and the Earth rotation angle at the precise (Time) UT1.
    const EarthOrientationEvaluatedAnglesAndRates orientation =
            earthRotationModel->getAnglesCalculator( )->getRotationAnglesAndTimesFromItrsToGcrs( testTime, tdb_scale );
    Eigen::Vector5d rotationAngles = orientation.angles;

    Eigen::Matrix3d expectedItrsToTirs =
            calculateRotationFromItrsToTirs(
                    rotationAngles[ 3 ], rotationAngles[ 4 ], getApproximateTioLocator( static_cast< double >( orientation.tt ) ) )
                    .toRotationMatrix( );
    Eigen::Matrix3d expectedTirsToCirs =
            calculateRotationFromTirsToCirs( sofa_interface::calculateEarthRotationAngleTemplated< Time >( orientation.ut1 ) )
                    .toRotationMatrix( );
    Eigen::Matrix3d expectedCirsToIcrs =
            calculateRotationFromCirsToGcrs( rotationAngles[ 0 ], rotationAngles[ 1 ], rotationAngles[ 2 ] ).toRotationMatrix( );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            earthRotationModel
                    ->getRotationBetweenIntermediateFrames(
                            EarthOrientationIntermediateFrame::itrs, EarthOrientationIntermediateFrame::tirs, testTime )
                    .toRotationMatrix( ),
            expectedItrsToTirs,
            ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            earthRotationModel
                    ->getRotationBetweenIntermediateFrames(
                            EarthOrientationIntermediateFrame::tirs, EarthOrientationIntermediateFrame::cirs, testTime )
                    .toRotationMatrix( ),
            expectedTirsToCirs,
            ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            earthRotationModel
                    ->getRotationBetweenIntermediateFrames(
                            EarthOrientationIntermediateFrame::cirs, EarthOrientationIntermediateFrame::icrs, testTime )
                    .toRotationMatrix( ),
            expectedCirsToIcrs,
            ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            earthRotationModel
                    ->getRotationBetweenIntermediateFrames(
                            EarthOrientationIntermediateFrame::itrs, EarthOrientationIntermediateFrame::icrs, testTime )
                    .toRotationMatrix( ),
            earthRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( ),
            ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );

    Eigen::Matrix3d expectedTirsToItrs = expectedItrsToTirs.transpose( );
    Eigen::Matrix3d tirsToItrs =
            earthRotationModel
                    ->getRotationBetweenIntermediateFrames(
                            EarthOrientationIntermediateFrame::tirs, EarthOrientationIntermediateFrame::itrs, testTime )
                    .toRotationMatrix( );
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( tirsToItrs( i, j ) - expectedTirsToItrs( i, j ) ), 1.0E-15 );
        }
    }
}

//! Test ITRS <-> GCRS rotation by comparing results from different base frames.
BOOST_AUTO_TEST_CASE( test_ItrsToGcrsRotationConsistency )
{
    // 1. Load SPICE kernels
    spice_interface::loadStandardSpiceKernels( );

    // 2. Define test time
    double testTime = 1.0E8;

    // 3. Scenario 1: J2000 as global frame
    Eigen::Vector3d marsPositionInItrsFromJ2000;
    {
        // Define bodies to be created
        std::vector< std::string > bodiesToCreate = { "Earth", "Mars" };

        // Define body settings
        BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, "SSB", "J2000" );

        // Create system of bodies
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Set Earth rotation model with J2000 as base frame
        auto earthRotationModelJ2000 =
                std::make_shared< GcrsToItrsRotationModel >( createStandardEarthOrientationCalculator( ), tdb_scale, "J2000" );
        bodies.at( "Earth" )->setRotationalEphemeris( earthRotationModelJ2000 );

        // Calculate Mars position w.r.t. Earth in J2000
        Eigen::Vector3d marsPositionInJ2000 = bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( testTime ).segment< 3 >( 0 ) -
                bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( testTime ).segment< 3 >( 0 );

        // Get rotation from J2000 to ITRS
        Eigen::Quaterniond rotationJ2000ToItrs = earthRotationModelJ2000->getRotationToTargetFrame( testTime );

        // Transform Mars position to ITRS
        marsPositionInItrsFromJ2000 = rotationJ2000ToItrs.toRotationMatrix( ) * marsPositionInJ2000;
    }

    // 4. Scenario 2: ECLIPJ2000 as global frame
    Eigen::Vector3d marsPositionInItrsFromEclipj2000;
    {
        // Define bodies to be created
        std::vector< std::string > bodiesToCreate = { "Earth", "Mars" };

        // Define body settings
        BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, "SSB", "ECLIPJ2000" );

        // Create system of bodies
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Set Earth rotation model with ECLIPJ2000 as base frame
        auto earthRotationModelEclipj2000 =
                std::make_shared< GcrsToItrsRotationModel >( createStandardEarthOrientationCalculator( ), tdb_scale, "ECLIPJ2000" );
        bodies.at( "Earth" )->setRotationalEphemeris( earthRotationModelEclipj2000 );

        // Calculate Mars position w.r.t. Earth in ECLIPJ2000
        Eigen::Vector3d marsPositionInEclipj2000 = bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( testTime ).segment< 3 >( 0 ) -
                bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( testTime ).segment< 3 >( 0 );

        // Get rotation from ECLIPJ2000 to ITRS
        Eigen::Quaterniond rotationEclipj2000ToItrs = earthRotationModelEclipj2000->getRotationToTargetFrame( testTime );

        // Transform Mars position to ITRS
        marsPositionInItrsFromEclipj2000 = rotationEclipj2000ToItrs.toRotationMatrix( ) * marsPositionInEclipj2000;
    }

    // 5. Comparison in ITRS
    // Check that the two positions in ITRS are consistent
    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( marsPositionInItrsFromJ2000( i ), marsPositionInItrsFromEclipj2000( i ), 1.0E-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
