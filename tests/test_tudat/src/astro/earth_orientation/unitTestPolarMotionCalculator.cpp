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

#include <limits>
#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/earth_orientation/polarMotionCalculator.h"
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::earth_orientation;

BOOST_AUTO_TEST_SUITE( test_polar_motion_calculator )

//! Check correct combination of data from the two contributions of polar motion: short-period terms and IERS daily values.
BOOST_AUTO_TEST_CASE( testPolarMotionCalculator )
{
    // Retrieve polar motion calculator
    std::shared_ptr< EarthOrientationAnglesCalculator > standardEarthRotationModel = createStandardEarthOrientationCalculator( );
    std::shared_ptr< PolarMotionCalculator > standardPolarMotionCalculator = standardEarthRotationModel->getPolarMotionCalculator( );

    // Get constituent polar motion calculation objects.
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator =
            getDefaultPolarMotionCorrectionCalculator( );

    // Define test time
    const double testEphemerisTime = 1.0E8;
    const auto converter = standardEarthRotationModel->getTerrestrialTimeScaleConverter( );
    const Time tai = converter->getTaiTime( basic_astrodynamics::tt_scale, testEphemerisTime );
    const Eigen::Vector2d dailyPolarMotion = converter->getEarthOrientationDataEvaluator( )->evaluate( tai ).polarMotion;

    // From TAI and given fundamental arguments: daily values plus short-period terms for these arguments.
    const Eigen::Vector6d fundamentalArguments =
            sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst( testEphemerisTime );
    const Eigen::Vector2d polarMotionFromArguments =
            standardPolarMotionCalculator->getPositionOfCipInItrsFromTai( tai, fundamentalArguments );
    const Eigen::Vector2d expectedFromArguments =
            dailyPolarMotion + shortPeriodPolarMotionCalculator->getCorrectionsFromFundamentalArgument( fundamentalArguments );
    BOOST_CHECK_EQUAL( polarMotionFromArguments.x( ), expectedFromArguments.x( ) );
    BOOST_CHECK_EQUAL( polarMotionFromArguments.y( ), expectedFromArguments.y( ) );

    // From TT: short-period terms use arguments computed from TT and the modelled UT1.
    const Eigen::Vector2d polarMotionFromTime = standardPolarMotionCalculator->getPositionOfCipInItrsAtTt( testEphemerisTime );
    const auto arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst(
            testEphemerisTime, testEphemerisTime, static_cast< double >( converter->getUt1FromTai( tai ) ) );
    const Eigen::Vector2d expectedFromTime =
            dailyPolarMotion + shortPeriodPolarMotionCalculator->getCorrectionsFromFundamentalArgument( arguments );
    BOOST_CHECK_EQUAL( polarMotionFromTime.x( ), expectedFromTime.x( ) );
    BOOST_CHECK_EQUAL( polarMotionFromTime.y( ), expectedFromTime.y( ) );
    BOOST_CHECK_GT( ( polarMotionFromTime - polarMotionFromArguments ).norm( ), 1e-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
