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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <boost/test/unit_test.hpp>


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
    std::shared_ptr< EarthOrientationAnglesCalculator > standardEarthRotationModel =
            createStandardEarthOrientationCalculator( );
    std::shared_ptr< PolarMotionCalculator > standardPolarMotionCalculator =
            standardEarthRotationModel->getPolarMotionCalculator( );

    // Get constituent polar motion calculation objects.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > >
        dailyPolarMotionValueInterpolator = standardPolarMotionCalculator->getDailyIersValueInterpolator( );
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator =
            getDefaultPolarMotionCorrectionCalculator( );

    // Define test time
    double testEphemerisTime = 1.0E8;
    double testUtc = sofa_interface::convertTTtoUTC( testEphemerisTime );

    // Compute fundamental arguments
    Eigen::Vector6d fundamentalArguments = sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst( testEphemerisTime );

    // Compute polar motion from both interfaces (time and arguments)
    Eigen::Vector2d totalPolarMotionFromTime =
            standardPolarMotionCalculator->getPositionOfCipInItrs( testEphemerisTime, testUtc );
    Eigen::Vector2d totalPolarMotionFromArguments =
            standardPolarMotionCalculator->getPositionOfCipInItrs( fundamentalArguments, testUtc );

    // Check combinations fo calculations
    BOOST_CHECK_EQUAL( totalPolarMotionFromTime.x( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).x( ) );
    BOOST_CHECK_EQUAL( totalPolarMotionFromTime.y( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).y( ) );

    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.x( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).x( ) );
    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.y( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).y( ) );

    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.x( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrectionsFromFundamentalArgument( fundamentalArguments ) ).x( ) );
    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.y( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrectionsFromFundamentalArgument( fundamentalArguments ) ).y( ) );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





