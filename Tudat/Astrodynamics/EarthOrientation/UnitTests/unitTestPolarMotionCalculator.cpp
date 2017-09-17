/*    Copyright (c) 2010-2017, Delft University of Technology
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
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/EarthOrientation/polarMotionCalculator.h"
#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::earth_orientation;

BOOST_AUTO_TEST_SUITE( test_polar_motion_calculator )

// Check correct combination of data from the two contributions of calculation,
// correctness of these contributions is calculated in other unit tests.
BOOST_AUTO_TEST_CASE( testPolarMotionCalculator )
{
    boost::shared_ptr< EarthOrientationAnglesCalculator > standardEarthRotationModel =
            createStandardEarthOrientationCalculator( );

    boost::shared_ptr< PolarMotionCalculator > standardPolarMotionCalculator =
            standardEarthRotationModel->getPolarMotionCalculator( );

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > >
        dailyPolarMotionValueInterpolator = standardPolarMotionCalculator->getDailyIersValueInterpolator( );

    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator =
            getDefaultPolarMotionCorrectionCalculator( );

    boost::shared_ptr< PolarMotionCalculator > testPolarMotionCalculator =
            boost::make_shared< PolarMotionCalculator >( dailyPolarMotionValueInterpolator, shortPeriodPolarMotionCalculator );

    double testEphemerisTime = 0.0;

    double testUtc = sofa_interface::convertTTtoUTC( testEphemerisTime );

    Eigen::Vector6d fundamentalArguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( testEphemerisTime );

    Eigen::Vector2d totalPolarMotionFromTime = testPolarMotionCalculator->getPositionOfCipInItrs( testEphemerisTime, testUtc );
    Eigen::Vector2d totalPolarMotionFromArguments = testPolarMotionCalculator->getPositionOfCipInItrs( fundamentalArguments, testUtc );

    BOOST_CHECK_EQUAL( totalPolarMotionFromTime.x( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).x( ) );
    BOOST_CHECK_EQUAL( totalPolarMotionFromTime.y( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).y( ) );

    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.x( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).x( ) );
    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.y( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).y( ) );

    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.x( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( fundamentalArguments ) ).x( ) );
    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.y( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( fundamentalArguments ) ).y( ) );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





