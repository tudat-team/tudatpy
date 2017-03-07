/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *
 *    References
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

    Eigen::Vector6d doodsonArguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( testEphemerisTime );

    Eigen::Vector2d totalPolarMotionFromTime = testPolarMotionCalculator->getPositionOfCipInItrs( testEphemerisTime, testUtc );
    Eigen::Vector2d totalPolarMotionFromArguments = testPolarMotionCalculator->getPositionOfCipInItrs( doodsonArguments, testUtc );

    BOOST_CHECK_EQUAL( totalPolarMotionFromTime.x( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).x( ) );
    BOOST_CHECK_EQUAL( totalPolarMotionFromTime.y( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).y( ) );

    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.x( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).x( ) );
    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.y( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( testEphemerisTime ) ).y( ) );

    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.x( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( doodsonArguments ) ).x( ) );
    BOOST_CHECK_EQUAL( totalPolarMotionFromArguments.y( ), ( dailyPolarMotionValueInterpolator->interpolate( testUtc ) +
                       shortPeriodPolarMotionCalculator->getCorrections( doodsonArguments ) ).y( ) );




}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





