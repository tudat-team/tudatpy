#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <iostream>

#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/ground_stations/poleTideDeformation.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"

namespace tudat
{
namespace unit_tests
{

using namespace basic_astrodynamics;
using namespace earth_orientation;
using namespace interpolators;
using namespace coordinate_conversions;
using namespace ground_stations;

BOOST_AUTO_TEST_SUITE( test_pole_tide_site_displacement )

BOOST_AUTO_TEST_CASE( testPoleTideSiteDisplacement )
{

    // Retrieve polar motion calculator
    std::shared_ptr< EarthOrientationAnglesCalculator > standardEarthRotationModel =
        createStandardEarthOrientationCalculator( );
    std::shared_ptr< PolarMotionCalculator > standardPolarMotionCalculator =
        standardEarthRotationModel->getPolarMotionCalculator( );

    std::bind( static_cast< Eigen::Vector2d( OneDimensionalInterpolator< double, Eigen::Vector2d >::* )( const double )>
                 ( &OneDimensionalInterpolator< double, Eigen::Vector2d >::interpolate ),
               standardPolarMotionCalculator->getDailyIersValueInterpolator( ), std::placeholders::_1 );

    PoleTideDeformation poleTideDeformationModel = PoleTideDeformation( standardPolarMotionCalculator );

    double radiansToArcSeconds = unit_conversions::convertArcSecondsToRadians( 1.0 );

    for( int testSite = 0; testSite < 3; testSite++ )
    {
        std::shared_ptr<ground_stations::GroundStationState> nominalSiteState;
        if( testSite == 0 )
        {
            Eigen::Vector3d sphericalState = ( Eigen::Vector3d( ) << 6378.0E3, 0.0, 0.0 ).finished( );
            nominalSiteState = std::make_shared<ground_stations::GroundStationState>(
                sphericalState, spherical_position );
        }
        else if( testSite == 1 )
        {
                Eigen::Vector3d sphericalState = ( Eigen::Vector3d( ) << 6378.0E3, 0.0, mathematical_constants::PI / 2.0 ).finished( );
                nominalSiteState = std::make_shared<ground_stations::GroundStationState>(
                    sphericalState, spherical_position );
        }
        else if( testSite == 2 )
        {
            Eigen::Vector3d sphericalState = ( Eigen::Vector3d( ) << 6378.0E3, mathematical_constants::PI / 4.0, 0.0 ).finished( );
            nominalSiteState = std::make_shared<ground_stations::GroundStationState>(
                sphericalState, spherical_position );
        }

        for ( int i = 0; i < 10; i++ )
        {
            double currentTime = physical_constants::JULIAN_YEAR * static_cast< double >( i );
            Eigen::Vector2d offsets =
                unit_conversions::convertRadiansToDegrees(
                    standardEarthRotationModel->getPolarMotionCalculator( )->getPositionOfCipInItrs(
                        currentTime, sofa_interface::convertTTtoUTC( currentTime ) ) ) * 3600.0 * 1000.0 -
                    earth_orientation::getSecularPolePositionInMas( currentTime );
            offsets( 1 ) *= -1.0;

            // Check that results are in sensible range
            BOOST_CHECK_SMALL( std::fabs( offsets( 0 ) ), 300.0 );
            BOOST_CHECK_SMALL( std::fabs( offsets( 1 ) ), 300.0 );

            Eigen::Vector3d enuDisplacement = poleTideDeformationModel.calculateEnuDisplacement( currentTime, nominalSiteState );
            if( testSite == 0 )
            {
                double expectedNorthDisplacement = -9.0 * offsets( 0 ) * radiansToArcSeconds;
                double expectedEastDisplacement = -9.0 * offsets( 1 ) * radiansToArcSeconds;

                BOOST_CHECK_SMALL( std::fabs( expectedNorthDisplacement - enuDisplacement( 1 ) ), 1.0E-9 );
                BOOST_CHECK_SMALL( std::fabs( expectedEastDisplacement - enuDisplacement( 0 ) ), 1.0E-9 );
                BOOST_CHECK_SMALL( std::fabs( enuDisplacement( 2 ) ), 1.0E-9 );
            }
            else if( testSite == 1 )
            {
                double expectedNorthDisplacement = -9.0 * offsets( 1 ) * radiansToArcSeconds;
                double expectedEastDisplacement = 9.0 * offsets( 0 ) * radiansToArcSeconds;

                BOOST_CHECK_SMALL( std::fabs( expectedNorthDisplacement - enuDisplacement( 1 ) ), 1.0E-9 );
                BOOST_CHECK_SMALL( std::fabs( expectedEastDisplacement - enuDisplacement( 0 ) ), 1.0E-9 );
                BOOST_CHECK_SMALL( std::fabs( enuDisplacement( 2 ) ), 1.0E-9 );
            }
            else if( testSite == 2 )
            {
                double expectedNorthDisplacement = 0.0;
                double expectedEastDisplacement = -9.0 / std::sqrt( 2.0 ) * offsets( 1 ) * radiansToArcSeconds;
                double expectedUpDisplacement = -33.0 * offsets( 0 ) * radiansToArcSeconds;

                BOOST_CHECK_SMALL( std::fabs( expectedNorthDisplacement - enuDisplacement( 1 ) ), 1.0E-9 );
                BOOST_CHECK_SMALL( std::fabs( expectedEastDisplacement - enuDisplacement( 0 ) ), 1.0E-9 );
                BOOST_CHECK_SMALL( std::fabs( expectedUpDisplacement - enuDisplacement( 2 ) ), 1.0E-9 );
            }
        }
    }



}

BOOST_AUTO_TEST_SUITE_END( )

}

}
