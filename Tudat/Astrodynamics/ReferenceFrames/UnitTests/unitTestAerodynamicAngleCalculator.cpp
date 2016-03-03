#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{
namespace unit_tests
{

using namespace unit_conversions;
using namespace reference_frames;

BOOST_AUTO_TEST_SUITE( test_aerodynamic_angle_calculator )

void testAerodynamicAngleCalculation(
        const basic_mathematics::Vector6d& testState,
        double testHeadingAngle,
        double testFlightPathAngle,
        double testLatitude,
        double testLongitude,
        double angleOfAttack,
        double angleOfSideslip,
        double bankAngle )
{
    AerodynamicAngleCalculator aerodynamicAngleCalculator(
                boost::lambda::constant( testState ),
                boost::lambda::constant( angleOfAttack ),
                boost::lambda::constant( angleOfSideslip),
                boost::lambda::constant( bankAngle ), 1 );

    aerodynamicAngleCalculator.update( );


    BOOST_CHECK_SMALL(
                std::fabs( ( aerodynamicAngleCalculator.getAerodynamicAngle( heading_angle ) -
                             testHeadingAngle ) ),2.0E-15 );
    BOOST_CHECK_SMALL(
                std::fabs( ( aerodynamicAngleCalculator.getAerodynamicAngle( flight_path_angle ) -
                             testFlightPathAngle ) ),2.0E-15 );
    BOOST_CHECK_SMALL(
                std::fabs( ( aerodynamicAngleCalculator.getAerodynamicAngle( latitude_angle ) -
                             testLatitude ) ),2.0E-15 );
    BOOST_CHECK_SMALL(
                std::fabs( ( aerodynamicAngleCalculator.getAerodynamicAngle( longitude_angle ) -
                             testLongitude ) ),2.0E-15 );

    Eigen::Matrix3d aerodynamicToBodyFrameMatrix =
            aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                aerodynamic_frame, body_frame ).toRotationMatrix( );
    Eigen::Matrix3d testAerodynamicToBodyFrameMatrix =
            getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
                angleOfAttack, angleOfSideslip );

    Eigen::Matrix3d trajectoryToAerodynamicFrameMatrix =
            aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                trajectory_frame, aerodynamic_frame ).toRotationMatrix( );
    Eigen::Matrix3d testTrajectoryToAerodynamicFrameMatrix =
            getTrajectoryToAerodynamicFrameTransformationMatrix(
                bankAngle );

    Eigen::Matrix3d verticalToTrajectoryFrameMatrix =
            aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                vertical_frame, trajectory_frame ).toRotationMatrix( );
    Eigen::Matrix3d testVerticalToTrajectoryFrameMatrix =
            getLocalVerticalFrameToTrajectoryTransformationMatrix(
                testFlightPathAngle, testHeadingAngle );

    Eigen::Matrix3d corotatingToVerticalFrameMatrix =
            aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                corotating_frame, vertical_frame ).toRotationMatrix( );
    Eigen::Matrix3d testCorotatingToVerticalFrameMatrix =
            getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(
                testLongitude, testLatitude );

    for( unsigned int l = 0; l < 3; l++ )
    {
        for( unsigned int m = 0; m < 3; m++ )
        {
            BOOST_CHECK_SMALL( std::fabs( aerodynamicToBodyFrameMatrix( l, m ) - testAerodynamicToBodyFrameMatrix( l, m ) ), 2.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( trajectoryToAerodynamicFrameMatrix( l, m ) - testTrajectoryToAerodynamicFrameMatrix( l, m ) ), 2.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( verticalToTrajectoryFrameMatrix( l, m ) - testVerticalToTrajectoryFrameMatrix( l, m ) ), 2.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( corotatingToVerticalFrameMatrix( l, m ) - testCorotatingToVerticalFrameMatrix( l, m ) ), 2.0E-15 );

        }

    }

    for( unsigned int i = 0; i < 5; i++ )
    {
        for( unsigned int j = 0; j < 5; j++ )
        {
            Eigen::Matrix3d directRotationMatrix =
                    aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                        static_cast< AerodynamicsReferenceFrames >( i ),
                        static_cast< AerodynamicsReferenceFrames >( j ) ).toRotationMatrix( );
            for( unsigned int k = 0; k < 5; k++ )
            {
                Eigen::Matrix3d indirectRotationMatrix =
                        aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                            static_cast< AerodynamicsReferenceFrames >( k ),
                            static_cast< AerodynamicsReferenceFrames >( j ) ).toRotationMatrix( ) *
                        aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                            static_cast< AerodynamicsReferenceFrames >( i ),
                            static_cast< AerodynamicsReferenceFrames >( k ) ).toRotationMatrix( );
                for( unsigned int l = 0; l < 3; l++ )
                {
                    for( unsigned int m = 0; m < 3; m++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( directRotationMatrix( l, m ) - indirectRotationMatrix( l, m ) ), 2.0E-15 );
                    }

                }
            }
        }
    }

}

// Test inertial to rotating planetocentric frame transformations.
BOOST_AUTO_TEST_CASE( testAerodynamicAngleCalculator )
{
    {
        basic_mathematics::Vector6d testState;
        testState<<-1656517.23153109, -5790058.28764025, -2440584.88186829,
                6526.30784888051, -2661.34558272018, 2377.09572383163;

        double testHeadingAngle = 1.229357188236127;
        double testFlightPathAngle = -0.024894033070522;
        double testLatitude = -0.385027359562548;
        double testLongitude = -1.849449608688977;

        double angleOfAttack = 1.232;
        double angleOfSideslip = -0.00322;
        double bankAngle = 2.323432;\

        testAerodynamicAngleCalculation( testState, testHeadingAngle, testFlightPathAngle, testLatitude,
                                         testLongitude, angleOfAttack, angleOfSideslip, bankAngle );
    }

    {
        basic_mathematics::Vector6d testState;
        testState<<0.0, 6498098.09700000, 0.0, 0.0, 0.0, 7.438147520000000e+03;

        double testHeadingAngle = 0.0;
        double testFlightPathAngle = 0.0;
        double testLatitude = 0.0;
        double testLongitude = mathematical_constants::PI / 2.0;

        double angleOfAttack = 1.232;
        double angleOfSideslip = -0.00322;
        double bankAngle = 2.323432;

        testAerodynamicAngleCalculation( testState, testHeadingAngle, testFlightPathAngle, testLatitude,
                                         testLongitude, angleOfAttack, angleOfSideslip, bankAngle );
    }

    {
        basic_mathematics::Vector6d testState;
        testState<<0.0, 0.0, 6.498098097000000e3, -7.438147520000000e3, 0.0, 0.0;

        double testHeadingAngle = 0.0;
        double testFlightPathAngle = 0.0;
        double testLatitude = mathematical_constants::PI / 2.0;
        double testLongitude = 0.0;

        double angleOfAttack = 1.232;
        double angleOfSideslip = -0.00322;
        double bankAngle = 2.323432;

        testAerodynamicAngleCalculation( testState, testHeadingAngle, testFlightPathAngle, testLatitude,
                                         testLongitude, angleOfAttack, angleOfSideslip, bankAngle );
    }



}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

