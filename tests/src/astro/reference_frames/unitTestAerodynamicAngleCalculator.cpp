/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
namespace tudat
{
namespace unit_tests
{

using namespace unit_conversions;
using namespace reference_frames;

class ManualAerodynamicAngleInterface: public BodyFixedAerodynamicAngleInterface
{
public:
    ManualAerodynamicAngleInterface(
            const std::function< Eigen::Vector3d( const double ) > manualAngleFuncion ):
    BodyFixedAerodynamicAngleInterface( custom_body_fixed_angles ),manualAngleFuncion_( manualAngleFuncion ){ }

    virtual ~ManualAerodynamicAngleInterface( ){ }

    Eigen::Vector3d getAngles( const double time,
                               const Eigen::Matrix3d& trajectoryToInertialFrame )
    {
        return manualAngleFuncion_( time );
    }

private:

    std::function< Eigen::Vector3d( const double ) > manualAngleFuncion_;
};

BOOST_AUTO_TEST_SUITE( test_aerodynamic_angle_calculator )

//! Function to test the aerodynamic angle calculator.
/*!
 * Function to test the aerodynamic angle calculator from a current body-fixed Cartesian state, the
 *  current orientation angles, and the expected spherical 6-dimensional state.
 * \param testState Body-fixed cartesian state.
 * \param testHeadingAngle Heading angle corresponding to testState
 * \param testFlightPathAngle Flight path angle corresponding to testState
 * \param testLatitude Latitude angle corresponding to testState
 * \param testLongitude Longitude angle corresponding to testState
 * \param angleOfAttack Angle of attack of vehicle
 * \param angleOfSideslip Angle of sideslip of vehicle
 * \param bankAngle Bank angle of vehicle
 */
void testAerodynamicAngleCalculation(
        const Eigen::Vector6d& testState,
        double testHeadingAngle,
        double testFlightPathAngle,
        double testLatitude,
        double testLongitude,
        double angleOfAttack,
        double angleOfSideslip,
        double bankAngle )
{
    // Create angle calculator
    AerodynamicAngleCalculator aerodynamicAngleCalculator(
                [ & ]( ){ return testState; },
                [ & ]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); }, "", 1 );

    aerodynamicAngleCalculator.setBodyFixedAngleInterface(
                std::make_shared< ManualAerodynamicAngleInterface >(
                    [=](const double)
    {
        return( Eigen::Vector3d( ) << angleOfAttack, angleOfSideslip, bankAngle ).finished( );
    } ) );

    // Update angle calculator.
    aerodynamicAngleCalculator.update( 0.0, true );

    // Compare expected against computed angles.
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

    // Compute rotation matrices manually and from AerodynamicAngleCalculator
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

    // Compare rotation matrices.
    for( unsigned int l = 0; l < 3; l++ )
    {
        for( unsigned int m = 0; m < 3; m++ )
        {
            BOOST_CHECK_SMALL( std::fabs( aerodynamicToBodyFrameMatrix( l, m ) -
                                          testAerodynamicToBodyFrameMatrix( l, m ) ), 2.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( trajectoryToAerodynamicFrameMatrix( l, m ) -
                                          testTrajectoryToAerodynamicFrameMatrix( l, m ) ), 2.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( verticalToTrajectoryFrameMatrix( l, m ) -
                                          testVerticalToTrajectoryFrameMatrix( l, m ) ), 2.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( corotatingToVerticalFrameMatrix( l, m ) -
                                          testCorotatingToVerticalFrameMatrix( l, m ) ), 2.0E-15 );

        }

    }

    // Test rotation matrix obtained directly and from combination of two rotation matrices from
    // AerodynamicAngleCalculator
    for( unsigned int i = 0; i < 5; i++ )
    {
        for( unsigned int j = 0; j < 5; j++ )
        {
            // Calculate direct rotation matrix
            Eigen::Matrix3d directRotationMatrix =
                    aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                        static_cast< AerodynamicsReferenceFrames >( i ),
                        static_cast< AerodynamicsReferenceFrames >( j ) ).toRotationMatrix( );

            // Calculate each possible intermediate rotation matrix.
            for( unsigned int k = 0; k < 5; k++ )
            {
                Eigen::Matrix3d indirectRotationMatrix =
                        aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                            static_cast< AerodynamicsReferenceFrames >( k ),
                            static_cast< AerodynamicsReferenceFrames >( j ) ).toRotationMatrix( ) *
                        aerodynamicAngleCalculator.getRotationQuaternionBetweenFrames(
                            static_cast< AerodynamicsReferenceFrames >( i ),
                            static_cast< AerodynamicsReferenceFrames >( k ) ).toRotationMatrix( );

                // Compare matrices
                for( unsigned int l = 0; l < 3; l++ )
                {
                    for( unsigned int m = 0; m < 3; m++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( directRotationMatrix( l, m )
                                                      - indirectRotationMatrix( l, m ) ), 2.0E-15 );
                    }

                }
            }
        }
    }

}

//! Test inertial to rotating planetocentric frame transformations, using Matlab script of Erwin
//! Mooij to generate reference data.
BOOST_AUTO_TEST_CASE( testAerodynamicAngleCalculator )
{
    // Test case 1: arbitrary rotation
    {
        std::cout<<"case 1"<<std::endl;
        Eigen::Vector6d testState;
        testState << -1656517.23153109, -5790058.28764025, -2440584.88186829,
                6526.30784888051, -2661.34558272018, 2377.09572383163;

        double testHeadingAngle = 1.229357188236127;
        double testFlightPathAngle = -0.024894033070522;
        double testLatitude = -0.385027359562548;
        double testLongitude = -1.849449608688977;

        double angleOfAttack = 1.232;
        double angleOfSideslip = -0.00322;
        double bankAngle = 2.323432;\

        testAerodynamicAngleCalculation( testState, testHeadingAngle, testFlightPathAngle,
                                         testLatitude, testLongitude,
                                         angleOfAttack, angleOfSideslip, bankAngle );
    }

    // Test case 2: rotation with zero and half pi angles.
    {
        std::cout<<"case 2"<<std::endl;
        Eigen::Vector6d testState;
        testState << 0.0, 6498098.09700000, 0.0, 0.0, 0.0, 7.438147520000000e+03;

        double testHeadingAngle = 0.0;
        double testFlightPathAngle = 0.0;
        double testLatitude = 0.0;
        double testLongitude = mathematical_constants::PI / 2.0;

        double angleOfAttack = 1.232;
        double angleOfSideslip = -0.00322;
        double bankAngle = 2.323432;

        testAerodynamicAngleCalculation( testState, testHeadingAngle, testFlightPathAngle,
                                         testLatitude, testLongitude,
                                         angleOfAttack, angleOfSideslip, bankAngle );
    }

    // Test case 3: rotation with zero and half pi angles.
    {
        std::cout<<"case 3"<<std::endl;
        Eigen::Vector6d testState;
        testState << 0.0, 0.0, 6.498098097000000e3, -7.438147520000000e3, 0.0, 0.0;

        double testHeadingAngle = 0.0;
        double testFlightPathAngle = 0.0;
        double testLatitude = mathematical_constants::PI / 2.0;
        double testLongitude = 0.0;

        double angleOfAttack = 1.232;
        double angleOfSideslip = -0.00322;
        double bankAngle = 2.323432;

        testAerodynamicAngleCalculation( testState, testHeadingAngle, testFlightPathAngle,
                                         testLatitude, testLongitude,
                                         angleOfAttack, angleOfSideslip, bankAngle );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

