/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Musegaas, P. (2012). Optimization of Space Trajectories Including Multiple Gravity Assists
 *          and Deep Space Maneuvers. MSc Thesis, Delft University of Technology, Delft,
 *          The Netherlands.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <vector>

#include <boost/make_shared.hpp>
#include <memory>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <tudat/astro/basic_astro/physicalConstants.h>
#include <tudat/basics/testMacros.h>
#include <tudat/math/basic/mathematicalConstants.h>
#include "tudat/astro/basic_astro/unitConversions.h"
#include <tudat/astro/basic_astro/orbitalElementConversions.h>


#include <tudat/io/basicInputOutput.h>

#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/mission_segments/createTransferTrajectory.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodies.h"

namespace tudat
{
namespace unit_tests
{

// Additional namespaces to be used.
using namespace tudat::input_output;
using namespace tudat::simulation_setup;
using namespace mission_segments;


//! Test implementation of trajectory class
BOOST_AUTO_TEST_SUITE( test_trajectory )

BOOST_AUTO_TEST_CASE( testMGATrajectory_New )
{
    // Expected test result based on the ideal Cassini 1 trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 4930.72686847243;

    // Create environment
    tudat::simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Set transfer order
    std::vector< std::string > bodyOrder = {
        "Earth", "Venus", "Venus", "Earth", "Jupiter", "Saturn" };
    int numberOfNodes = bodyOrder.size( );

    std::shared_ptr< TransferTrajectory > transferTrajectory;
    double nominalDeltaV = TUDAT_NAN;
    double nominalCaptureDeltaV = TUDAT_NAN;

    for( unsigned int creationType = 0; creationType < 3; creationType++ )
    {
        // Create leg settings (all unpowered)
        std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
        std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;

        // Define minimum periapsis altitudes for flybys;
        std::map< std::string, double > minimumPeriapses;
        minimumPeriapses[ "Venus" ] = 6351800.0;
        minimumPeriapses[ "Earth" ] = 6678000.0;
        minimumPeriapses[ "Jupiter" ] =  600000000.0;
        minimumPeriapses[ "Saturn" ] = 65000000.0;

        if( creationType == 0 )
        {
            transferLegSettings.resize( numberOfNodes - 1 );
            transferLegSettings[ 0 ] = unpoweredLeg( );
            transferLegSettings[ 1 ] = unpoweredLeg( );
            transferLegSettings[ 2 ] = unpoweredLeg( );
            transferLegSettings[ 3 ] = unpoweredLeg( );
            transferLegSettings[ 4 ] = unpoweredLeg( );



            transferNodeSettings.resize( numberOfNodes );
            transferNodeSettings[ 0 ] = escapeAndDepartureNode( std::numeric_limits< double >::infinity( ), 0.0 );
            transferNodeSettings[ 1 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 1 ) ) );
            transferNodeSettings[ 2 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 2 ) ) );
            transferNodeSettings[ 3 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 3 ) ) );
            transferNodeSettings[ 4 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 4 ) ) );
            transferNodeSettings[ 5 ] = captureAndInsertionNode( 1.0895e8 / 0.02, 0.98 );
        }
        else if( creationType == 1 )
        {
            getMgaTransferTrajectorySettingsWithoutDsm(
                        transferLegSettings, transferNodeSettings, bodyOrder,
                        std::make_pair( std::numeric_limits< double >::infinity( ), 0.0 ),
                        std::make_pair( 1.0895e8 / 0.02, 0.98 ),
                        minimumPeriapses );
        }
        else if( creationType == 2 )
        {
            getMgaTransferTrajectorySettingsWithoutDsm(
                        transferLegSettings, transferNodeSettings, bodyOrder,
                        std::make_pair( std::numeric_limits< double >::infinity( ), 0.0 ),
                        std::make_pair( TUDAT_NAN, TUDAT_NAN ),
                        minimumPeriapses );

        }
        transferTrajectory = createTransferTrajectory(
                    bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun" );




        // Create list of node times
        double JD = physical_constants::JULIAN_DAY;
        std::vector< double > nodeTimes;
        // Subtract 0.5 to go from MJD2000 to J2000
        nodeTimes.push_back( ( -789.8117 - 0.5 ) * JD );
        nodeTimes.push_back( nodeTimes.at( 0 ) + 158.302027105278 * JD );
        nodeTimes.push_back( nodeTimes.at( 1 ) + 449.385873819743 * JD );
        nodeTimes.push_back( nodeTimes.at( 2 ) + 54.7489684339665 * JD );
        nodeTimes.push_back( nodeTimes.at( 3 ) + 1024.36205846918 * JD );
        nodeTimes.push_back( nodeTimes.at( 4 ) + 4552.30796805542 * JD );

        std::vector< Eigen::VectorXd > transferLegFreeParameters;
        for( int i = 0; i < numberOfNodes - 1; i++ )
        {
            transferLegFreeParameters.push_back( Eigen::VectorXd( 0 ) );
        }
        std::vector< Eigen::VectorXd > transferNodeFreeParameters;
        for( int i = 0; i < numberOfNodes; i++ )
        {
            transferNodeFreeParameters.push_back( Eigen::VectorXd( 0 ) );
        }
        if( creationType == 2 )
        {
            transferNodeFreeParameters[ numberOfNodes - 1 ] = ( Eigen::Vector3d( )<<65000000.0, 0.0, 0.0 ).finished( );
        }

        printTransferParameterDefinition( transferLegSettings, transferNodeSettings );

        transferTrajectory->evaluateTrajectory(
                    nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );

        if( creationType < 2 )
        {
            BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, transferTrajectory->getTotalDeltaV( ), 1.0E-3 );

            if( creationType == 0 )
            {
                nominalDeltaV = transferTrajectory->getTotalDeltaV( );
                nominalCaptureDeltaV = transferTrajectory->getNodeDeltaV( 5 );
            }
            else
            {
                BOOST_CHECK_CLOSE_FRACTION( nominalDeltaV, transferTrajectory->getTotalDeltaV( ), 1.0E-12 );
            }
        }
        else if( creationType == 2 )
        {

            BOOST_CHECK_CLOSE_FRACTION( nominalCaptureDeltaV, nominalDeltaV - transferTrajectory->getTotalDeltaV( ), 1.0E-12 );
        }

        std::vector< std::map< double, Eigen::Vector6d > > statesAlongTrajectoryPerLeg;
        transferTrajectory-> getStatesAlongTrajectoryPerLeg(
                    statesAlongTrajectoryPerLeg, 10 );

        double sunGravitationalParameter = bodies.at( "Sun" )->getGravitationalParameter( );
        for( unsigned int i = 0; i < statesAlongTrajectoryPerLeg.size( ); i++ )
        {
            double legStartTime = nodeTimes.at( i );
            double legEndTime = nodeTimes.at( i + 1 );

            std::map< double, Eigen::Vector6d > statesAlongSingleLeg = statesAlongTrajectoryPerLeg.at( i );

            // Check initial and final time on output list
            BOOST_CHECK_CLOSE_FRACTION( statesAlongSingleLeg.begin( )->first, legStartTime, 1.0E-14 );
            BOOST_CHECK_CLOSE_FRACTION( statesAlongSingleLeg.rbegin( )->first, legEndTime, 1.0E-14 );

            // Check if Keplerian state (slow elements) is the same for each output point
            Eigen::Vector6d previousKeplerianState = Eigen::Vector6d::Constant( TUDAT_NAN );
            for( auto it : statesAlongSingleLeg )
            {
                Eigen::Vector6d currentCartesianState = it.second;
                Eigen::Vector6d currentKeplerianState = tudat::orbital_element_conversions::convertCartesianToKeplerianElements(
                            currentCartesianState, sunGravitationalParameter );
                if( previousKeplerianState == previousKeplerianState )
                {
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                ( previousKeplerianState.segment( 0, 5 ) ),
                                ( currentKeplerianState.segment( 0, 5 ) ),
                                1.0E-10 );

                }
                previousKeplerianState = currentKeplerianState;
            }

            // Check if output meets boundary conditions
            Eigen::Vector3d depatureBodyPosition = bodies.at( bodyOrder.at( i ) )->getStateInBaseFrameFromEphemeris(
                        legStartTime ).segment( 0, 3 );
            Eigen::Vector3d arrivalBodyPosition = bodies.at( bodyOrder.at( i + 1 ) )->getStateInBaseFrameFromEphemeris(
                        legEndTime ).segment( 0, 3 );

            for( int j = 0; j < 3; j++ )
            {
                //TODO: Find out why tolerance needs to be so big for one of the legs
                BOOST_CHECK_SMALL( std::fabs( statesAlongSingleLeg.begin( )->second( j ) - depatureBodyPosition( j ) ), 20.0E3 );
                BOOST_CHECK_SMALL( std::fabs( statesAlongSingleLeg.rbegin( )->second( j ) - arrivalBodyPosition( j ) ), 20.0E3 );
            }
        }
    }
}


//! Test delta-V computation for MGA-1DSM Velocity Formulation trajectory model.
BOOST_AUTO_TEST_CASE( testMGA1DSMVFTrajectory1 )
{
    // Expected test result based on the ideal Messenger trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 8630.83256199051;

    // Create environment
    simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Set transfer order
    std::vector< std::string > bodyOrder = {
        "Earth", "Earth", "Venus", "Venus", "Mercury" };
    int numberOfNodes = bodyOrder.size( );

    // Create leg settings
    std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
    transferLegSettings.resize( numberOfNodes - 1 );
    transferLegSettings[ 0 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 1 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 2 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 3 ] = dsmVelocityBasedLeg( );

    // Create node settings
    std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;
    transferNodeSettings.resize( numberOfNodes );
    transferNodeSettings[ 0 ] = escapeAndDepartureNode( std::numeric_limits< double >::infinity( ), 0.0 );
    transferNodeSettings[ 1 ] = swingbyNode( );
    transferNodeSettings[ 2 ] = swingbyNode( );
    transferNodeSettings[ 3 ] = swingbyNode( );
    transferNodeSettings[ 4 ] = captureAndInsertionNode( std::numeric_limits< double >::infinity( ), 0.0 );

    std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun" );

    // Add the time of flight and start epoch, which are in JD.
    double JD = physical_constants::JULIAN_DAY;
    std::vector< double > nodeTimes;
    nodeTimes.push_back( ( 1171.64503236 - 0.5 ) * JD );
    nodeTimes.push_back( nodeTimes.at( 0 ) + 399.999999715 * JD );
    nodeTimes.push_back( nodeTimes.at( 1 ) + 178.372255301 * JD );
    nodeTimes.push_back( nodeTimes.at( 2 ) + 299.223139512  * JD );
    nodeTimes.push_back( nodeTimes.at( 3 ) + 180.510754824 * JD );

    std::vector< Eigen::VectorXd > transferLegFreeParameters;
    transferLegFreeParameters.resize( numberOfNodes - 1 );
    transferLegFreeParameters[ 0 ] = ( Eigen::VectorXd( 1 ) << 0.234594654679 ).finished( );
    transferLegFreeParameters[ 1 ] = ( Eigen::VectorXd( 1 ) << 0.0964769387134 ).finished( );
    transferLegFreeParameters[ 2 ] = ( Eigen::VectorXd( 1 ) << 0.829948744508 ).finished( );
    transferLegFreeParameters[ 3 ] = ( Eigen::VectorXd( 1 ) << 0.317174785637 ).finished( );

    std::vector< Eigen::VectorXd > transferNodeFreeParameters;
    transferNodeFreeParameters.resize( numberOfNodes );
    transferNodeFreeParameters[ 0 ] = ( Eigen::VectorXd( 3 ) <<1408.99421278, 0.37992647165 * 2 * 3.14159265358979, std::acos(  2 * 0.498004040298 - 1. ) - 3.14159265358979 / 2 ).finished( );
    transferNodeFreeParameters[ 1 ] = ( Eigen::VectorXd( 3 ) <<1.80629232251 * 6.378e6, 1.35077257078, 0.0 ).finished( );
    transferNodeFreeParameters[ 2 ] = ( Eigen::VectorXd( 3 ) <<3.04129845698 * 6.052e6, 1.09554368115, 0.0 ).finished( );
    transferNodeFreeParameters[ 3 ] = ( Eigen::VectorXd( 3 ) <<1.10000000891 * 6.052e6, 1.34317576594, 0.0 ).finished( );


    transferTrajectory->evaluateTrajectory(
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    printTransferParameterDefinition( transferLegSettings, transferNodeSettings );

    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, transferTrajectory->getTotalDeltaV( ), 1.0E-3 );
}

//! Test delta-V computation for another MGA-1DSM Velocity Formulation trajectory model.
BOOST_AUTO_TEST_CASE( testMGA1DSMVFTrajectory2 )
{
    // Expected test result based on the ideal Cassini 2 trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 8385.15784516116;

    // Create environment
    simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Set transfer order
    std::vector< std::string > bodyOrder = { "Earth", "Venus", "Venus",  "Earth", "Jupiter", "Saturn" };
    int numberOfNodes = bodyOrder.size( );

    // Create leg settings
    std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
    transferLegSettings.resize( numberOfNodes - 1 );
    transferLegSettings[ 0 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 1 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 2 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 3 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 4 ] = dsmVelocityBasedLeg( );


    std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;
    transferNodeSettings.resize( numberOfNodes );
    transferNodeSettings[ 0 ] = escapeAndDepartureNode( std::numeric_limits< double >::infinity( ), 0.0 );
    transferNodeSettings[ 1 ] = swingbyNode( );
    transferNodeSettings[ 2 ] = swingbyNode( );
    transferNodeSettings[ 3 ] = swingbyNode( );
    transferNodeSettings[ 4 ] = swingbyNode( );
    transferNodeSettings[ 5 ] = captureAndInsertionNode( std::numeric_limits< double >::infinity( ), 0.0 );


    std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun" );

    // Add the time of flight and start epoch, which are in JD.
    double JD = physical_constants::JULIAN_DAY;
    std::vector< double > nodeTimes;
    nodeTimes.push_back( ( -779.046753814506 - 0.5 ) * JD );
    nodeTimes.push_back( nodeTimes.at( 0 ) + 167.378952534645 * JD );
    nodeTimes.push_back( nodeTimes.at( 1 ) + 424.028254165204 * JD );
    nodeTimes.push_back( nodeTimes.at( 2 ) + 53.2897409769205  * JD );
    nodeTimes.push_back( nodeTimes.at( 3 ) + 589.766954923325 * JD );
    nodeTimes.push_back( nodeTimes.at( 4 ) + 2200.00000000000 * JD );

    std::vector< Eigen::VectorXd > transferLegFreeParameters;
    transferLegFreeParameters.resize( numberOfNodes - 1 );
    transferLegFreeParameters[ 0 ] = ( Eigen::VectorXd( 1 ) << 0.769483451363201 ).finished( );
    transferLegFreeParameters[ 1 ] = ( Eigen::VectorXd( 1 ) << 0.513289529822621 ).finished( );
    transferLegFreeParameters[ 2 ] = ( Eigen::VectorXd( 1 ) << 0.0274175362264024 ).finished( );
    transferLegFreeParameters[ 3 ] = ( Eigen::VectorXd( 1 ) << 0.263985256705873 ).finished( );
    transferLegFreeParameters[ 4 ] = ( Eigen::VectorXd( 1 ) << 0.599984695281461 ).finished( );

    std::vector< Eigen::VectorXd > transferNodeFreeParameters;
    transferNodeFreeParameters.resize( numberOfNodes );
    transferNodeFreeParameters[ 0 ] = ( Eigen::VectorXd( 3 ) <<3259.11446832345, 0.525976214695235 * 2 * 3.14159265358979, std::acos(  2 * 0.38086496458657 - 1 ) - 3.14159265358979 / 2 ).finished( );; // 1st leg.
    transferNodeFreeParameters[ 1 ] = ( Eigen::VectorXd( 3 ) <<1.34877968657176 * 6.052e6, -1.5937371121191, 0.0  ).finished( ); // 2nd leg.
    transferNodeFreeParameters[ 2 ] = ( Eigen::VectorXd( 3 ) <<1.05 * 6.052e6, -1.95952512232447, 0.0  ).finished( );// 3rd leg.
    transferNodeFreeParameters[ 3 ] = ( Eigen::VectorXd( 3 ) <<1.30730278372017 * 6.378e6, -1.55498859283059, 0.0  ).finished( ); //4th leg.
    transferNodeFreeParameters[ 4 ] = ( Eigen::VectorXd( 3 ) <<69.8090142993495 * 7.1492e7, -1.5134625299674, 0.0  ).finished( );//4th leg.
    transferNodeFreeParameters[ 5 ] = Eigen::VectorXd( 0 );


    transferTrajectory->evaluateTrajectory(
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    printTransferParameterDefinition( transferLegSettings, transferNodeSettings );


    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, transferTrajectory->getTotalDeltaV( ), 1.0E-3);
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
