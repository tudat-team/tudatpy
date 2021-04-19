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

namespace tudat
{
namespace unit_tests
{

// Additional namespaces to be used.
using namespace tudat::input_output;
using namespace tudat::simulation_setup;
using namespace mission_segments;

simulation_setup::SystemOfBodies getApproximatePlanetBodyMap( )
{
    using namespace ephemerides;
    using namespace gravitation;


    SystemOfBodies bodies;
    bodies.createEmptyBody( "Sun" );
    bodies.createEmptyBody( "Mercury" );
    bodies.createEmptyBody( "Venus" );
    bodies.createEmptyBody( "Earth" );
    bodies.createEmptyBody( "Jupiter" );
    bodies.createEmptyBody( "Saturn" );

    bodies.getBody( "Sun" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ) ) );
    bodies.getBody( "Mercury" )->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                            ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury ) );
    bodies.getBody( "Venus" )->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                          ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus ) );
    bodies.getBody( "Earth" )->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                          ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter ) );
    bodies.getBody( "Jupiter" )->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                            ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter ) );
    bodies.getBody( "Saturn" )->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                           ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn ) );

    bodies.getBody( "Sun" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( 1.32712428e20 ) );
    bodies.getBody( "Mercury" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( 2.2321e13 ) );
    bodies.getBody( "Venus" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( 3.24860e14 ) );
    bodies.getBody( "Earth" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( 3.9860119e14 ) );
    bodies.getBody( "Jupiter" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( 1.267e17 ) );
    bodies.getBody( "Saturn" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( 3.79e16 ) );

    return bodies;
}


//! Test implementation of trajectory class
BOOST_AUTO_TEST_SUITE( test_trajectory )

BOOST_AUTO_TEST_CASE( testMGATrajectory_New )
{

    // Set tolerance. Due to ephemeris differences, higher accuracy cannot be achieved.
    const double tolerance = 15.0e-2;

    // Expected test result based on the ideal Cassini 1 trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 4930.72686847243;

    // Create environment
    simulation_setup::SystemOfBodies bodies = getApproximatePlanetBodyMap( );

    // Set transfer order
    std::vector< std::string > bodyOrder = {
        "Earth", "Venus", "Venus", "Earth", "Jupiter", "Saturn" };
    int numberOfNodes = bodyOrder.size( );

    // Create leg settings (all unpowered)
    std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
    transferLegSettings.resize( numberOfNodes - 1 );
    transferLegSettings[ 0 ] = unpoweredLeg( );
    transferLegSettings[ 1 ] = unpoweredLeg( );
    transferLegSettings[ 2 ] = unpoweredLeg( );
    transferLegSettings[ 3 ] = unpoweredLeg( );
    transferLegSettings[ 4 ] = unpoweredLeg( );

    // Define minimum periapsis altitudes for flybys;
    std::map< std::string, double >minimumPeriapses;
    minimumPeriapses[ "Venus" ] = 6351800.0;
    minimumPeriapses[ "Earth" ] = 6778000.0;
    minimumPeriapses[ "Jupiter" ] =  600000000.0;

    std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;
    transferNodeSettings.resize( numberOfNodes );
    transferNodeSettings[ 0 ] = escapeAndDepartureNode( std::numeric_limits< double >::infinity( ), 0.0 );
    transferNodeSettings[ 1 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 1 ) ) );
    transferNodeSettings[ 2 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 2 ) ) );
    transferNodeSettings[ 3 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 3 ) ) );
    transferNodeSettings[ 4 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 4 ) ) );
    transferNodeSettings[ 5 ] = captureAndInsertionNode( 1.0895e8 / 0.02, 0.98 );

    // Create list of node times
    double JD = physical_constants::JULIAN_DAY;
    std::vector< double > nodeTimes;
    nodeTimes.push_back( -789.8117 * JD );
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

    std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun",
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    transferTrajectory->evaluateTrajectory(
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    printTransferParameterDefinition( transferLegSettings, transferNodeSettings );

    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, transferTrajectory->getTotalDeltaV( ), tolerance );
}


//! Test delta-V computation for MGA-1DSM Velocity Formulation trajectory model.
BOOST_AUTO_TEST_CASE( testMGA1DSMVFTrajectory1 )
{
    // Set tolerance. Due to ephemeris differences, higher accuracy cannot be achieved.
    const double tolerance = 2.0e-3;

    // Expected test result based on the ideal Messenger trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 8630.83256199051;

    // Create environment
    simulation_setup::SystemOfBodies bodies = getApproximatePlanetBodyMap( );

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

    // Add the time of flight and start epoch, which are in JD.
    double JD = physical_constants::JULIAN_DAY;
    std::vector< double > nodeTimes;
    nodeTimes.push_back( 1171.64503236 * JD );
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

    std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun",
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    transferTrajectory->evaluateTrajectory(
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    printTransferParameterDefinition( transferLegSettings, transferNodeSettings );

    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, transferTrajectory->getTotalDeltaV( ), tolerance );

}

//! Test delta-V computation for another MGA-1DSM Velocity Formulation trajectory model.
BOOST_AUTO_TEST_CASE( testMGA1DSMVFTrajectory2 )
{
    // Set tolerance. Due to ephemeris differences, higher accuracy cannot be achieved.
    const double tolerance = 4.0e-2;

    // Expected test result based on the ideal Cassini 2 trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 8385.15784516116;

    // Create environment
    simulation_setup::SystemOfBodies bodies = getApproximatePlanetBodyMap( );

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


    // Add the time of flight and start epoch, which are in JD.
    double JD = physical_constants::JULIAN_DAY;
    std::vector< double > nodeTimes;
    nodeTimes.push_back( -779.046753814506 * JD );
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

    std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun",
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    transferTrajectory->evaluateTrajectory(
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    printTransferParameterDefinition( transferLegSettings, transferNodeSettings );


    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, transferTrajectory->getTotalDeltaV( ), tolerance );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
