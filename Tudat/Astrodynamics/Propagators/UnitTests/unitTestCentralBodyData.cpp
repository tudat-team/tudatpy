/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::propagators;

BOOST_AUTO_TEST_SUITE( test_central_body_data )

//! Test update order of central body data
//! NOTE: The transformatons done by CentralBodyData are effectively handled by
//! testCowellPopagatorCentralBodies
BOOST_AUTO_TEST_CASE( testCentralBodyData )
{
    // Define list of propagated  bodies and origins (with multi-levelled hierarchy).
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Moon" );
    bodiesToIntegrate.push_back( "Satellite1" );
    bodiesToIntegrate.push_back( "Phobos" );
    bodiesToIntegrate.push_back( "Earth" );
    bodiesToIntegrate.push_back( "Sun" );
    bodiesToIntegrate.push_back( "Satellite2" );
    bodiesToIntegrate.push_back( "Mars" );
    bodiesToIntegrate.push_back( "Deimos" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Earth" );
    centralBodies.push_back( "Moon" );
    centralBodies.push_back( "Mars" );
    centralBodies.push_back( "Sun" );
    centralBodies.push_back( "SSB" );
    centralBodies.push_back( "Deimos" );
    centralBodies.push_back( "Sun" );
    centralBodies.push_back( "Mars" );

    // Create dummy state functions.
    std::map< std::string,
              boost::function< Eigen::Matrix< double, 6, 1 >( const double ) > > stateFunctions;
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        stateFunctions[ bodiesToIntegrate[ i ] ]
            =  boost::lambda::constant( Eigen::Matrix< double, 6, 1 >::Zero( ) );
    }

    // Create central bodies object.
    boost::shared_ptr< CentralBodyData< double > > centralBodyData
        = boost::make_shared< CentralBodyData< double > >(
                centralBodies, bodiesToIntegrate, stateFunctions, boost::lambda::constant( Eigen::Vector6d::Zero( ) ), "SSB" );

    // Get update order.
    std::vector< int > updateOrder = centralBodyData->getUpdateOrder( );
    std::map< std::string, int > updateOrderMap;
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        updateOrderMap[ bodiesToIntegrate[ updateOrder[ i ] ] ] = i;
    }

    // Check whether the update order is consistent with input data.
    BOOST_CHECK_EQUAL( updateOrderMap[ "Moon" ] < updateOrderMap[ "Earth" ], 0 );
    BOOST_CHECK_EQUAL( updateOrderMap[ "Satellite1" ] < updateOrderMap[ "Moon" ], 0  );
    BOOST_CHECK_EQUAL( updateOrderMap[ "Phobos" ] < updateOrderMap[ "Mars" ], 0 );
    BOOST_CHECK_EQUAL( updateOrderMap[ "Deimos" ] < updateOrderMap[ "Mars" ], 0 );
    BOOST_CHECK_EQUAL( updateOrderMap[ "Satellite2" ] < updateOrderMap[ "Deimos" ], 0 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
