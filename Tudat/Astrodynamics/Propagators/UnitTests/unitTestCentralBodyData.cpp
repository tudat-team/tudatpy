#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include "Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Astrodynamics/Bodies/body.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::propagators;
using namespace tudat::bodies;

BOOST_AUTO_TEST_SUITE( test_hardware_creation )

BOOST_AUTO_TEST_CASE( testHardwareCreation )
{

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

    std::map< std::string, boost::function< Eigen::Matrix< double, 6, 1 >( const double ) > > stateFunctions;

    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        stateFunctions[ bodiesToIntegrate[ i ] ] =  boost::lambda::constant( Eigen::Matrix< double, 6, 1 >::Zero( ) );
    }

    boost::shared_ptr< CentralBodyData< double > > centralBodyData = boost::make_shared< CentralBodyData< double > >(
                centralBodies, bodiesToIntegrate, stateFunctions );

    std::vector< int > updateOrder = centralBodyData->getUpdateOrder( );

    std::map< std::string, int > updateOrderMap;
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        updateOrderMap[ bodiesToIntegrate[ updateOrder[ i ] ] ] = i;
    }

    BOOST_CHECK_EQUAL( updateOrderMap[ "Moon" ] < updateOrderMap[ "Earth" ], 0 );
    BOOST_CHECK_EQUAL( updateOrderMap[ "Satellite1" ] < updateOrderMap[ "Moon" ], 0  );
    BOOST_CHECK_EQUAL( updateOrderMap[ "Phobos" ] < updateOrderMap[ "Mars" ], 0 );
    BOOST_CHECK_EQUAL( updateOrderMap[ "Deimos" ] < updateOrderMap[ "Mars" ], 0 );
    BOOST_CHECK_EQUAL( updateOrderMap[ "Satellite2" ] < updateOrderMap[ "Deimos" ], 0 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
