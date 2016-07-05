#define BOOST_TEST_MAIN

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::spice_interface;

BOOST_AUTO_TEST_SUITE( test_frame_manager )

BOOST_AUTO_TEST_CASE( test_FrameManager )
{

    std::string kernelsPath = input_output::getSpiceKernelPath( );

    //Load spice kernels.
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    std::map< std::string, boost::shared_ptr< Ephemeris > > ephemerisList;

    basic_mathematics::Vector6d barycentricSunState = getBodyCartesianStateAtEpoch( "Sun", getBaseFrameName( ), "ECLIPJ2000", "NONE", 0.0 );
    ephemerisList[ "Sun" ] = boost::make_shared< ConstantEphemeris >( barycentricSunState, getBaseFrameName( ), "ECLIPJ2000" );

    basic_mathematics::Vector6d sunCentricEarthState = getBodyCartesianStateAtEpoch( "Earth", "Sun", "ECLIPJ2000", "NONE", 0.0 );
    ephemerisList[ "Earth" ] = boost::make_shared< ConstantEphemeris >( sunCentricEarthState, "Sun", "ECLIPJ2000" );

    basic_mathematics::Vector6d earthCentricMoonState = getBodyCartesianStateAtEpoch( "Moon", "Earth", "ECLIPJ2000", "NONE", 0.0 );
    ephemerisList[ "Moon" ] = boost::make_shared< ConstantEphemeris >( earthCentricMoonState, "Earth", "ECLIPJ2000" );

    basic_mathematics::Vector6d earthCentricLageosState = basic_mathematics::Vector6d::Zero( );
    earthCentricLageosState( 1 ) = 2.5E6;
    earthCentricLageosState( 2 ) = 4.0E6;
    ephemerisList[ "LAGEOS" ] = boost::make_shared< ConstantEphemeris >( earthCentricLageosState, "Earth", "ECLIPJ2000" );

    basic_mathematics::Vector6d moonCentricLroState = basic_mathematics::Vector6d::Zero( );
    moonCentricLroState( 0 ) = 1.0E6;
    moonCentricLroState( 1 ) = 2.0E6;

    ephemerisList[ "LRO" ] = boost::make_shared< ConstantEphemeris >( moonCentricLroState, "Moon", "ECLIPJ2000" );

    basic_mathematics::Vector6d sunCentricMarsState = getBodyCartesianStateAtEpoch( "Mars", "Sun", "ECLIPJ2000", "NONE", 0.0 );
    ephemerisList[ "Mars" ] = boost::make_shared< ConstantEphemeris >( sunCentricMarsState, "Sun", "ECLIPJ2000" );

    basic_mathematics::Vector6d marsCentricPhobosState = basic_mathematics::Vector6d::Zero( );
    marsCentricPhobosState( 0 ) = 2.3E5;
    marsCentricPhobosState( 1 ) = 2.9E4;
    marsCentricPhobosState( 2 ) = 600;

    ephemerisList[ "Phobos" ] = boost::make_shared< ConstantEphemeris >( marsCentricPhobosState, "Mars", "ECLIPJ2000" );

    boost::shared_ptr< ReferenceFrameManager > frameManager = boost::make_shared< ReferenceFrameManager >( ephemerisList );

    std::map< std::string, int > expectedFrameLevel;
    expectedFrameLevel[ "Sun" ] = 0;
    expectedFrameLevel[ "Earth" ] = 1;
    expectedFrameLevel[ "Mars" ] = 1;
    expectedFrameLevel[ "LAGEOS" ] = 2;
    expectedFrameLevel[ "Moon" ] = 2;
    expectedFrameLevel[ "Phobos" ] = 2;
    expectedFrameLevel[ "LRO" ] = 3;

    for( std::map< std::string, int >::iterator it = expectedFrameLevel.begin( ); it != expectedFrameLevel.end( ); it++ )
    {
        if( frameManager->getFrameLevel( it->first ).first != it->second )
        {
            throw std::runtime_error(
                        "Error when identifying frame level of " + it->first + " found " +
                        boost::lexical_cast< std::string >( frameManager->getFrameLevel( it->first ).first ) + " expected" +
                        boost::lexical_cast< std::string >( it->second ) );
        }
    }

    std::vector< std::string > frames;
    frames.push_back( "Earth" );
    frames.push_back( "Mars" );
    std::pair< std::string, int > commonFrame = frameManager->getNearestCommonFrame( frames );
    BOOST_CHECK_EQUAL( commonFrame.first, "Sun" );
    BOOST_CHECK_EQUAL( commonFrame.second, 0 );

    frames.clear( );
    frames.push_back( "Moon" );
    frames.push_back( "Sun" );

    commonFrame = frameManager->getNearestCommonFrame( frames );

    BOOST_CHECK_EQUAL( commonFrame.first, "Sun" );
    BOOST_CHECK_EQUAL( commonFrame.second, 0 );

    frames.clear( );
    frames.push_back( "Earth" );
    frames.push_back( "LRO" );

    commonFrame = frameManager->getNearestCommonFrame( frames );

    BOOST_CHECK_EQUAL( commonFrame.first, "Earth" );
    BOOST_CHECK_EQUAL( commonFrame.second, 1 );

    basic_mathematics::Vector6d testState = frameManager->getEphemeris< >( "Moon", "LAGEOS" )->getCartesianStateFromEphemeris( 0.0 );
    basic_mathematics::Vector6d expectedState = earthCentricLageosState - earthCentricMoonState;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testState, expectedState, std::numeric_limits< double >::epsilon( ) );

    testState = frameManager->getEphemeris( "Phobos", "Sun" )->getCartesianStateFromEphemeris( 0.0 );
    expectedState = sunCentricMarsState + marsCentricPhobosState;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testState, ( -1.0 * expectedState ), std::numeric_limits< double >::epsilon( ) );

    testState = frameManager->getEphemeris( "Sun", "Phobos" )->getCartesianStateFromEphemeris( 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testState, expectedState, std::numeric_limits< double >::epsilon( ) );

    testState = frameManager->getEphemeris( "Sun", "Earth" )->getCartesianStateFromEphemeris( 0.0 );
    expectedState = sunCentricEarthState;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testState, expectedState, std::numeric_limits< double >::epsilon( ) );

    testState = frameManager->getEphemeris( getBaseFrameName( ), "Earth" )->getCartesianStateFromEphemeris( 0.0 );
    expectedState = barycentricSunState + sunCentricEarthState;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testState, expectedState, std::numeric_limits< double >::epsilon( ) );

    testState = frameManager->getEphemeris( "Earth", getBaseFrameName( ) )->getCartesianStateFromEphemeris( 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testState, ( -1.0 * expectedState ), std::numeric_limits< double >::epsilon( ) );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
