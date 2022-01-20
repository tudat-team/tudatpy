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

#include <limits>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"
#include "tudat/simulation/simulation.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace simulation_setup;
using namespace input_output;

BOOST_AUTO_TEST_SUITE( test_synchronous_rotational_ephemeris )

BOOST_AUTO_TEST_CASE( test_SynchronousRotationModel )
{

    // Load Spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    // Create relevant bodies
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Jupiter" );

    std::vector< std::string > bodiesToTest;
    bodiesToTest.push_back( "Io" );
    bodiesToTest.push_back( "Europa" );
    bodiesToTest.push_back( "Ganymede" );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames );
    for( unsigned int i = 0; i < bodiesToTest.size( ); i++ )
    {
        bodySettings.at( bodiesToTest.at( i ) )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                    "Jupiter", "ECLIPJ2000", "IAU_" + bodiesToTest.at( i ) );
        bodySettings.at( bodiesToTest.at( i ) )->ephemerisSettings->resetFrameOrigin( "Jupiter" );
    }
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    // Test rotation model when body is in propagation, and outside of propagation
    for( unsigned int areBodiesInPropagation = 0; areBodiesInPropagation < 2; areBodiesInPropagation++ )
    {
        setAreBodiesInPropagation( bodies, areBodiesInPropagation );

        // Test different bodies
        for( unsigned int i = 0; i < bodiesToTest.size( ); i++ )
        {
            // Get synchronous rotation model and ephemeris
            std::shared_ptr< SynchronousRotationalEphemeris > currentRotationalEphemeris =
                    std::dynamic_pointer_cast< SynchronousRotationalEphemeris >(
                        bodies.at( bodiesToTest.at( i ) )->getRotationalEphemeris( ) );
            std::shared_ptr< Ephemeris > currentEphemeris = bodies.at( bodiesToTest.at( i ) )->getEphemeris( );

            // Test for various times (and anomalies)
            for( int j = 0; j < 100; j ++ )
            {
                // Determine current test time
                double testTime = 86400.0 * static_cast< double >( j );

                // Compute current satellite state
                Eigen::Vector6d currentSatelliteState;
                if( !areBodiesInPropagation )
                {
                    currentSatelliteState = currentEphemeris->getCartesianState( testTime );
                }
                else
                {
                    bodies.at( bodiesToTest.at( i ) )->setStateFromEphemeris( testTime );
                    bodies.at( "Jupiter" )->setStateFromEphemeris( testTime );
                    currentSatelliteState = bodies.at( bodiesToTest.at( i ) )->getState( ) -
                             bodies.at( "Jupiter" )->getState(  );
                }

                // Retrieve rotation matrix
                Eigen::Matrix3d currentRotationToBodyFixedFrame =
                        currentRotationalEphemeris->getRotationToTargetFrame( testTime ).toRotationMatrix( );

                // Transform r, v and r x v vectors to body-fixed frame
                Eigen::Vector3d bodyFixedRadialVector =
                        currentRotationToBodyFixedFrame * currentSatelliteState.segment( 0, 3 ).normalized( );
                Eigen::Vector3d bodyFixedVelocityVector =
                        currentRotationToBodyFixedFrame * currentSatelliteState.segment( 3, 3 ).normalized( );
                Eigen::Vector3d bodyFixedOrbitalAngularMomentumVelocityVector =
                        currentRotationToBodyFixedFrame *
                        ( Eigen::Vector3d( currentSatelliteState.segment( 0, 3 ) ).cross(
                              Eigen::Vector3d( currentSatelliteState.segment( 3, 3 ) ) ) ).normalized( );

                // Test vector magnitudes
                BOOST_CHECK_SMALL( std::fabs( bodyFixedRadialVector( 0 ) + 1.0 ), 10.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( bodyFixedRadialVector( 1 ) ), 10.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( bodyFixedRadialVector( 2 ) ), 10.0 * std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_SMALL( std::fabs( bodyFixedVelocityVector( 0 ) ), 0.01 );
                BOOST_CHECK_SMALL( std::fabs( bodyFixedVelocityVector( 1 ) + 1.0 ), 0.01 );
                BOOST_CHECK_SMALL( std::fabs( bodyFixedVelocityVector( 2 ) ), 10.0 * std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_SMALL( std::fabs( bodyFixedOrbitalAngularMomentumVelocityVector( 0 ) ), 10.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( bodyFixedOrbitalAngularMomentumVelocityVector( 1 ) ), 10.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( bodyFixedOrbitalAngularMomentumVelocityVector( 2 ) - 1.0 ), 10.0 * std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
