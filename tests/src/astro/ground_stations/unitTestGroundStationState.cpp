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

#include "tudat/basics/testMacros.h"

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/createLightTimeCalculator.h"
#include "tudat/simulation/simulation.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::coordinate_conversions;
using namespace tudat::basic_astrodynamics;
using namespace tudat::ground_stations;
using namespace tudat::simulation_setup;
using namespace tudat::unit_conversions;
using namespace tudat::spice_interface;


BOOST_AUTO_TEST_SUITE( test_ground_station_state )

//! Test if ground stations are correctly created.
BOOST_AUTO_TEST_CASE( test_GroundStationState )
{
    // Create Earth object
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    SystemOfBodies bodies;
    bodies.addBody( earth, "Earth" );

    // Central body characteristics (WGS84 Earth ellipsoid).
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    std::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSpheroidModel =
            std::make_shared< basic_astrodynamics::OblateSpheroidBodyShapeModel >(
                equatorialRadius, flattening );

    earth->setShapeModel( oblateSpheroidModel );


    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testGeodeticPosition( -63.667,
                                                convertDegreesToRadians( -7.26654999 ),
                                                convertDegreesToRadians( 72.36312094 ) );

    // Manually compute associated spherical position
    Eigen::Vector3d testSphericalPosition = coordinate_conversions::convertCartesianToSpherical(
                testCartesianPosition );
    testSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - testSphericalPosition( 1 );

    // Creatre ground stations: same position, but different representation
    createGroundStation( earth, "Station1", testCartesianPosition, cartesian_position );
    createGroundStation( earth, "Station2", testSphericalPosition, spherical_position );
    createGroundStation( earth, "Station3", testGeodeticPosition, geodetic_position );

    std::shared_ptr< GroundStationState > station1 = earth->getGroundStation( "Station1" )->getNominalStationState( );
    std::shared_ptr< GroundStationState > station2 = earth->getGroundStation( "Station2" )->getNominalStationState( );
    std::shared_ptr< GroundStationState > station3 = earth->getGroundStation( "Station3" )->getNominalStationState( );

    // Check if ground station representations are correctly converted
    Eigen::Vector3d position1, position2;

    {
        position1 = station1->getNominalCartesianPosition( );
        position2 = station2->getNominalCartesianPosition( );

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( position1( i ) - position2( i ), 1.0E-3 );
        }

        position1 = station1->getNominalCartesianPosition( );
        position2 = station3->getNominalCartesianPosition( );

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( position1( i ) - position2( i ), 1.0E-3 );
        }

        position1 = station2->getNominalCartesianPosition( );
        position2 = station3->getNominalCartesianPosition( );

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( position1( i ) - position2( i ), 1.0E-3 );
        }
    }

    {
        position1 = station1->getNominalGeodeticPosition( );
        position2 = station2->getNominalGeodeticPosition( );

        BOOST_CHECK_SMALL( position1( 0 ) - position2( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( position1( 1 ) - position2( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( position1( 2 ) - position2( 2 ), 1.0E-3 / equatorialRadius );

        position1 = station1->getNominalGeodeticPosition( );
        position2 = station3->getNominalGeodeticPosition( );

        BOOST_CHECK_SMALL( position1( 0 ) - position2( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( position1( 1 ) - position2( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( position1( 2 ) - position2( 2 ), 1.0E-3 / equatorialRadius );

        position1 = station2->getNominalGeodeticPosition( );
        position2 = station3->getNominalGeodeticPosition( );

        BOOST_CHECK_SMALL( position1( 0 ) - position2( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( position1( 1 ) - position2( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( position1( 2 ) - position2( 2 ), 1.0E-3 / equatorialRadius );
    }

    {
        position1 = station1->getNominalSphericalPosition( );
        position2 = station2->getNominalSphericalPosition( );

        BOOST_CHECK_SMALL( position1( 0 ) - position2( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( position1( 1 ) - position2( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( position1( 2 ) - position2( 2 ), 1.0E-3 / equatorialRadius );

        position1 = station1->getNominalSphericalPosition( );
        position2 = station3->getNominalSphericalPosition( );

        BOOST_CHECK_SMALL( position1( 0 ) - position2( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( position1( 1 ) - position2( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( position1( 2 ) - position2( 2 ), 1.0E-3 / equatorialRadius );

        position1 = station2->getNominalSphericalPosition( );
        position2 = station3->getNominalSphericalPosition( );

        BOOST_CHECK_SMALL( position1( 0 ) - position2( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( position1( 1 ) - position2( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( position1( 2 ) - position2( 2 ), 1.0E-3 / equatorialRadius );
    }
}

//! Test if global state function for ground station is correctly created.
BOOST_AUTO_TEST_CASE( test_GroundStationGlobalState )
{
    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Create Earth object
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    SystemOfBodies bodies;
    bodies.addBody( earth, "Earth" );

    // Central body characteristics (WGS84 Earth ellipsoid).
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    std::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSpheroidModel =
            std::make_shared< basic_astrodynamics::OblateSpheroidBodyShapeModel >(
                equatorialRadius, flattening );

    earth->setShapeModel( oblateSpheroidModel );
    earth->setEphemeris( std::make_shared< ephemerides::SpiceEphemeris >(
                             "Earth", "SSB", false, true, true, "ECLIPJ2000" ) );
    earth->setRotationalEphemeris( std::make_shared< ephemerides::SpiceRotationalEphemeris >(
                                       "ECLIPJ2000", "IAU_Earth" ) );


    // Define ground station state
    const Eigen::Vector3d groundStationPosition( 1917032.190, 6029782.349, -801376.113 );
    Eigen::Vector6d groundStationState;
    groundStationState << groundStationPosition, 0.0, 0.0, 0.0;

    // Create ground station
    createGroundStation( earth, "Station1", groundStationPosition, cartesian_position );

    // Make state function of ground station w.r.t. SSB in inertial frame
    std::function< Eigen::Matrix< double, 6, 1 >( const double ) > stateFunction =
            simulation_setup::getLinkEndCompleteEphemerisFunction(
                observation_models::LinkEndId( "Earth", "Station1" ), bodies );

    // Compare state function with manual computation.
    Eigen::Vector6d currentGlobalState, currentGlobalStateFromFunction;
    for( double testTime = 1.0E7; testTime < 5.0E7; testTime += 2.5E6 )
    {
        currentGlobalState = earth->getEphemeris( )->getCartesianState( testTime ) +
                ephemerides::transformStateToInertialOrientation( groundStationState, testTime, earth->getRotationalEphemeris( ) );
        currentGlobalStateFromFunction = stateFunction( testTime );
        for( unsigned int i = 0; i < 6; i++ )
        {
            BOOST_CHECK_EQUAL( currentGlobalState( i ), currentGlobalStateFromFunction( i ) );
        }
    }
}


Eigen::Vector3d getCustomDisplacementFunction( const double time )
{
    return ( Eigen::Vector3d( )<<4.3, 1.3, -1.3 ).finished( ) *
            std::sin( time / ( 10.0 * physical_constants::JULIAN_DAY ) );
}

Eigen::Vector3d getBodyFixedSpecialRelativisticDisplacementFunction(
    const double time, const simulation_setup::SystemOfBodies& bodies, const Eigen::Vector3d& nominalStationState )
{
    Eigen::Quaterniond bodyFixedToInertial = bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToBaseFrame( time );
    Eigen::Vector3d inertialStationPosition = bodyFixedToInertial * nominalStationState;
    Eigen::Vector3d earthBarycentricVelocity = bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( time ).segment( 3, 3 );

    return bodyFixedToInertial.inverse( ) * ( inertialStationPosition.dot( earthBarycentricVelocity ) * earthBarycentricVelocity * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT );
}

Eigen::Vector3d getExpectedPiecewiseConstantDisplacement( const int j )
{
    Eigen::Vector3d expectedStateDifference;
    switch( j )
    {
    case 0:
        expectedStateDifference = Eigen::Vector3d::Zero( );
        break;
    case 1:
        expectedStateDifference = ( Eigen::Vector3d( ) << 5.3, -1.2, 7.6 ).finished( );
        break;
    case 2:
        expectedStateDifference = ( Eigen::Vector3d( ) << 4.5, 3.2, -4.6 ).finished( );
        break;
    case 3:
        expectedStateDifference = ( Eigen::Vector3d( ) << 0.2, 6.4, 12.4 ).finished( );
        break;
    }
    return expectedStateDifference;
}

//! Test if global state function for ground station is correctly created.
BOOST_AUTO_TEST_CASE( test_GroundStationTimeVaryingState )
{
    //Load spice kernels.
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

    //Define environment settings
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    for( unsigned int bodyDeformationType = 0; bodyDeformationType < 2; bodyDeformationType++ )
    {
        for( unsigned int i = 0; i < 5; i++ )
        {
            BodyListSettings bodySettings =
                    getDefaultBodySettings( bodyNames );
            if( bodyDeformationType == 1 )
            {
                std::map< int, std::pair< double, double > > loveNumbers;
                loveNumbers[ 2 ] = std::make_pair( 0.7, 0.0 );
                bodySettings.at( "Earth" )->bodyDeformationSettings.push_back(
                            std::make_shared< BasicSolidBodyDeformationSettings >(
                                std::vector< std::string >( { "Moon" } ), loveNumbers ) );
                bodySettings.at( "Earth" )->bodyDeformationSettings.push_back(
                            std::make_shared< BasicSolidBodyDeformationSettings >(
                                std::vector< std::string >( { "Sun" } ), loveNumbers ) );
            }
            SystemOfBodies bodies = createSystemOfBodies( bodySettings );
            const Eigen::Vector3d nominalStationState( 1917032.190, 6029782.349, -801376.113 );

            std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings;

            //! Test case 0
            double referenceTime = 10.0 * physical_constants::JULIAN_YEAR;
            Eigen::Vector3d linearMotion = ( Eigen::Vector3d( ) << 2.1, 0.5 , -1.5 ).finished( ) / physical_constants::JULIAN_YEAR;

            //! Test case 1
            std::map< double, Eigen::Vector3d > positionDiscontinuities;
            positionDiscontinuities[ -2.4 * physical_constants::JULIAN_YEAR ] = ( Eigen::Vector3d( ) << 4.3, -9.2, 4.1 ).finished( );
            positionDiscontinuities[ -20.0 * physical_constants::JULIAN_DAY ] = ( Eigen::Vector3d( ) << 5.3, -1.2, 7.6 ).finished( );
            positionDiscontinuities[ -10.0 * physical_constants::JULIAN_DAY ] = ( Eigen::Vector3d( ) << 4.5, 3.2, -4.6 ).finished( );
            positionDiscontinuities[ 10.0 * physical_constants::JULIAN_YEAR ] = ( Eigen::Vector3d( ) << 0.2, 6.4, 12.4 ).finished( );

            if( i == 0 || i == 3 )
            {
                stationMotionSettings.push_back(
                            std::make_shared< LinearGroundStationMotionSettings >( linearMotion, referenceTime ) );
            }
            if( i == 1 || i == 3 )
            {
                stationMotionSettings.push_back(
                            std::make_shared< PiecewiseConstantGroundStationMotionSettings >( positionDiscontinuities ) );
            }
            if( i == 2 || i == 3 )
            {
                stationMotionSettings.push_back(
                            std::make_shared< CustomGroundStationMotionSettings >( &getCustomDisplacementFunction ) );
            }
            if( i == 4 )
            {
                stationMotionSettings.push_back(
                    std::make_shared< BodyCentricToBarycentricGroundStationMotionSettings >( "", false ) );
            }
            stationMotionSettings.push_back( std::make_shared< BodyDeformationStationMotionSettings >( true ) );

            std::shared_ptr< GroundStationSettings > stationSettings = std::make_shared< GroundStationSettings >(
                        "Station1", nominalStationState, coordinate_conversions::cartesian_position, stationMotionSettings );
            createGroundStation( bodies.at( "Earth" ), stationSettings );

            std::shared_ptr< tudat::ground_stations::GroundStationState > stationState = bodies.at( "Earth" )->getGroundStation(
                        "Station1" )->getNominalStationState( );

            std::vector< double > testTimes ={ -3.0 * physical_constants::JULIAN_YEAR, -20.0 * physical_constants::JULIAN_DAY, 0.0,
                                               20.0 * physical_constants::JULIAN_YEAR };

            std::vector< std::shared_ptr< tudat::basic_astrodynamics::BodyDeformationModel > > deformationModels;
            if( bodyDeformationType == 1 )
            {
                deformationModels = bodies.at( "Earth" )->getBodyDeformationModels( );
            }
            if( i == 0 )
            {
                for( unsigned int j = 0; j < testTimes.size( ); j++ )
                {
                    Eigen::Vector3d currentState = stationState->getCartesianPositionInTime( testTimes.at( j ), bodies.getFrameOrigin( ) );
                    Eigen::Vector3d stateDifference = currentState - nominalStationState;
                    Eigen::Vector3d expectedStateDifference = ( linearMotion * ( testTimes.at( j ) - referenceTime ) );
                    if( bodyDeformationType == 1 )
                    {
                        expectedStateDifference += deformationModels.at( 0 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                        expectedStateDifference += deformationModels.at( 1 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                    }
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateDifference, expectedStateDifference, 1.0E-7 );

                }
            }
            else if( i == 1 )
            {
                for( unsigned int j = 0; j < testTimes.size( ); j++ )
                {
                    Eigen::Vector3d currentState = stationState->getCartesianPositionInTime( testTimes.at( j ), bodies.getFrameOrigin( ) );
                    Eigen::Vector3d stateDifference = currentState - nominalStationState;
                    Eigen::Vector3d expectedStateDifference = getExpectedPiecewiseConstantDisplacement( j );
                    if( bodyDeformationType == 1 )
                    {
                        expectedStateDifference += deformationModels.at( 0 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                        expectedStateDifference += deformationModels.at( 1 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                    }
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateDifference, expectedStateDifference, 1.0E-7 );

                }
            }
            else if( i == 2 )
            {
                for( unsigned int j = 0; j < testTimes.size( ); j++ )
                {
                    Eigen::Vector3d currentState = stationState->getCartesianPositionInTime( testTimes.at( j ), bodies.getFrameOrigin( ) );
                    Eigen::Vector3d stateDifference = currentState - nominalStationState;
                    Eigen::Vector3d expectedStateDifference = getCustomDisplacementFunction( testTimes.at( j ) );
                    if( bodyDeformationType == 1 )
                    {
                        expectedStateDifference += deformationModels.at( 0 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                        expectedStateDifference += deformationModels.at( 1 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                    }
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateDifference, expectedStateDifference, 1.0E-7 );
                }
            }
            else if( i == 3 )
            {
                for( unsigned int j = 0; j < testTimes.size( ); j++ )
                {
                    Eigen::Vector3d currentState = stationState->getCartesianPositionInTime( testTimes.at( j ), bodies.getFrameOrigin( ) );
                    Eigen::Vector3d stateDifference = currentState - nominalStationState;
                    Eigen::Vector3d expectedStateDifference =
                            getCustomDisplacementFunction( testTimes.at( j ) ) + getExpectedPiecewiseConstantDisplacement( j ) +
                            ( linearMotion * ( testTimes.at( j ) - referenceTime ) );
                    if( bodyDeformationType == 1 )
                    {
                        expectedStateDifference += deformationModels.at( 0 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                        expectedStateDifference += deformationModels.at( 1 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                    }
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateDifference, expectedStateDifference, 1.0E-6 );
                }
            }
            else if( i == 4 )
            {
                for( unsigned int j = 0; j < testTimes.size( ); j++ )
                {
                    {
                        Eigen::Vector3d currentState = stationState->getCartesianPositionInTime( testTimes.at( j ), bodies.getFrameOrigin( ) );
                        Eigen::Vector3d stateDifference = currentState - nominalStationState;
                        Eigen::Vector3d expectedStateDifference = getBodyFixedSpecialRelativisticDisplacementFunction( testTimes.at( j ), bodies, nominalStationState );
                        if( bodyDeformationType == 1 )
                        {
                            expectedStateDifference += deformationModels.at( 0 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                            expectedStateDifference += deformationModels.at( 1 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                        }

                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateDifference, expectedStateDifference, 1.0E-6 );
                    }

                    {
                        Eigen::Vector3d currentState = stationState->getCartesianPositionInTime( testTimes.at( j ), "Earth" );
                        Eigen::Vector3d stateDifference = currentState - nominalStationState;
                        Eigen::Vector3d expectedStateDifference = Eigen::Vector3d::Zero( );
                        if( bodyDeformationType == 1 )
                        {
                            expectedStateDifference += deformationModels.at( 0 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                            expectedStateDifference += deformationModels.at( 1 )->calculateDisplacement( testTimes.at( j ), nominalStationState );
                        }

                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateDifference, expectedStateDifference, 1.0E-6 );
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
