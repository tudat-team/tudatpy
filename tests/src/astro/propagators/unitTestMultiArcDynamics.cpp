/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"

namespace tudat
{

namespace unit_tests
{

//! Using declarations.
using namespace interpolators;
using namespace numerical_integrators;
using namespace spice_interface;
using namespace ephemerides;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace propagators;

BOOST_AUTO_TEST_SUITE( test_multi_arc_dynamics )

BOOST_AUTO_TEST_CASE( testKeplerMultiArcDynamics )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    for( unsigned testCase = 0; testCase < 3; testCase++ )
    {
        std::vector< std::string > bodyNames;
        bodyNames.push_back( "Earth" );
        bodyNames.push_back( "Moon" );

        // Specify initial time
        double initialEphemerisTime = 1.0E7;
        double finalEphemerisTime = 2.0E7;
        double maximumTimeStep = 3600.0;
        double buffer = 5.0 * maximumTimeStep;

        // Create bodies needed in simulation
        BodyListSettings bodySettings =
                getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
        std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( bodySettings.at( "Moon" )->ephemerisSettings )->
                resetFrameOrigin( "Earth" );
        bodySettings.at( "Moon" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
        bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                    Eigen::Vector6d::Zero( ) );

        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        

        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
        accelerationsOfMoon[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationMap[ "Moon" ] = accelerationsOfMoon;

        std::vector< std::string > bodiesToIntegrate, centralBodies;
        bodiesToIntegrate.push_back( "Moon" );
        centralBodies.push_back( "SSB" );

        std::vector< double > integrationArcStarts, integrationArcEnds;

        double integrationStartTime = initialEphemerisTime + 1.0E4;
        double integrationEndTime = finalEphemerisTime - 1.0E4;
        double arcDuration = 1.0E6;
        double arcOverlap = 1.0E4;

        double currentStartTime = integrationStartTime;
        double currentEndTime = integrationStartTime + arcDuration;

        do
        {
            integrationArcStarts.push_back( currentStartTime );
            integrationArcEnds.push_back( currentEndTime );

            currentStartTime = currentEndTime - arcOverlap;
            currentEndTime = currentStartTime + arcDuration;
        }
        while( currentEndTime < integrationEndTime );

        unsigned int numberOfIntegrationArcs = integrationArcStarts.size( );

        std::vector< Eigen::VectorXd > systemInitialStates;
        std::vector< Eigen::Vector6d > initialKeplerElements;

        systemInitialStates.resize( numberOfIntegrationArcs );
        initialKeplerElements.resize( numberOfIntegrationArcs );

        double earthGravitationalParameter =  bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        for(  unsigned int j = 0; j < numberOfIntegrationArcs; j++ )
        {
            systemInitialStates[ j ]  = spice_interface::getBodyCartesianStateAtEpoch(
                        bodiesToIntegrate[ 0 ], "Earth", "ECLIPJ2000", "NONE", integrationArcStarts[ j ] );
            initialKeplerElements[ j ] = (
                        orbital_element_conversions::convertCartesianToKeplerianElements(
                            Eigen::Vector6d( systemInitialStates[ j ] ),
                            earthGravitationalParameter ) );
        }

        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToIntegrate, centralBodies );

        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
        for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
        {
            arcPropagationSettingsList.push_back(
                        std::make_shared< TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                          systemInitialStates.at( i ), integrationArcEnds.at( i ) ) );
        }

        // For case 0: test multi-arc estimation with same integration settings for each arc
        if( testCase == 0 )
        {
            std::shared_ptr< IntegratorSettings< > > integratorSettings =
                    std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, initialEphemerisTime, 120.0 );
            MultiArcDynamicsSimulator< > dynamicsSimulator(
                        bodies, integratorSettings, std::make_shared< MultiArcPropagatorSettings< double > >(
                            arcPropagationSettingsList ), integrationArcStarts );
        }
        // For case 1: test multi-arc estimation with different integration settings object for each arc
        else if( testCase == 1 )
        {
            std::vector< std::shared_ptr< IntegratorSettings< > > > integratorSettingsList;
            for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
            {
                integratorSettingsList.push_back( std::make_shared< IntegratorSettings< > >
                                                  ( rungeKutta4, integrationArcStarts.at( i ), 120.0 ) );
            }
            MultiArcDynamicsSimulator< > dynamicsSimulator(
                        bodies, integratorSettingsList, std::make_shared< MultiArcPropagatorSettings< double > >(
                            arcPropagationSettingsList ) );
        }
        // For case 0: test multi-arc estimation with same integration settings for each arc, and arc initial state interpolated
        // from previous state
        else  if( testCase == 2 )
        {
            std::shared_ptr< IntegratorSettings< > > integratorSettings =
                    std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, initialEphemerisTime, 120.0 );
            MultiArcDynamicsSimulator< > dynamicsSimulator(
                        bodies, integratorSettings, std::make_shared< MultiArcPropagatorSettings< double > >(
                            arcPropagationSettingsList, true ), integrationArcStarts );
        }


        std::shared_ptr< Ephemeris > moonEphemeris = bodies.at( "Moon" )->getEphemeris( );

        double testStartTime, testEndTime;
        double testTimeStep = 10000.0;
        double timeBuffer = 1000.0;

        Eigen::Vector6d stateDifference;

        for( unsigned int i = 0; i < numberOfIntegrationArcs ; i++ )
        {
            if( i == 0 )
            {
                testStartTime = integrationArcStarts.at( i ) + timeBuffer;
            }
            else
            {
                testStartTime = integrationArcEnds.at( i - 1 ) + timeBuffer;
            }

            if( i == numberOfIntegrationArcs - 1 )
            {
                testEndTime = integrationArcEnds.at( i ) - timeBuffer;
            }
            else
            {
                testEndTime = integrationArcStarts.at( i + 1 ) - timeBuffer;
            }

            // Check if output corresponds to expected analytical solution
            if( testCase < 2 || i == 0 )
            {
                double currentTestTime = testStartTime;
                while( currentTestTime < testEndTime )
                {
                    stateDifference = ( moonEphemeris->getCartesianState( currentTestTime ) ) -
                            ( orbital_element_conversions::convertKeplerianToCartesianElements(
                                  propagateKeplerOrbit( initialKeplerElements.at( i ), currentTestTime - integrationArcStarts.at( i ),
                                                        earthGravitationalParameter ), earthGravitationalParameter ) );
                    for( int i = 0; i < 3; i++ )
                    {
                        BOOST_CHECK_SMALL( stateDifference( i ), 1.0E-4 );
                        BOOST_CHECK_SMALL( stateDifference( i + 3 ), 1.0E-10 );

                    }
                    currentTestTime += testTimeStep;
                }
            }
            // Check if arc information is properly passed to next arc
            else
            {
                double currentTestTime = testStartTime;
                while( currentTestTime < testEndTime )
                {
                    stateDifference = ( moonEphemeris->getCartesianState( currentTestTime ) ) -
                            ( orbital_element_conversions::convertKeplerianToCartesianElements(
                                  propagateKeplerOrbit( initialKeplerElements.at( 0 ), currentTestTime - integrationArcStarts.at( 0 ),
                                                        earthGravitationalParameter ), earthGravitationalParameter ) );
                    for( int i = 0; i < 3; i++ )
                    {
                        BOOST_CHECK_SMALL( stateDifference( i ), 1.0E-3);
                        BOOST_CHECK_SMALL( stateDifference( i + 3 ), 1.0E-9 );

                    }
                    currentTestTime += testTimeStep;
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
