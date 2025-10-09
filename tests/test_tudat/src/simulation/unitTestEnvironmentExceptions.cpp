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

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <functional>
#include <cmath>

#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{


// Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_environment_exceptions )


//! Unit test to check if desaturation deltaV values are estimated correctly
BOOST_AUTO_TEST_CASE( test_EphemerisInterpolationException )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 1.0e7;
    double finalEphemerisTime = initialEphemerisTime + 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime, finalEphemerisTime, "SSB", "ECLIPJ2000", 30.0  );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    {
      bool caughtStateException = false;
      try
      {
          bodies.at("Earth")->getState( );
      }
      catch( std::runtime_error& exception )
      {
          caughtStateException = true;
          std::cout<<exception.what( )<<std::endl;
      }
      BOOST_CHECK_EQUAL( caughtStateException, true );
    }


    std::shared_ptr< TabulatedCartesianEphemeris< double, double > > moonEphemeris =
            std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >( bodies.at( "Moon" )->getEphemeris( ) );
    moonEphemeris->resetInterpolator(
            createStateInterpolatorFromSpice( "Moon", initialEphemerisTime, finalEphemerisTime, 30.0, moonEphemeris->getReferenceFrameOrigin( ),
                                       moonEphemeris->getReferenceFrameOrientation( ), std::make_shared< LagrangeInterpolatorSettings >(
                                               8, false, huntingAlgorithm, lagrange_cubic_spline_boundary_interpolation, throw_exception_at_boundary ) ) );

    for( unsigned int test = 0; test < 2; test++ )
    {
        // Create ground stations.
        std::vector< std::string > groundStationNames;
        groundStationNames.push_back( "Station1" );

        createGroundStation( bodies.at( "Earth" ),
                             "Station1",
                             ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ),
                             coordinate_conversions::geodetic_position );


        // Define link ends.
        std::vector< LinkDefinition > observationLinkEnds;

        LinkDefinition linkEnds;
        linkEnds[ receiver ] = std::pair< std::string, std::string >( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Moon", "" );
        observationLinkEnds.push_back( linkEnds );

        std::map< ObservableType, std::vector< LinkDefinition > > linkEndsPerObservable;
        linkEndsPerObservable[ one_way_range ].push_back( observationLinkEnds[ 0 ] );


        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
        for( std::map< ObservableType, std::vector< LinkDefinition > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
             linkEndIterator != linkEndsPerObservable.end( );
             linkEndIterator++ )
        {
            ObservableType currentObservable = linkEndIterator->first;

            std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
            for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
            {
                observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                        currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ) ) );
            }
        }


        // Compute list of observation times.
        std::vector< double > baseTimeList;

        double buffer = ( test == 0 ) ? -300.0 : 3600.0;
        double observationTime = initialEphemerisTime + buffer;
        double observationInterval = 60.0;
        while( observationTime < finalEphemerisTime - buffer )
        {
            baseTimeList.push_back( observationTime );
            observationTime += observationInterval;
        }

        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
        for( std::map< ObservableType, std::vector< LinkDefinition > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
             linkEndIterator != linkEndsPerObservable.end( );
             linkEndIterator++ )
        {
            ObservableType currentObservable = linkEndIterator->first;
            std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
            for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
            {
                measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings<> >(
                        currentObservable, currentLinkEndsList[ i ], baseTimeList, receiver ) );
            }
        }

        std::vector< std::shared_ptr< ObservationSimulatorBase< > > > observationSimulators =
                createObservationSimulators(
                observationSettingsList, bodies );

        // Simulate observations.
        bool caughtException = false;
        try
        {
            std::shared_ptr< ObservationCollection<> > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, observationSimulators, bodies );
        }
        catch( std::runtime_error& exception )
        {
            caughtException = true;
            std::cout<<exception.what( )<<std::endl;
        }
        if( test == 0 )
        {
            BOOST_CHECK_EQUAL( caughtException, true );
        }
        else if( test == 1 )
        {
            BOOST_CHECK_EQUAL( caughtException, false );
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
