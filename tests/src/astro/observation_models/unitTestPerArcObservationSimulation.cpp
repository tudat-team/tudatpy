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
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
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
using namespace tudat::coordinate_conversions;
using namespace tudat::statistics;

BOOST_AUTO_TEST_SUITE( test_observation_noise_models )


std::vector< std::vector< double > > splitArcTimes(
        const std::vector< double >& observationTimes )
{

    std::vector< std::vector< double > > perArcObservationTimes;
    std::vector< double > currentArc;

    currentArc.push_back( observationTimes.at( 0 ) );
    for( unsigned int i = 1; i < observationTimes.size( ); i++ )
    {
        if( observationTimes.at( i ) - observationTimes.at( i - 1 ) > 60.0 && currentArc.size( ) > 0 )
        {
            perArcObservationTimes.push_back( currentArc );
            currentArc.clear( );
        }
        currentArc.push_back( observationTimes.at( i ) );

        if( i == observationTimes.size( ) - 1 )
        {
            perArcObservationTimes.push_back( currentArc );
        }
    }
    return perArcObservationTimes;
}

std::vector< double > getArcLengths(
        const std::vector< std::vector< double > >& perArcTimes )
{
    std::vector< double > arcLengths;
    for( unsigned int i = 0; i < perArcTimes.size( ); i++ )
    {
        arcLengths.push_back( perArcTimes.at( i ).at( perArcTimes.at( i ).size( ) - 1 ) -
                              perArcTimes.at( i ).at( 0 ) );
    }
    return arcLengths;
}

//! Test whether observation noise is correctly added when simulating noisy observations
BOOST_AUTO_TEST_CASE( testObservationNoiseModels )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialEphemerisTime = double( 1.0E7 );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Earth",
                spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                initialEphemerisTime, 2.0 * mathematical_constants::PI /
                ( physical_constants::JULIAN_DAY ) );

    Eigen::Vector6d spacecraftOrbitalElements;
    spacecraftOrbitalElements( semiMajorAxisIndex ) = 2000.0E3;
    spacecraftOrbitalElements( eccentricityIndex ) = 0.05;
    spacecraftOrbitalElements( inclinationIndex ) = 1.5;
    spacecraftOrbitalElements( argumentOfPeriapsisIndex ) = 0.0;
    spacecraftOrbitalElements( longitudeOfAscendingNodeIndex ) = 0.0;
    spacecraftOrbitalElements( trueAnomalyIndex ) = 0.0;
    bodySettings.addSettings( "LunarOrbiter" );
    bodySettings.at( "LunarOrbiter" )->ephemerisSettings =
            keplerEphemerisSettings( spacecraftOrbitalElements, 0.0, getBodyGravitationalParameter( "Moon" ), "Moon" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    // Creatre ground stations: same position, but different representation
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );

    // Define link ends.
    LinkEnds testLinkEnds;
    testLinkEnds[ transmitter ] = std::pair< std::string, std::string >( std::make_pair( "Earth", "Station1" ) );
    testLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "LunarOrbiter", "" );

    // Create observation settings
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                           one_way_range,  testLinkEnds ) );

    // Create observation simulators
    std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > >  observationSimulators =
            createObservationSimulators( observationSettingsList, bodies );


    // Define observation simulation settings (observation type, link end, times and reference link end)
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > idealMeasurementSimulationInput;
    idealMeasurementSimulationInput.push_back(
                std::make_shared< PerArcObservationSimulationSettings< double > >(
                    one_way_range, testLinkEnds,
                    physical_constants::JULIAN_YEAR,
                    physical_constants::JULIAN_YEAR + 28.0 * physical_constants::JULIAN_DAY, 60.0,
                    elevationAngleViabilitySettings( std::make_pair( "Earth", "Station1" ), 0.0 ) ) );


    std::shared_ptr< ObservationCollection< > > idealObservationsAndTimes = simulateObservations< double, double >(
                idealMeasurementSimulationInput, observationSimulators, bodies );
    std::vector< double > idealObservationTimes = idealObservationsAndTimes->getConcatenatedTimeVector( );
    std::vector< std::vector< double > > perArcIdealObservationTimes = splitArcTimes( idealObservationTimes );
    std::vector< double > idealArcLengths = getArcLengths( perArcIdealObservationTimes );


    std::shared_ptr< ObservationCollection< > > caseTwoObservationsAndTimes;

    for( int test = 0; test < 4; test++ )
    {
        double minimumArcDuration = TUDAT_NAN;
        double maximumArcDuration = TUDAT_NAN;
        std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > additionalViabilitySettingsList;

        if( test == 0 )
        {
            minimumArcDuration = 12.0 * 3600.0;
        }

        if( test == 1 )
        {
            maximumArcDuration = 12.0 * 3600.0;
        }

        if( test == 2 || test == 3 )
        {
            minimumArcDuration = 11.5 * 3600.0;
            maximumArcDuration = 12.5 * 3600.0;
        }

        if( test == 3 )
        {
            additionalViabilitySettingsList.push_back( bodyOccultationViabilitySettings(
                                                           std::make_pair( "LunarOrbiter", "" ), "Moon"  ) );

        }

        std::cout<<"******************** TEST ********************* "<<std::endl;
        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
        measurementSimulationInput.push_back(
                    std::make_shared< PerArcObservationSimulationSettings< double > >(
                        one_way_range, testLinkEnds,
                        physical_constants::JULIAN_YEAR,
                        physical_constants::JULIAN_YEAR + 28.0 * physical_constants::JULIAN_DAY, 60.0,
                        elevationAngleViabilitySettings( std::make_pair( "Earth", "Station1" ), 0.0 ),
                        minimumArcDuration, maximumArcDuration, TUDAT_NAN, observation_models::unidentified_link_end,
                        additionalViabilitySettingsList ) );
        std::shared_ptr< ObservationCollection< > > testObservationsAndTimes = simulateObservations< double, double >(
                    measurementSimulationInput, observationSimulators, bodies );
        if( test == 2 )
        {
            caseTwoObservationsAndTimes = testObservationsAndTimes;
        }
        std::vector< double > testObservationTimes = testObservationsAndTimes->getConcatenatedTimeVector( );
        std::vector< std::vector< double > > perArcTestObservationTimes = splitArcTimes( testObservationTimes );
        std::vector< double > testArcLengths = getArcLengths( perArcTestObservationTimes );

        int testArcCounter = 0;

        if( test == 1 )
        {
            BOOST_CHECK_EQUAL( perArcIdealObservationTimes.size( ),
                               perArcTestObservationTimes.size( ) );
        }

        for( unsigned int i = 0; i < perArcIdealObservationTimes.size( ); i++ )
        {
            if( test == 0 )
            {
                if( idealArcLengths.at( i ) >= 12.0 * 3600.0 )
                {
                    BOOST_CHECK_EQUAL( perArcIdealObservationTimes.at( i ).size( ),
                                       perArcTestObservationTimes.at( testArcCounter ).size( ) );
                    for( unsigned int j = 0; j < perArcIdealObservationTimes.at( i ).size( ); j++ )
                    {
                        BOOST_CHECK_EQUAL( perArcIdealObservationTimes.at( i ).at( j ),
                                           perArcTestObservationTimes.at( testArcCounter ).at( j ) );
                    }

                    testArcCounter++;
                }
            }
            if( test == 1 )
            {
                if( idealArcLengths.at( i ) < 12.0 * 3600.0 )
                {
                    BOOST_CHECK_EQUAL( perArcIdealObservationTimes.at( i ).size( ),
                                       perArcTestObservationTimes.at( i ).size( ) );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( testArcLengths.at( i ) - 12.0 * 3600.0 ), 60.0 );
                }

                for( unsigned int j = 0; j < perArcTestObservationTimes.at( i ).size( ); j++ )
                {
                    BOOST_CHECK_EQUAL( perArcIdealObservationTimes.at( i ).at( j ),
                                       perArcTestObservationTimes.at( i ).at( j ) );
                }
            }

            if( test == 2 )
            {
                if( idealArcLengths.at( i ) >= 11.5 * 3600.0 )
                {
                    if( idealArcLengths.at( i ) < 12.5 * 3600.0 )
                    {
                        BOOST_CHECK_EQUAL( perArcIdealObservationTimes.at( i ).size( ),
                                           perArcTestObservationTimes.at( testArcCounter ).size( ) );
                    }
                    else
                    {
                        BOOST_CHECK_SMALL( std::fabs( testArcLengths.at( testArcCounter ) - 12.5 * 3600.0 ), 60.0 );
                    }

                    for( unsigned int j = 0; j < perArcTestObservationTimes.at( testArcCounter ).size( ); j++ )
                    {
                        BOOST_CHECK_EQUAL( perArcIdealObservationTimes.at( i ).at( j ),
                                           perArcTestObservationTimes.at( testArcCounter ).at( j ) );
                    }
                    testArcCounter++;
                }
            }
        }

        if( test == 0 || test == 2 )
        {
            BOOST_CHECK_EQUAL( testArcCounter, perArcTestObservationTimes.size( ) );
        }

        if( test == 3 )
        {
            std::vector< double > testObservationTimes = testObservationsAndTimes->getConcatenatedTimeVector( );
            std::vector< double > referenceObservationTimes = caseTwoObservationsAndTimes->getConcatenatedTimeVector( );

            std::shared_ptr< observation_models::ObservationModel< 1 > > observationModel =
                    std::dynamic_pointer_cast< ObservationSimulator< 1 > >( observationSimulators.at( 0 ) )->getObservationModel(
                        testLinkEnds );

            unsigned int testIndex = 0;

            Eigen::VectorXd currentObservation;
            std::vector< Eigen::Vector6d > vectorOfStates;
            std::vector< double > vectorOfTimes;
            bool isObservationFeasible;

            std::vector< std::shared_ptr< observation_models::ObservationViabilityCalculator > > additionalViabilityCalculators =
                    observation_models::createObservationViabilityCalculators(
                        bodies, testLinkEnds, one_way_range, additionalViabilitySettingsList );

            BOOST_CHECK_EQUAL( ( testObservationTimes.size( ) < referenceObservationTimes.size( ) ), true );

            for( unsigned int i = 0; i < referenceObservationTimes.size( ); i++ )
            {
                currentObservation = observationModel->computeIdealObservationsWithLinkEndData(
                            referenceObservationTimes.at( i ), receiver, vectorOfTimes, vectorOfStates );
                isObservationFeasible = isObservationViable( vectorOfStates, vectorOfTimes, additionalViabilityCalculators );
                if( testObservationTimes.at( testIndex ) == referenceObservationTimes.at( i )  )
                {
                    BOOST_CHECK_EQUAL( isObservationFeasible, true );
                    testIndex++;
                }
                else
                {
                    BOOST_CHECK_EQUAL( isObservationFeasible, false );
                }

                if( testIndex == testObservationTimes.size( ) )
                {
                    break;
                }
            }

            BOOST_CHECK_EQUAL( testIndex, testObservationTimes.size( ) );

        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

}

}

