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

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;


BOOST_AUTO_TEST_SUITE( test_one_way_doppler_model )

double integrationTimeFunction( const double currentObservationTime )
{
    return 60.0 + 30.0 * ( currentObservationTime - 3.0 * 86400.0 ) / ( 7.0 * 86400.0 );
}

BOOST_AUTO_TEST_CASE( testOneWayDoppplerModel )
{
    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = std::make_pair( "Earth" , ""  );
    linkEnds[ receiver ] = std::make_pair( "Mars" , ""  );

    // Create range rate observation settings and model
    std::shared_ptr< ObservationModelSettings > rangeRateObservableSettings = std::make_shared<
            OneWayDifferencedRangeRateObservationSettings >( linkEnds, std::bind( &integrationTimeFunction, std::placeholders::_1 ) );
    std::shared_ptr< ObservationModel< 1, double, double> > rangeRateObservationModel =
            ObservationModelCreator< 1, double, double>::createObservationModel(
                rangeRateObservableSettings, bodies );

    // Create range rate observation settings and model
    std::shared_ptr< ObservationModelSettings > rangeObservableSettings = std::make_shared< ObservationModelSettings >
            ( one_way_range, linkEnds );
    std::shared_ptr< ObservationModel< 1, double, double> > rangeObservationModel =
            ObservationModelCreator< 1, double, double>::createObservationModel(
                rangeObservableSettings, bodies );

    // Test observable for both fixed link ends
    for( unsigned testCase = 0; testCase < 2; testCase++ )
    {
        for( double observationTime = 86400.0; observationTime <= 86400.0; observationTime += observationTime )
        {

            std::cout << "TEST: ************************************* " << testCase << " " << observationTime << std::endl;
            double dopplerCountInterval = integrationTimeFunction( observationTime );
            double arcStartObservationTime = observationTime - dopplerCountInterval;
            double arcEndObservationTime = observationTime;
            std::vector< double > rangeRateLinkEndTimes;
            std::vector< Eigen::Vector6d > rangeRateLinkEndStates;

            std::vector< double > rangeStartLinkEndTimes;
            std::vector< Eigen::Vector6d > rangeStartLinkEndStates;

            std::vector< double > rangeEndLinkEndTimes;
            std::vector< Eigen::Vector6d > rangeEndLinkEndStates;

            // Define link end
            LinkEndType referenceLinkEnd;
            if( testCase == 0 )
            {
                referenceLinkEnd = transmitter;
            }
            else
            {
                referenceLinkEnd = receiver;
            }

            // Compute observable
            double rangeRateObservable = rangeRateObservationModel->computeObservationsWithLinkEndData(
                        observationTime, referenceLinkEnd, rangeRateLinkEndTimes, rangeRateLinkEndStates )( 0 );

            double arcEndRange = rangeObservationModel->computeObservationsWithLinkEndData(
                        arcEndObservationTime, referenceLinkEnd, rangeEndLinkEndTimes, rangeEndLinkEndStates )( 0 );
            double arcStartRange = rangeObservationModel->computeObservationsWithLinkEndData(
                        arcStartObservationTime, referenceLinkEnd, rangeStartLinkEndTimes, rangeStartLinkEndStates )( 0 );

            for( unsigned linkCompare = 0; linkCompare < 2; linkCompare++ )
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeStartLinkEndStates.at( linkCompare ), rangeRateLinkEndStates.at( linkCompare ), 1.0E-15 );
                BOOST_CHECK_CLOSE_FRACTION( rangeStartLinkEndTimes.at( linkCompare ), rangeRateLinkEndTimes.at( linkCompare ), 1.0E-15 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeEndLinkEndStates.at( linkCompare ), rangeRateLinkEndStates.at( linkCompare + 2 ), 1.0E-15 );
                BOOST_CHECK_CLOSE_FRACTION( rangeEndLinkEndTimes.at( linkCompare ), rangeRateLinkEndTimes.at( linkCompare + 2 ), 1.0E-15 );
            }

            double manualDifferencedRange = ( arcEndRange - arcStartRange ) / dopplerCountInterval;

            // Test numerical derivative against Doppler observable
            BOOST_CHECK_SMALL( std::fabs( manualDifferencedRange - rangeRateObservable ), 1.0E-4 * dopplerCountInterval );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}


