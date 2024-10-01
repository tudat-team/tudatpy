/*    Copyright (c) 2010-2023, Delft University of Technology
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

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <boost/lambda/lambda.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"

#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

#include "tudat/support/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;
using namespace tudat::input_output;

BOOST_AUTO_TEST_SUITE( test_dsn_n_way_averaged_doppler_observation_partials)

std::vector< double > getRetransmissionDelays( const double evaluationTime, const int numberOfRetransmitters )
{
    std::vector< double > retransmissionDelays;

        for( int i = 0; i < numberOfRetransmitters; i++ )
        {
            retransmissionDelays.push_back( evaluationTime * 5.0E-17 * static_cast< double >( i + 1 ) );
        }
    return retransmissionDelays;
}

//! Test partial derivatives of DSN N-way Doppler observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testDsnNWayAveragedDopplerPartials )
{

    Eigen::VectorXd parameterPerturbationMultipliers =
            ( Eigen::VectorXd( 4 ) << 100.0, 100.0, 1.0, 100.0 ).finished( );

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "DSS-55" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

    double initialEphemerisTime = 544845633.0;
    double finalEphemerisTime = 544869060.0;
    double stateEvaluationTime = initialEphemerisTime + 8.0e3;

    // Read ODF file - used just for the automatic creation of ground station ramp frequency calculator
    std::shared_ptr< OdfRawFileContents > rawOdfFileContents =
            std::make_shared< OdfRawFileContents >( tudat::paths::getTudatTestDataPath( ) + "mromagr2017_097_1335xmmmv1.odf" );

    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, initialEphemerisTime,
                                                  finalEphemerisTime, stateEvaluationTime, true );

        // Process ODF file
        std::shared_ptr< ProcessedOdfFileContents< Time > > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents< Time > >( rawOdfFileContents, "MSL", true );
        // Create ground stations
        setTransmittingFrequenciesInGroundStations( processedOdfFileContents, bodies.getBody( "Earth" ) );
        // Set turnaround ratios in spacecraft (ground station)
        std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared< system_models::VehicleSystems >( );
        vehicleSystems->setTransponderTurnaroundRatio( &getDsnDefaultTurnaroundRatios );
        bodies.getBody( "Mars" )->getGroundStation( "MSL" )->setVehicleSystems( vehicleSystems );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 0 ];
        linkEnds[ retransmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate DSN n-way averaged doppler model
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList;
        lightTimeCorrectionsList.push_back(
                    std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) );

        std::shared_ptr< ObservationModel< 1, double, Time > > dsnNWayAveragedDopplerModel =
                observation_models::ObservationModelCreator< 1, double, Time >::createObservationModel(
                    std::make_shared< observation_models::DsnNWayAveragedDopplerObservationSettings >(
                        linkEnds,
                        lightTimeCorrectionsList ) , bodies );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, stateEvaluationTime );

        testObservationPartials< 1 >(
                dsnNWayAveragedDopplerModel, bodies, fullEstimatableParameterSet, linkEnds,
                dsn_n_way_averaged_doppler, 1.0E-4, true, true, 1000.0, parameterPerturbationMultipliers,
                getDsnNWayAveragedDopplerAncillarySettings(
                    std::vector< FrequencyBands >{ x_band, x_band }, x_band, 7.0e9, 60.0, getRetransmissionDelays( initialEphemerisTime, 1 ) ),
                    stateEvaluationTime );


    }

    // Test partials with real ephemerides (without test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, initialEphemerisTime,
                                                  finalEphemerisTime, stateEvaluationTime, true );

        // Process ODF file
        std::shared_ptr< ProcessedOdfFileContents< Time > > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents< Time > >( rawOdfFileContents, "MSL", true );
        // Create ground stations
        setTransmittingFrequenciesInGroundStations( processedOdfFileContents, bodies.getBody( "Earth" ) );
        // Set turnaround ratios in spacecraft (ground station)
        std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared< system_models::VehicleSystems >( );
        vehicleSystems->setTransponderTurnaroundRatio( &getDsnDefaultTurnaroundRatios );
        bodies.getBody( "Mars" )->getGroundStation( "MSL" )->setVehicleSystems( vehicleSystems );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 0 ];
        linkEnds[ retransmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way range model
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList;
        lightTimeCorrectionsList.push_back(
                    std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) );

       std::shared_ptr< ObservationModel< 1, double, Time > > dsnNWayAveragedDopplerModel =
                observation_models::ObservationModelCreator< 1, double, Time >::createObservationModel(
                    std::make_shared< observation_models::DsnNWayAveragedDopplerObservationSettings >(
                        linkEnds,
                        lightTimeCorrectionsList ) , bodies );

       // Create parameter objects.
       std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
               createEstimatableParameters( bodies, stateEvaluationTime );

       testObservationPartials< 1 >(
               dsnNWayAveragedDopplerModel, bodies, fullEstimatableParameterSet, linkEnds,
               dsn_n_way_averaged_doppler, 1.0E-4, false, true, 1000.0, parameterPerturbationMultipliers,
               getDsnNWayAveragedDopplerAncillarySettings(
                   std::vector< FrequencyBands >{ x_band, x_band }, x_band, 7.0e9, 60.0, getRetransmissionDelays( initialEphemerisTime, 1 ) ),
                   stateEvaluationTime );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
