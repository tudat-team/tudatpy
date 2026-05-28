/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
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
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
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

using namespace tudat;
using namespace tudat::unit_tests;
using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_differenced_frequency_of_arrival_partials )

//! Test partial derivatives of differenced frequency of arrival observable directly against sub-model partials.
BOOST_AUTO_TEST_CASE( testDifferencedFrequencyOfArrivalDirectPartial )
{
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
            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    defaultBodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Set frequency of transmitter (Mars) to X-Band: 8.4 GHz
    double transmitterFrequency = 8.4e9;
    bodies.getBody( "Mars" )->getVehicleSystems( )->setTransmittedFrequencyCalculator(
            std::make_shared< ground_stations::ConstantFrequencyInterpolator >( transmitterFrequency ) );

    // Define link ends for the two one-way measured frequency models
    LinkEnds linkEnds1;
    linkEnds1[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    linkEnds1[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    LinkEnds linkEnds2;
    linkEnds2[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-43" );
    linkEnds2[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    // Define link ends for differenced model
    LinkEnds differencedLinkEnds;
    differencedLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    differencedLinkEnds[ receiver2 ] = std::make_pair< std::string, std::string >( "Earth", "DSS-43" );
    differencedLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back(
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    // Create observation settings for the two sub-models and the differenced model
    std::shared_ptr< ObservationModelSettings > freqSettings1 =
            oneWayDopplerMeasuredFrequencySettings( linkEnds1, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > freqSettings2 =
            oneWayDopplerMeasuredFrequencySettings( linkEnds2, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > differencedSettings =
            differencedFrequencyOfArrivalObservationSettings( differencedLinkEnds, lightTimeCorrectionSettings );

    // Create observation models
    std::shared_ptr< ObservationModel< 1, double, Time > > firstFreqModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( freqSettings1, bodies );
    std::shared_ptr< ObservationModel< 1, double, Time > > secondFreqModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( freqSettings2, bodies );
    std::shared_ptr< ObservationModel< 1, double, Time > > differencedModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( differencedSettings, bodies );

    // Create parameters to estimate
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Mars",
            spice_interface::getBodyCartesianStateAtEpoch( "Mars", "SSB", bodies.getFrameOrientation( ), "None", initialEphemerisTime ),
            "SSB",
            bodies.getFrameOrientation( ) ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );

    // Create partials for all three models
    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
               std::shared_ptr< PositionPartialScaling > >
            firstFreqPartials = ObservationPartialCreator< 1, double, Time >::createObservationPartials(
                    firstFreqModel, bodies, parametersToEstimate );
    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
               std::shared_ptr< PositionPartialScaling > >
            secondFreqPartials = ObservationPartialCreator< 1, double, Time >::createObservationPartials(
                    secondFreqModel, bodies, parametersToEstimate );
    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
               std::shared_ptr< PositionPartialScaling > >
            differencedPartials = ObservationPartialCreator< 1, double, Time >::createObservationPartials(
                    differencedModel, bodies, parametersToEstimate );

    std::cout << "Partials sizes: " << differencedPartials.first.size( ) << " "
              << firstFreqPartials.first.size( ) << " " << secondFreqPartials.first.size( ) << std::endl;

    // Compute observations and partials
    double receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
    std::vector< double > linkEndTimesDifferenced, linkEndTimesFreq1, linkEndTimesFreq2;
    std::vector< Eigen::Vector6d > linkEndStatesDifferenced, linkEndStatesFreq1, linkEndStatesFreq2;

    Eigen::VectorXd differencedObs = differencedModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, linkEndTimesDifferenced, linkEndStatesDifferenced );
    differencedPartials.second->update( linkEndStatesDifferenced, linkEndTimesDifferenced, receiver, differencedObs );
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > differencedPartialValues =
            differencedPartials.first.begin( )->second->calculatePartial(
                    linkEndStatesDifferenced, linkEndTimesDifferenced, receiver );

    Eigen::VectorXd firstFreqObs = firstFreqModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, linkEndTimesFreq1, linkEndStatesFreq1 );
    firstFreqPartials.second->update( linkEndStatesFreq1, linkEndTimesFreq1, receiver, firstFreqObs );
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > firstFreqPartialValues =
            firstFreqPartials.first.begin( )->second->calculatePartial( linkEndStatesFreq1, linkEndTimesFreq1, receiver );

    Eigen::VectorXd secondFreqObs = secondFreqModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, linkEndTimesFreq2, linkEndStatesFreq2 );
    secondFreqPartials.second->update( linkEndStatesFreq2, linkEndTimesFreq2, receiver, secondFreqObs );
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > secondFreqPartialValues =
            secondFreqPartials.first.begin( )->second->calculatePartial( linkEndStatesFreq2, linkEndTimesFreq2, receiver );

    // Verify transmitter times match (both sub-models evaluated at same receiver time, so transmitter times should be close)
    BOOST_CHECK_CLOSE_FRACTION( firstFreqPartialValues.at( 0 ).second, differencedPartialValues.at( 0 ).second,
                                ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
    BOOST_CHECK_CLOSE_FRACTION( secondFreqPartialValues.at( 0 ).second, differencedPartialValues.at( 1 ).second,
                                ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );

    // Verify: differenced partial = first partial - second partial
    // The first entry in differencedPartialValues should equal firstFreqPartialValues
    // The second entry (negated) should equal secondFreqPartialValues
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            firstFreqPartialValues.at( 0 ).first,
            ( differencedPartialValues.at( 0 ).first ),
            ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            secondFreqPartialValues.at( 0 ).first,
            ( -differencedPartialValues.at( 1 ).first ),
            ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
}

//! Test partial derivatives of differenced frequency of arrival observable, using general test suite.
BOOST_AUTO_TEST_CASE( testDifferencedFrequencyOfArrivalPartials )
{
    Eigen::VectorXd parameterPerturbationMultipliers = ( Eigen::VectorXd( 4 ) << 100.0, 100.0, 1.0, 100.0 ).finished( );

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 3 );
    groundStations[ 0 ] = std::make_pair( "Mars", "MSL" );
    groundStations[ 1 ] = std::make_pair( "Earth", "DSS-13" );
    groundStations[ 2 ] = std::make_pair( "Earth", "DSS-35" );

    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true, 1.0, false, true );

        // Set transmitter frequency on Mars ground station
        bodies.getBody( "Mars" )->getGroundStation( "MSL" )->setTransmittingFrequencyCalculator(
                std::make_shared< ground_stations::ConstantFrequencyInterpolator >( 8.4e9 ) );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 0 ];
        linkEnds[ receiver ] = groundStations[ 1 ];
        linkEnds[ receiver2 ] = groundStations[ 2 ];

        // Generate differenced frequency of arrival model
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::shared_ptr< ObservationModel< 1 > > differencedFreqModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        differencedFrequencyOfArrivalObservationSettings(
                                linkEnds,
                                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >(
                                        { std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                                basic_astrodynamics::tdb_scale ),
                        bodies );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7 );

        testObservationPartials< 1 >( differencedFreqModel,
                                      bodies,
                                      fullEstimatableParameterSet,
                                      linkEnds,
                                      differenced_frequency_of_arrival,
                                      1.5E-3,
                                      true,
                                      true,
                                      10.0,
                                      parameterPerturbationMultipliers );
    }

    // Test partials with real ephemerides (without test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false, 1.0, false, true );

        // Set transmitter frequency on Mars ground station
        bodies.getBody( "Mars" )->getGroundStation( "MSL" )->setTransmittingFrequencyCalculator(
                std::make_shared< ground_stations::ConstantFrequencyInterpolator >( 8.4e9 ) );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 0 ];
        linkEnds[ receiver ] = groundStations[ 1 ];
        linkEnds[ receiver2 ] = groundStations[ 2 ];

        // Generate differenced frequency of arrival model
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::shared_ptr< ObservationModel< 1 > > differencedFreqModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        differencedFrequencyOfArrivalObservationSettings(
                                linkEnds,
                                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >(
                                        { std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                                basic_astrodynamics::tdb_scale ),
                        bodies );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7 );

        testObservationPartials< 1 >( differencedFreqModel,
                                      bodies,
                                      fullEstimatableParameterSet,
                                      linkEnds,
                                      differenced_frequency_of_arrival,
                                      1.5E-3,
                                      false,
                                      true,
                                      10.0,
                                      parameterPerturbationMultipliers );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
