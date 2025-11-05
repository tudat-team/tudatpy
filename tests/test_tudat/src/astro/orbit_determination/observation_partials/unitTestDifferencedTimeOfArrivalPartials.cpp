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

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <boost/lambda/lambda.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

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

using namespace tudat;
using namespace tudat::unit_tests;
using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_differenced_time_of_arrival_partials )

BOOST_AUTO_TEST_CASE( testDifferencedTimeOfArrivalDirectPartial )
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

    // Define link ends for observations.
    LinkEnds linkEnds1;
    linkEnds1[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    linkEnds1[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    LinkEnds linkEnds2;
    linkEnds2[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-43" );
    linkEnds2[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    LinkEnds differencedLinkEnds;
    differencedLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    differencedLinkEnds[ receiver2 ] = std::make_pair< std::string, std::string >( "Earth", "DSS-43" );
    differencedLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back(
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    // Create observation settings
    std::shared_ptr< ObservationModelSettings > rangeObservableSettings1 = std::make_shared< ObservationModelSettings >(
            one_way_range,
            linkEnds1,
            lightTimeCorrectionSettings );

    std::shared_ptr< ObservationModelSettings > rangeObservableSettings2 = std::make_shared< ObservationModelSettings >(
            one_way_range,
            linkEnds2,
            lightTimeCorrectionSettings );

    std::shared_ptr< ObservationModelSettings > differencedObservableSettings = std::make_shared< ObservationModelSettings >(
            differenced_time_of_arrival,
            differencedLinkEnds,
            lightTimeCorrectionSettings );

    // Create observation model.
    std::shared_ptr< ObservationModel< 1, double, Time > > firstRangeObservationModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( rangeObservableSettings1, bodies );
    std::shared_ptr< ObservationModel< 1, double, Time > > secondRangeObservationModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( rangeObservableSettings2, bodies );
    std::shared_ptr< ObservationModel< 1, double, Time > > differencedTimeOfArrivalObservationModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( differencedObservableSettings, bodies );

    for( int parameterTest = 0; parameterTest < 1; parameterTest++ )
    {
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                "Mars",
                spice_interface::getBodyCartesianStateAtEpoch( "Mars", "SSB", bodies.getFrameOrientation( ), "None", initialEphemerisTime ),
                "SSB",
                bodies.getFrameOrientation( ) ) );

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate( parameterNames, bodies );

        std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
                   std::shared_ptr< PositionPartialScaling > >
                firstRangePartials = ObservationPartialCreator< 1, double, Time >::createObservationPartials(
                        firstRangeObservationModel, bodies, parametersToEstimate );
        std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
                   std::shared_ptr< PositionPartialScaling > >
                secondRangePartials = ObservationPartialCreator< 1, double, Time >::createObservationPartials(
                        secondRangeObservationModel, bodies, parametersToEstimate );
        std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
                   std::shared_ptr< PositionPartialScaling > >
                differencedTimeOfArrivalPartials = ObservationPartialCreator< 1, double, Time >::createObservationPartials(
                        differencedTimeOfArrivalObservationModel, bodies, parametersToEstimate );

        std::cout << "Partials sizes: " << differencedTimeOfArrivalPartials.first.size( ) << " " << firstRangePartials.first.size( ) << " "
                  << secondRangePartials.first.size( ) << " " << std::endl;

        // Compute observation separately with two functions.
        double receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
        std::vector< double > linkEndTimesDifferenced, linkEndTimesRange1, linkEndTimesRange2;
        std::vector< Eigen::Vector6d > linkEndStatesDifferenced, linkEndStatesRange1, linkEndStatesRange2;

        Eigen::VectorXd differencedTimeOfArrival = differencedTimeOfArrivalObservationModel->computeObservationsWithLinkEndData(
                receiverObservationTime, receiver, linkEndTimesDifferenced, linkEndStatesDifferenced );
        differencedTimeOfArrivalPartials.second->update(
                linkEndStatesDifferenced, linkEndTimesDifferenced, receiver, differencedTimeOfArrival );
        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > differencedTimeOfArrivalPartialValues =
                differencedTimeOfArrivalPartials.first.begin( )->second->calculatePartial(
                        linkEndStatesDifferenced, linkEndTimesDifferenced, receiver );

        Eigen::VectorXd firstRange = firstRangeObservationModel->computeObservationsWithLinkEndData(
                receiverObservationTime, receiver, linkEndTimesRange1, linkEndStatesRange1 );
        firstRangePartials.second->update( linkEndStatesRange1, linkEndTimesRange1, receiver, firstRange );
        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > firstRangePartialValues =
                firstRangePartials.first.begin( )->second->calculatePartial( linkEndStatesRange1, linkEndTimesRange1, receiver );

        Eigen::VectorXd secondRange = secondRangeObservationModel->computeObservationsWithLinkEndData(
                linkEndTimesRange1.at( 0 ), transmitter, linkEndTimesRange2, linkEndStatesRange2 );
        secondRangePartials.second->update( linkEndStatesRange2, linkEndTimesRange2, transmitter, secondRange );
        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > secondRangePartialValues =
                secondRangePartials.first.begin( )->second->calculatePartial( linkEndStatesRange2, linkEndTimesRange2, transmitter );

        BOOST_CHECK_CLOSE_FRACTION( firstRangePartialValues.at( 0 ).second, differencedTimeOfArrivalPartialValues.at( 0 ).second,
                                           ( 4.0 * std::numeric_limits< double >::epsilon( ) ));
        BOOST_CHECK_CLOSE_FRACTION( secondRangePartialValues.at( 0 ).second, differencedTimeOfArrivalPartialValues.at( 0 ).second,
                                           ( 4.0 * std::numeric_limits< double >::epsilon( ) ));

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( 
                firstRangePartialValues.at( 0 ).first,
                ( differencedTimeOfArrivalPartialValues.at( 0 ).first * physical_constants::SPEED_OF_LIGHT ),
                ( 4.0 * std::numeric_limits< double >::epsilon( ) ));
        
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( 
                secondRangePartialValues.at( 0 ).first,
                ( -differencedTimeOfArrivalPartialValues.at( 1 ).first * physical_constants::SPEED_OF_LIGHT ), 
                ( 4.0 * std::numeric_limits< double >::epsilon( ) ));
//        std::cout << std::setprecision( 15 ) << std::endl;
//        std::cout << firstRangePartialValues.at( 0 ).first << std::endl;
//        std::cout << differencedTimeOfArrivalPartialValues.at( 0 ).first * physical_constants::SPEED_OF_LIGHT << std::endl;
//
//        std::cout << secondRangePartialValues.at( 0 ).first << std::endl;
//        std::cout << -differencedTimeOfArrivalPartialValues.at( 1 ).first * physical_constants::SPEED_OF_LIGHT << std::endl;

        std::cout << "Done" << std::endl;
    }

}
//! Test partial derivatives of one-way differenced range observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testDifferencedTimeOfArrival )
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
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];
        linkEnds[ receiver2 ] = groundStations[ 2 ];

        // Generate one-way differenced range model
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::shared_ptr< ObservationModel< 1 > > oneWayDifferencedRangeModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        std::make_shared< observation_models::ObservationModelSettings >(
                                observation_models::differenced_time_of_arrival,
                                linkEnds,
                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) ),
                        bodies );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

        testObservationPartials< 1 >( oneWayDifferencedRangeModel,
                                      bodies,
                                      fullEstimatableParameterSet,
                                      linkEnds,
                                      differenced_time_of_arrival,
                                      1.0E-4,
                                      true,
                                      true,
                                      1000.0,
                                      parameterPerturbationMultipliers,
                                      getAveragedDopplerAncilliarySettings( 60.0 ) );
    }

//    // Test partials with real ephemerides (without test of position partials)
//    {
//        // Create environment
//        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );
//
//        // Set link ends for observation model
//        LinkEnds linkEnds;
//        linkEnds[ transmitter ] = groundStations[ 1 ];
//        linkEnds[ receiver ] = groundStations[ 0 ];
//
//        // Generate one-way range model
//        std::vector< std::string > perturbingBodies;
//        perturbingBodies.push_back( "Earth" );
//        std::shared_ptr< ObservationModel< 1 > > oneWayDifferencedRangeModel =
//                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
//                        std::make_shared< observation_models::ObservationModelSettings >(
//                                observation_models::one_way_differenced_range,
//                                linkEnds,
//                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) ),
//                        bodies );
//
//        // Create parameter objects.
//        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );
//
//        testObservationPartials< 1 >( oneWayDifferencedRangeModel,
//                                      bodies,
//                                      fullEstimatableParameterSet,
//                                      linkEnds,
//                                      one_way_differenced_range,
//                                      1.0E-4,
//                                      false,
//                                      true,
//                                      1000.0,
//                                      parameterPerturbationMultipliers,
//                                      getAveragedDopplerAncilliarySettings( 60.0 ) );
//    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
