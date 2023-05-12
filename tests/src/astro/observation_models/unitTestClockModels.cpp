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


#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

namespace tudat
{
namespace unit_tests
{


using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ground_stations;
using namespace tudat::system_models;

BOOST_AUTO_TEST_SUITE( test_clock_range_influence )

BOOST_AUTO_TEST_CASE( test_ClockRangeInfluence )
{
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 1.2E7;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define celestial bodies in simulation.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames );
    SystemOfBodies bodyMap = createSystemOfBodies( bodySettings );
    bodyMap.createEmptyBody( "LAGEOS" );
    std::shared_ptr< Body > lageos = bodyMap.at( "LAGEOS" );

    // Create dummy constant state of satellite.
    Eigen::Vector6d lageosKeplerianElements;
    lageosKeplerianElements[ semiMajorAxisIndex ] = 8000.0E3;
    lageosKeplerianElements[ eccentricityIndex ] = 0.0044;
    lageosKeplerianElements[ inclinationIndex ] = 109.89 * mathematical_constants::PI / 180.0;
    lageosKeplerianElements[ argumentOfPeriapsisIndex ] = 259.35 * mathematical_constants::PI / 180.0;
    lageosKeplerianElements[ longitudeOfAscendingNodeIndex ] = 31.56 * mathematical_constants::PI / 180.0;
    lageosKeplerianElements[ trueAnomalyIndex ] = 1.0;
    Eigen::Vector6d dummyLageosState = convertKeplerianToCartesianElements(
                lageosKeplerianElements, getBodyGravitationalParameter("Earth" ) );
    lageos->setEphemeris( std::make_shared< ConstantEphemeris >(
                              dummyLageosState, "Earth", "ECLIPJ2000" ) );

    // Create ground stations
    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Earth", "Graz" ) ] =
            ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );
    createGroundStations( bodyMap, groundStationsToCreate );

    std::shared_ptr< GroundStation > grazStation = bodyMap.at( "Earth" )->getGroundStation( "Graz" );

    // Create satellite time error arcs.
    std::vector< Time > satelliteClockErrorArcTimes;
    double satelliteSingleArcLength = 0.5E5;
    Time currentTime = Time( initialEphemerisTime );
    while( currentTime < Time( finalEphemerisTime ) )
    {
        satelliteClockErrorArcTimes.push_back( currentTime );
        currentTime += satelliteSingleArcLength;
    }
    satelliteClockErrorArcTimes.push_back( finalEphemerisTime + satelliteSingleArcLength );

    // Set errors for all arcs of satellite.
    std::vector< double > lageosArcPolynomialErrors;
    lageosArcPolynomialErrors.push_back( 0.0 );
    lageosArcPolynomialErrors.push_back( 0.0 );
    lageosArcPolynomialErrors.push_back( 0.0 );

    // Create and set timing system of satellite.
    std::shared_ptr< TimingSystem > lageosTimingSystem = std::make_shared< TimingSystem >( satelliteClockErrorArcTimes, lageosArcPolynomialErrors );
    std::shared_ptr< VehicleSystems > lageosSystemHardware = std::make_shared< VehicleSystems >( );
    lageosSystemHardware->setTimingSystem( lageosTimingSystem );
    lageos->setVehicleSystems( lageosSystemHardware );

    // Create graz time error arcs.
    std::vector< Time > grazClockErrorArcTimes;
    double grazSingleArcLength = 1.2E5;
    currentTime = Time( initialEphemerisTime );
    while( currentTime < Time( finalEphemerisTime ) )
    {
        grazClockErrorArcTimes.push_back( currentTime );
        currentTime += grazSingleArcLength;
    }
    grazClockErrorArcTimes.push_back( finalEphemerisTime + grazSingleArcLength );

    // Set errors for all arcs of ground station.
    std::vector< double > grazArcPolynomialErrors;
    grazArcPolynomialErrors.push_back( 0.0 );
    grazArcPolynomialErrors.push_back( 0.0 );
    grazArcPolynomialErrors.push_back( 0.0 );

    // Create and set timing system of graz.
    std::shared_ptr< TimingSystem > grazTimingSystem = std::make_shared< TimingSystem >( grazClockErrorArcTimes, grazArcPolynomialErrors );
    grazStation->setTimingSystem( grazTimingSystem );

    // Set link ends for test model.
    LinkEnds linkEndList;
    linkEndList[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Graz" ) );
    linkEndList[ receiver ] = LinkEndId( std::make_pair( "LAGEOS", "" ) );

    // Link timing systems to observation biases
    std::vector< LinkEnds > allLinkEndsList;
    allLinkEndsList.push_back( linkEndList );
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ observation_models::one_way_range ] = allLinkEndsList;


    // Create range model.
    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
    biasSettingsList.push_back( std::make_shared< TiminigSystemBiasSettings >( "LAGEOS", "" ) );
    biasSettingsList.push_back( std::make_shared< TiminigSystemBiasSettings >( "Earth", "Graz" ) );
    std::shared_ptr< ObservationModelSettings > rangeSettings = std::make_shared< ObservationModelSettings >(
            one_way_range, linkEndList, nullptr, std::make_shared< MultipleObservationBiasSettings >( biasSettingsList ) );

    std::shared_ptr< ObservationModel< 1, double, double > > rangeModel = ObservationModelCreator< 1, double, double >::createObservationModel(
            rangeSettings, bodyMap );

    // Get range without timing errors.
    double arcStartTime = satelliteClockErrorArcTimes[ 10 ];
    double timeIntoArc = 2.0E4;
    double testTime = arcStartTime + timeIntoArc;
    std::cout<<"No time error, pre"<<std::endl;
    Eigen::VectorXd rangeWithoutTimeError = rangeModel->computeObservations( testTime, receiver );
    std::cout<<"No time error, post"<<std::endl;

    // Set polynomial offsets of satellite clock.
    std::map< int, double > newClockCorrections;
    newClockCorrections[ 0 ] = 1.0E-3;
    newClockCorrections[ 1 ] = -2.0E-7;
    newClockCorrections[ 2 ] = 5.0E-13;
    lageosTimingSystem->setGlobalPolynomialClockCorrections( newClockCorrections );

    // Get range with satellite timing errors.
    std::cout<<"With satellite error, pre"<<std::endl;
    Eigen::VectorXd rangeWithSatelliteTimeError = rangeModel->computeObservations( testTime, receiver );
    std::cout<<"With satellite error, post"<<std::endl;

    // Manually calculate satellite timing error
    double expectedSatelliteTimingError = newClockCorrections[ 0 ] + newClockCorrections[ 1 ] * timeIntoArc +
            newClockCorrections[ 2 ] * timeIntoArc * timeIntoArc;

    BOOST_CHECK_CLOSE_FRACTION( expectedSatelliteTimingError * physical_constants::SPEED_OF_LIGHT,
                                ( rangeWithSatelliteTimeError - rangeWithoutTimeError ).x( ),
                                ( 2.0 * rangeWithoutTimeError.x( ) * std::numeric_limits< double >::epsilon( ) ) );

    // Reset satellite clock errors to zero.
    newClockCorrections[ 0 ] = 0.0;
    newClockCorrections[ 1 ] = 0.0;
    newClockCorrections[ 2 ] = 0.0;
    lageosTimingSystem->setGlobalPolynomialClockCorrections( newClockCorrections );

    // Set station clock errors.
    newClockCorrections[ 0 ] = -2.3E-4;
    newClockCorrections[ 1 ] = 5.2E-9;
    newClockCorrections[ 2 ] = 3.0E-13;
    grazTimingSystem->setGlobalPolynomialClockCorrections( newClockCorrections );

    std::pair< Time, int > timeIntoCurrentStationArc = grazTimingSystem->getTimeIntoCurrentArcAndArcIndex(
                testTime - rangeWithoutTimeError.x( ) / physical_constants::SPEED_OF_LIGHT );

    std::cout<<"With station error, pre"<<std::endl;
    Eigen::VectorXd rangeWithStationTimeError = rangeModel->computeObservations(
                testTime, receiver );
    std::cout<<"With station error, post"<<std::endl;

    double expectedStationTimingError = newClockCorrections[ 0 ] + newClockCorrections[ 1 ] * static_cast< double >( timeIntoCurrentStationArc.first ) +
            newClockCorrections[ 2 ] * static_cast< double >( timeIntoCurrentStationArc.first ) *
            static_cast< double >( timeIntoCurrentStationArc.first );

    BOOST_CHECK_CLOSE_FRACTION( -expectedStationTimingError * physical_constants::SPEED_OF_LIGHT,
                                ( rangeWithStationTimeError - rangeWithoutTimeError ).x( ),
                                ( 2.0 * rangeWithoutTimeError.x( ) * std::numeric_limits< double >::epsilon( ) ) );


    newClockCorrections[ 0 ] = 1.0E-3;
    newClockCorrections[ 1 ] = -2.0E-7;
    newClockCorrections[ 2 ] = 5.0E-13;
    lageosTimingSystem->setGlobalPolynomialClockCorrections( newClockCorrections );

    Eigen::VectorXd rangeWithStationAndSatelliteTimeError = rangeModel->computeObservations(
                testTime, receiver );

    BOOST_CHECK_CLOSE_FRACTION( ( expectedSatelliteTimingError - expectedStationTimingError ) * physical_constants::SPEED_OF_LIGHT,
                                ( rangeWithStationAndSatelliteTimeError - rangeWithoutTimeError ).x( ),
                                ( 2.0 * rangeWithoutTimeError.x( ) * std::numeric_limits< double >::epsilon( ) ) );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
