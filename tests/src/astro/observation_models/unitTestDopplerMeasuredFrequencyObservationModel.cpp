/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
//
//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>


#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/astro/basic_astro/unitConversions.h"

//namespace tudat
//{
//namespace unit_tests
//{

using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::coordinate_conversions;
using namespace tudat::unit_conversions;


//BOOST_AUTO_TEST_SUITE(test_doppler_measured_frequency)

//BOOST_AUTO_TEST_CASE(testSimpleCase)

int main()
{
    // Load Spice kernels
    std::string kernelsPath = paths::getSpiceKernelPath();
    spice_interface::loadStandardSpiceKernels();

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mars" };

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings = getDefaultBodySettings(bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer);

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies(defaultBodySettings);

    // Create ground stations
    const Eigen::Vector3d stationCartesianPosition(1917032.190, 6029782.349, -801376.113);
    createGroundStation(bodies.at("Earth"), "Station1", stationCartesianPosition, cartesian_position);

    // // Set station with unrealistic position to force stronger proper time effect
    // const Eigen::Vector3d stationCartesianPosition2(4324532.0, 157372.0, -9292843.0);
    // createGroundStation(bodies.at("Earth"), "Station2", stationCartesianPosition2, cartesian_position);

    // Create Spacecraft
    Eigen::Vector6d spacecraftOrbitalElements;
    spacecraftOrbitalElements(semiMajorAxisIndex) = 10000.0E3;
    spacecraftOrbitalElements(eccentricityIndex) = 0.33;
    spacecraftOrbitalElements(inclinationIndex) = convertDegreesToRadians(65.3);
    spacecraftOrbitalElements(argumentOfPeriapsisIndex) = convertDegreesToRadians(235.7);
    spacecraftOrbitalElements(longitudeOfAscendingNodeIndex) = convertDegreesToRadians(23.4);
    spacecraftOrbitalElements(trueAnomalyIndex) = convertDegreesToRadians(0.0);
    double earthGravitationalParameter = bodies.at("Earth")->getGravityFieldModel()->getGravitationalParameter();

    bodies.createEmptyBody("Spacecraft");;
    bodies.at("Spacecraft")->setEphemeris(createBodyEphemeris(std::make_shared< KeplerEphemerisSettings >(spacecraftOrbitalElements, 0.0, earthGravitationalParameter, "Earth"), "Spacecraft"));
    bodies.processBodyFrameDefinitions();


    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[transmitter] = std::make_pair< std::string, std::string >("Earth", "");
    linkEnds[retransmitter] = std::make_pair< std::string, std::string >("Mars", "");
    linkEnds[receiver] = std::make_pair< std::string, std::string >("Earth", "");


    LinkEnds uplinkLinkEnds;
    uplinkLinkEnds[transmitter] = std::make_pair< std::string, std::string >("Earth", "");
    uplinkLinkEnds[receiver] = std::make_pair< std::string, std::string >("Mars", "");

    LinkEnds downlinkLinkEnds;
    downlinkLinkEnds[transmitter] = std::make_pair< std::string, std::string >("Mars", "");
    downlinkLinkEnds[receiver] = std::make_pair< std::string, std::string >("Earth", "");

    // Create observation settings

    std::shared_ptr< TwoWayDopplerObservationModel< double, double> > dopplerFrequencyObservationModel =
        std::dynamic_pointer_cast<TwoWayDopplerObservationModel< double, double>>(
            ObservationModelCreator< 1, double, double>::createObservationModel(
                std::make_shared< ObservationModelSettings >(doppler_measured_frequency, linkEnds), bodies));

    // std::shared_ptr< ObservationModel< 1, double, double > > uplinkDopplerObservationModel = ObservationModelCreator< 1, double, double>::createObservationModel(
    //     std::make_shared< ObservationModelSettings >(one_way_doppler, uplinkLinkEnds), bodies);
    // std::shared_ptr< ObservationModel< 1, double, double > > downlinkDopplerObservationModel = ObservationModelCreator< 1, double, double>::createObservationModel(std::make_shared< ObservationModelSettings >(one_way_doppler, downlinkLinkEnds), bodies);


    // // Creare independent light time calculator objects
    // std::shared_ptr< LightTimeCalculator< double, double > > uplinkLightTimeCalculator =
    //     createLightTimeCalculator(linkEnds, transmitter, retransmitter, bodies);
    // std::shared_ptr< LightTimeCalculator< double, double > > downlinkLightTimeCalculator =createLightTimeCalculator(linkEnds, retransmitter, receiver, bodies);

    // Test observable for both fixed link ends
    for (unsigned testCase = 0; testCase < 3; testCase++)
    {

        double observationTime = (finalEphemerisTime + initialEphemerisTime) / 2.0;
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;

        // Define link end
        LinkEndType referenceLinkEnd = receiver;
      
        // Compute observables
        double dopplerObservable = dopplerFrequencyObservationModel->computeObservationsWithLinkEndData(
            observationTime, referenceLinkEnd, linkEndTimes, linkEndStates)(0); 

        std::cout << "TEST: Doppler observable: " << dopplerObservable << "std::endl";
    }
    // std::dynamic_pointer_cast<OneWayDopplerObservationModel< double, double>>(
    //     uplinkDopplerObservationModel)->setNormalizeWithSpeedOfLight(0);
    // std::dynamic_pointer_cast<OneWayDopplerObservationModel< double, double>>(
    //     downlinkDopplerObservationModel)->setNormalizeWithSpeedOfLight(0);



}

//
//BOOST_AUTO_TEST_SUITE_END()
//
//}
//
//}
