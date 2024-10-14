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
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/io/readTrackingTxtFile.h"

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

namespace tio = tudat::input_output;


//BOOST_AUTO_TEST_SUITE(test_doppler_measured_frequency)

//BOOST_AUTO_TEST_CASE(testSimpleCase)


const static std::string juiceDataFile = "/home/simon/lib/tudat-bundle/tudatpy/examples/estimation/data/Fdets.jui2023.09.14.Hb.r2i.txt";

std::shared_ptr<tio::TrackingTxtFileContents> readJuiceFdetsFile(const std::string& fileName)
{
    std::vector<std::string>
        columnTypes({ "utc_datetime_string", "signal_to_noise_ratio", "normalised_spectral_max", "doppler_measured_frequency_hz", "doppler_noise_hz", });

    auto rawFileContents = tio::createTrackingTxtFileContents(fileName, columnTypes, '#', ", \t");
    rawFileContents->addMetaData(tio::TrackingDataType::file_name, "JUICE Fdets Test File");
    return rawFileContents;
}


int main()
{
    // Load Spice kernels
    std::string kernelsPath = paths::getSpiceKernelPath();
    spice_interface::loadStandardSpiceKernels();

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Jupiter" };
    std::string globalFrameOrigin = "Sun";
    std::string globalFrameOrientation = "J2000";


    // Specify simulation time settings
    std::string initialEphemerisDateString = "2023-09-14T06:18:05.000"; // FDETS FILE START
    double initialEphemerisTime = convertDateStringToEphemerisTime(initialEphemerisDateString);
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    double initialTimeEnvironment = initialEphemerisTime - buffer;
    double finalTimeEnvironment = finalEphemerisTime + buffer;

    // Create bodies settings needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings(
        bodiesToCreate, initialTimeEnvironment, finalTimeEnvironment, globalFrameOrigin, globalFrameOrientation);

    // Create Spacecraft
    const std::string spacecraftName = "JUICE";
    bodiesToCreate.push_back(spacecraftName);
    bodySettings.addSettings(spacecraftName);
    bodySettings.get(spacecraftName)->ephemerisSettings = directSpiceEphemerisSettings("Sun", "J2000", false);

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies(bodySettings);

    // Create ground station
    std::string stationName = "HOBART12";
    const Eigen::Vector3d stationCartesianPosition = Eigen::Vector3d(-3949990.106, 2522421.118, -4311708.734);
    createGroundStation(bodies.at("Earth"), stationName, stationCartesianPosition, cartesian_position);


    // Set turnaround ratios in spacecraft (ground station)
    std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared< system_models::VehicleSystems >();
    vehicleSystems->setTransponderTurnaroundRatio(&getDsnDefaultTurnaroundRatios);
    bodies.at(spacecraftName)->setVehicleSystems(vehicleSystems);

    bodies.processBodyFrameDefinitions();


    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[transmitter] = std::make_pair< std::string, std::string >("Earth", static_cast<std::string>(stationName));
    linkEnds[retransmitter] = std::make_pair< std::string, std::string >(static_cast<std::string>(spacecraftName), "");
    linkEnds[receiver] = std::make_pair< std::string, std::string >("Earth", static_cast<std::string>(stationName));


    std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator =
        std::make_shared< ground_stations::ConstantFrequencyInterpolator >(7.18E9);

    bodies.at("Earth")->getGroundStation(stationName)->setTransmittingFrequencyCalculator(transmittingFrequencyCalculator);

    // Create observation settings
    std::shared_ptr< DopplerMeasuredFrequencyObservationModel< double, double> > dopplerFrequencyObservationModel =
        std::dynamic_pointer_cast<DopplerMeasuredFrequencyObservationModel< double, double>>(
            ObservationModelCreator< 1, double, double>::createObservationModel(
                std::make_shared< ObservationModelSettings >(doppler_measured_frequency, linkEnds), bodies));

    // Test observable for both fixed link ends

    double observationTime = initialEphemerisTime;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;

    // Define link end
    LinkEndType referenceLinkEnd = receiver;

    // Ancillary Settings
    std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings =
        std::make_shared< ObservationAncilliarySimulationSettings >();
    ancillarySettings->setAncilliaryDoubleVectorData(frequency_bands, { x_band, x_band });

    // Compute observables
    double dopplerObservable = dopplerFrequencyObservationModel->computeObservationsWithLinkEndData(
        observationTime, referenceLinkEnd, linkEndTimes, linkEndStates, ancillarySettings)(0);

    std::cout << "TEST: Doppler observable: " << dopplerObservable - 8.435E9 << std::endl;

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
