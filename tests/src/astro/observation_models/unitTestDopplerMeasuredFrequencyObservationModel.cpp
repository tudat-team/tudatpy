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
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/io/readTrackingTxtFile.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::coordinate_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::input_output;



BOOST_AUTO_TEST_SUITE(test_doppler_measured_frequency)

BOOST_AUTO_TEST_CASE(testJuiceMeasuredFrequency)
{
    // Load Spice kernels
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( tudat::paths::getTudatTestDataPath( ) + "juice_orbc_000074_230414_310721_v01.bsp" );

    Eigen::VectorXd originalResidual;
    for( int testType = 0; testType < 2 ; testType++ )
    {
        // Define bodies to use.
        std::vector<std::string> bodiesToCreate = { "Earth", "Moon", "Sun", "Jupiter" };
        std::string globalFrameOrigin = "SSB";
        std::string globalFrameOrientation = "J2000";

        // Create bodies settings needed in simulation
        BodyListSettings bodySettings = getDefaultBodySettings(
            bodiesToCreate, globalFrameOrigin, globalFrameOrientation );
        bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
        bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
            basic_astrodynamics::iau_2006, globalFrameOrientation );
        bodySettings.at( "Earth" )->bodyDeformationSettings.push_back( iers2010TidalBodyShapeDeformation( ) );

        std::shared_ptr<GroundStationSettings> nnorciaSettings = std::make_shared<GroundStationSettings>(
            "NWNORCIA", getCombinedApproximateGroundStationPositions( ).at( "NWNORCIA" ));
        nnorciaSettings->addStationMotionSettings(
            std::make_shared<LinearGroundStationMotionSettings>(
                ( Eigen::Vector3d( ) << -45.00, 10.00, 47.00 ).finished( ) / 1.0E3 / physical_constants::JULIAN_YEAR,
                0.0 ));

        std::shared_ptr<GroundStationSettings> yarragadeeSettings = std::make_shared<GroundStationSettings>(
            "YARRAGAD", getCombinedApproximateGroundStationPositions( ).at( "YARRAGAD" ));
        yarragadeeSettings->addStationMotionSettings(
            std::make_shared<LinearGroundStationMotionSettings>(
                ( Eigen::Vector3d( ) << -47.45, 9.12, 51.76 ).finished( ) / 1.0E3 / physical_constants::JULIAN_YEAR, 0.0 ));

        bodySettings.at( "Earth" )->groundStationSettings.push_back( nnorciaSettings );
        bodySettings.at( "Earth" )->groundStationSettings.push_back( yarragadeeSettings );


        // Create Spacecraft
        const std::string spacecraftName = "JUICE";
        bodiesToCreate.push_back( spacecraftName );
        bodySettings.addSettings( spacecraftName );
        bodySettings.get( spacecraftName )->ephemerisSettings = directSpiceEphemerisSettings( "Earth", "J2000", false );

        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Set turnaround ratios in spacecraft (ground station)
        std::shared_ptr<system_models::VehicleSystems> vehicleSystems = std::make_shared<system_models::VehicleSystems>( );
        vehicleSystems->setTransponderTurnaroundRatio( &getDsnDefaultTurnaroundRatios );
        bodies.at( spacecraftName )->setVehicleSystems( vehicleSystems );

        bodies.processBodyFrameDefinitions( );


        // Define link ends for observations.
        LinkEnds linkEnds;
        linkEnds[ transmitter ] =
            std::make_pair<std::string, std::string>( "Earth", static_cast<std::string>( "NWNORCIA" ));
        linkEnds[ retransmitter ] =
            std::make_pair<std::string, std::string>( static_cast<std::string>(spacecraftName), "" );
        linkEnds[ receiver ] = std::make_pair<std::string, std::string>( "Earth", static_cast<std::string>( "YARRAGAD" ));


        std::shared_ptr<ground_stations::StationFrequencyInterpolator> transmittingFrequencyCalculator =
            std::make_shared<ground_stations::ConstantFrequencyInterpolator>( 7180127320 );

        bodies.at( "Earth" )->getGroundStation( "NWNORCIA" )->setTransmittingFrequencyCalculator(
            transmittingFrequencyCalculator );

        // Create observation settings
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList;
        lightTimeCorrectionsList.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >( { "Sun", "Moon", "Earth" } ) ) );
        std::shared_ptr<DopplerMeasuredFrequencyObservationModel<double, Time> > dopplerFrequencyObservationModel =
            std::dynamic_pointer_cast<DopplerMeasuredFrequencyObservationModel<double, Time>>(
                ObservationModelCreator<1, double, Time>::createObservationModel(
                    std::make_shared<ObservationModelSettings>( doppler_measured_frequency, linkEnds, lightTimeCorrectionsList   ), bodies ));
        
        // Define link end
        LinkEndType referenceLinkEnd = receiver;

        // Ancillary Settings
        std::shared_ptr<observation_models::ObservationCollection<double, Time> > observationCollection;
        std::string juiceDataFile = tudat::paths::getTudatTestDataPath( ) + "Fdets.jui2024.08.20.Yg.r2i.txt";
        std::vector< std::string > columnTypes =  { "scan_number", "utc_datetime_string", "signal_to_noise_ratio", "normalised_spectral_max",
                                                    "doppler_measured_frequency_hz", "doppler_noise_hz" };
        if( testType == 0 )
        {

            std::shared_ptr<TrackingTxtFileContents> fdetsFileContents = readFdetsFile(
                juiceDataFile, columnTypes );
            fdetsFileContents->addMetaData( TrackingDataType::receiving_station_name, "YARRAGAD" );
            fdetsFileContents->addMetaData( TrackingDataType::transmitting_station_name, "NWNORCIA" );
            fdetsFileContents->addMetaData( TrackingDataType::doppler_base_frequency, 8422.49E6 );

            observationCollection =
                createTrackingTxtFileObservationCollection< double, Time >(
                    fdetsFileContents, "JUICE" );
        }
        else
        {
            observationCollection = createFdetsObservedObservationCollectionFromFile(
                juiceDataFile, 8422.49E6, columnTypes, "JUICE", "NWNORCIA", "YARRAGAD", x_band, x_band );
        }

        auto observationTimes = observationCollection->getConcatenatedObservationTimes( );
        auto observations = observationCollection->getConcatenatedObservations( );

        // Compute observables
        std::vector<double> linkEndTimes;
        std::vector<Eigen::Vector6d> linkEndStates;

        Eigen::VectorXd observableVector = Eigen::VectorXd::Zero( observationTimes.size( ) );
        Eigen::VectorXd residualVector = Eigen::VectorXd::Zero( observationTimes.size( ) );
        std::shared_ptr<ObservationAncilliarySimulationSettings> ancillarySettings =
            std::make_shared<ObservationAncilliarySimulationSettings>( );
        ancillarySettings->setAncilliaryDoubleVectorData( frequency_bands, { x_band, x_band } );
        for( unsigned int i = 0; i < observationTimes.size( ); i++ )
        {
            double dopplerObservable = dopplerFrequencyObservationModel->computeObservationsWithLinkEndData(
                observationTimes.at( i ), referenceLinkEnd, linkEndTimes, linkEndStates, ancillarySettings )( 0 );
            observableVector( i ) = dopplerObservable;
            residualVector( i ) = dopplerObservable - observations( i );
        }

        BOOST_CHECK_SMALL( linear_algebra::getVectorEntryRootMeanSquare( residualVector.segment( 0, 7000 ) ), 2.0 );
        BOOST_CHECK_SMALL( std::fabs( linear_algebra::getVectorEntryMean( residualVector.segment( 0, 7000 ) ) ), 2.0 );

        if( testType == 0 )
        {
            originalResidual = residualVector;
        }
        else
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( originalResidual, residualVector, std::numeric_limits< double >::epsilon( ) );
        }
//    std::cout<<linear_algebra::getVectorEntryRootMeanSquare( residualVector.segment( 0, 7000 ) )<<std::endl;
//    std::cout<<linear_algebra::getVectorEntryMean( residualVector.segment( 0, 7000 ) )<<std::endl;

//
//    input_output::writeMatrixToFile( observableVector, "pride_doppler.dat", 16 );
//    input_output::writeMatrixToFile( residualVector, "pride_residual.dat", 16 );
//    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( utilities::staticCastVector< double, Time >( observationTimes ) ), "pride_times.dat", 16 );
    }

}


BOOST_AUTO_TEST_SUITE_END()

}

}
