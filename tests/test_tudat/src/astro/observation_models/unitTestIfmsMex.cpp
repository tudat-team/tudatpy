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
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

#include "tudat/io/preProcessIfmsFile.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/simulation/estimation_setup/createObservationDataset.h"
#include "tudat/simulation/estimation_setup/compressDopplerObservationCollection.h"

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"

using namespace tudat::propagators;
using namespace tudat::estimatable_parameters;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::input_output;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::reference_frames;
using namespace tudat;

namespace tudat
{
namespace unit_tests
{

//! Starting the entire test suite
BOOST_AUTO_TEST_SUITE( test_ifms_observation );

//! Test reading of mars express IFMS files.
BOOST_AUTO_TEST_CASE( testIfmsObservationMex )
{
    double initialTimeEnvironment = DateTime( 2013, 12, 29, 0, 0, 0.0 ).epoch< double >( );
    double finalTimeEnvironment = DateTime( 2013, 12, 30, 0, 0, 0.0 ).epoch< double >( );

    /****************************************************************************************
     ************************** CREATE EVIRONMENT
     *****************************************************************************************/

    std::vector< std::string > ifmsFileNames;
    ifmsFileNames.push_back( paths::getTudatTestDataPath( ) +
                             "/estrack_n_way_doppler_observation_model/M32ICL3L02_D2S_133621904_00.TAB.txt" );

    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( paths::getTudatTestDataPath( ) +
                                             "/estrack_n_way_doppler_observation_model/MEX_ROB_130101_131231_001_shortened.bsp" );

    for( int i = 0; i < 1; i++ )
    {
        FrequencyBands currentReceptionBand = x_band;
        if( i < 2 )
        {
            currentReceptionBand = s_band;
        }

        // Create settings for default bodies
        std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Moon", "Mars" };
        std::string globalFrameOrigin = "SSB";
        std::string globalFrameOrientation = "J2000";
        BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, globalFrameOrigin, globalFrameOrientation );

        // Add high-accuracy Earth settings
        bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
        bodySettings.at( "Earth" )->rotationModelSettings =
                gcrsToItrsRotationModelSettings( basic_astrodynamics::iau_2006, globalFrameOrientation );
        bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );
        bodySettings.at( "Earth" )->bodyDeformationSettings.push_back( iers2010TidalBodyShapeDeformation( ) );
        bodySettings.at( "Earth" )->groundStationSettings = getRadioTelescopeStationSettings( );

        // Add spacecraft settings
        std::string spacecraftName = "MeX";
        std::string spacecraftCentralBody = "Mars";
        bodySettings.addSettings( spacecraftName );
        bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< DirectSpiceEphemerisSettings >( spacecraftCentralBody, globalFrameOrientation );

        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

        // Set turnaround ratios in spacecraft (ground station)
        std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared< system_models::VehicleSystems >( );
        vehicleSystems->setTransponderTurnaroundRatio( &getDsnDefaultTurnaroundRatios );
        bodies.at( "MeX" )->setVehicleSystems( vehicleSystems );

        /****************************************************************************************
         ************************** LOAD IFMS FILES
         *****************************************************************************************/

        auto trackingDataAndSupplementaryData = readIfmsFiles< long double, Time >(
                { ifmsFileNames.at( i ) },
                "MeX",
                "NWNORCIA",
                "Earth",
                true,
                true,
                std::vector< std::string >( { getFrequencyBandString( x_band ), getFrequencyBandString( currentReceptionBand ) } ),
                getFrequencyBandString( currentReceptionBand ),
                0.0 );

        setTrackingSupplementaryDataInBodies( bodies, trackingDataAndSupplementaryData.second );

        std::shared_ptr< observation_models::ObservationDataset< long double, Time > > observedUncompressedObservationDataset =
                createObservationDatasetFromTrackingData< long double, Time >( trackingDataAndSupplementaryData.first, bodies );

        std::shared_ptr< observation_models::ObservationDataset< long double, Time > > observedObservationDataset =
                createCompressedDopplerDataset( observedUncompressedObservationDataset, 60.0 );

        /****************************************************************************************
         ************************** CREATE OBSERVATION MODEL SETTINGS
         *****************************************************************************************/

        // Define observation model settings
        std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
        std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
        lightTimeCorrectionSettings.push_back( firstOrderRelativisticLightTimeCorrectionSettings( { "Sun" } ) );
        for( unsigned int setId = 0; setId < observedObservationDataset->getNumberOfObservationSets( ); ++setId )
        {
            const ObservationSetMetadata< long double, Time >& metadata = observedObservationDataset->getObservationSetMetadata( setId );
            if( metadata.observableType_ == observation_models::dsn_n_way_averaged_doppler )
            {
                observationModelSettingsList.push_back( observation_models::dsnNWayAveragedDopplerObservationSettings(
                        observedObservationDataset->getLinkDefinition( metadata.linkDefinitionId_ ),
                        lightTimeCorrectionSettings,
                        constantAbsoluteBias( Eigen::Vector1d::Zero( ) ),
                        std::make_shared< LightTimeConvergenceCriteria >( true ),
                        false ) );
            }
        }

        // Create observation simulator
        std::vector< std::shared_ptr< ObservationSimulatorBase< long double, Time > > > observationSimulators =
                createObservationSimulators< long double, Time >( observationModelSettingsList, bodies );

        /****************************************************************************************
         ************************** SIMULATE OBSERVATIONS AND COMPUTE RESIDUALS
         *****************************************************************************************/

        std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > observationSimulationSettings =
                getObservationSimulationSettingsFromObservationDataset( observedObservationDataset, bodies );
        std::shared_ptr< observation_models::ObservationDataset< long double, Time > > computedObservationDataset =
                simulateObservationDataset( observationSimulationSettings, observationSimulators, bodies );

        Eigen::Matrix< long double, Eigen::Dynamic, 1 > residualVector =
                observedObservationDataset->createEstimationFlattenedObservationData( ).getObservationVector( ) -
                computedObservationDataset->createEstimationFlattenedObservationData( ).getObservationVector( );
        double rmsResidual = linear_algebra::getVectorEntryRootMeanSquare( residualVector.cast< double >( ) );
        double meanResidual = linear_algebra::getVectorEntryMean( residualVector.cast< double >( ) );

        std::cout << residualVector.rows( ) << " " << rmsResidual << " " << meanResidual << std::endl;
        BOOST_CHECK_SMALL( rmsResidual, 3.5E-3 );
        BOOST_CHECK_SMALL( meanResidual, 1.0E-3 );

        //        input_output::writeMatrixToFile( observedObservationCollection->getObservationVector( ), "ifms_doppler_" +
        //        std::to_string( i ) + ".dat", 16 ); input_output::writeMatrixToFile( residualVector, "ifms_residuals_" +
        //        std::to_string( i ) + ".dat", 16 ); input_output::writeMatrixToFile(
        //            utilities::convertStlVectorToEigenVector(
        //                utilities::staticCastVector< double, Time >( observedObservationCollection->getConcatenatedTimeVector() ) ),
        //                "ifms_times_" + std::to_string( i ) + ".dat", 16 );
    }
}
// End test suite
BOOST_AUTO_TEST_SUITE_END( );

}  // namespace unit_tests

}  // namespace tudat
