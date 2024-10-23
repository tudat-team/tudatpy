/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/simulation/estimation_setup.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"

using namespace tudat::propagators;
using namespace tudat::estimatable_parameters;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::input_output;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::numerical_integrators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::reference_frames;

using namespace tudat;

int main( )
{
    double initialTimeEnvironment = DateTime( 2012, 03, 02, 0, 0, 0.0 ).epoch< double >();
    double finalTimeEnvironment = DateTime( 2012, 05, 29, 0, 0, 0.0 ).epoch< double >();

//    DateTime initialDateTime = basic_astrodynamics::getCalendarDateFromTime< double >( initialTimeEnvironment );
//    DateTime finalDateTime = basic_astrodynamics::getCalendarDateFromTime< double >( finalTimeEnvironment );

    /****************************************************************************************
     ************************** CREATE EVIRONMENT
     *****************************************************************************************/

    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( get_spice_kernels_path( ) + + "/moon_de440_200625.tf");
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/grail_v07.tf" );
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/gra_sclkscet_00013.tsc" );
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/gra_sclkscet_00014.tsc" );
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/grail_120301_120529_sci_v02.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/gra_rec_120402_120408.bc" );
    spice_interface::loadSpiceKernelInTudat( get_spice_kernels_path( ) + + "/moon_pa_de440_200625.bpc");

    // Create settings for default bodies
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon" };
    std::string globalFrameOrigin = "SSB";
    std::string globalFrameOrientation = "J2000";
    BodyListSettings bodySettings = getDefaultBodySettings(
            bodiesToCreate, initialTimeEnvironment, finalTimeEnvironment, globalFrameOrigin, globalFrameOrientation, 120.0 );

    // Add high-accuracy Earth settings
    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
            basic_astrodynamics::iau_2006, globalFrameOrientation,
            std::make_shared< interpolators::InterpolatorGenerationSettings< double > >(
                    interpolators::cubicSplineInterpolation( ), initialTimeEnvironment, finalTimeEnvironment, 3600.0 ),
            std::make_shared< interpolators::InterpolatorGenerationSettings< double > >(
                    interpolators::cubicSplineInterpolation( ), initialTimeEnvironment, finalTimeEnvironment, 3600.0 ),
            std::make_shared< interpolators::InterpolatorGenerationSettings< double > >(
                    interpolators::cubicSplineInterpolation( ), initialTimeEnvironment, finalTimeEnvironment, 60.0 ));
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );
    bodySettings.at( "Moon" )->rotationModelSettings = spiceRotationModelSettings(
            bodySettings.getFrameOrientation( ), "MOON_PA_DE440", "MOON_PA_DE440" );
    bodySettings.at( "Moon" )->gravityFieldSettings =
            std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( gggrx1200, 500 );
    std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
            bodySettings.at( "Moon" )->gravityFieldSettings )->resetAssociatedReferenceFrame( "MOON_PA_DE440" );
    bodySettings.at( "Moon" )->gravityFieldVariationSettings.push_back(
            fixedSingleDegreeLoveNumberGravityFieldVariationSettings( "Earth", 0.02405, 2 ) );
    bodySettings.at( "Moon" )->gravityFieldVariationSettings.push_back(
            fixedSingleDegreeLoveNumberGravityFieldVariationSettings( "Sun", 0.02405, 2 ) );
    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( "Earth" );

    std::vector< std::shared_ptr< PanelRadiosityModelSettings > > panelRadiosityModels;
    panelRadiosityModels.push_back(angleBasedThermalPanelRadiosityModelSettings( 95.0, 385.0, 0.95, "Sun" ) );
    panelRadiosityModels.push_back(albedoPanelRadiosityModelSettings( SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1, "Sun" ) );

    std::map< std::string, std::vector< std::string > > originalSourceToSourceOccultingBodies;
    originalSourceToSourceOccultingBodies[ "Sun" ].push_back( "Earth" );
    bodySettings.at( "Moon" )->radiationSourceModelSettings =
            extendedRadiationSourceModelSettingsWithOccultationMap(
                    panelRadiosityModels, { 4, 8, 12, 16  }, originalSourceToSourceOccultingBodies );

    // Add spacecraft settings
    std::string spacecraftName = "GRAIL-A";
    std::string spacecraftCentralBody = "Moon";
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings =
            std::make_shared< InterpolatedSpiceEphemerisSettings >(
                    initialTimeEnvironment, finalTimeEnvironment, 10.0, spacecraftCentralBody, globalFrameOrientation );

    // Create spacecraft
    bodySettings.at( spacecraftName )->constantMass = 150.0;

    // Create radiation pressure settings
    double referenceAreaRadiation = 5.0;
    double radiationPressureCoefficient = 1.5;
    std::map< std::string, std::vector< std::string > > sourceToTargetOccultingBodies;
    sourceToTargetOccultingBodies[ "Sun" ].push_back( "Moon" );
    bodySettings.at( spacecraftName )->radiationPressureTargetModelSettings = cannonballRadiationPressureTargetModelSettingsWithOccultationMap(
            referenceAreaRadiation, radiationPressureCoefficient, sourceToTargetOccultingBodies );
    bodySettings.at( spacecraftName )->rotationModelSettings = spiceRotationModelSettings(
            globalFrameOrientation, spacecraftName + "_SPACECRAFT", "" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );
    bodies.at( "GRAIL-A" )->getVehicleSystems( )->setReferencePointPosition(
            "Antenna", ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( ) );
    std::cout<<"Number of reference points: "<<bodies.at( "GRAIL-A" )->getVehicleSystems( )->getReferencePoints( ).size( )<<std::endl;

    /****************************************************************************************
     ************************** LOAD ODF FILES
     *****************************************************************************************/

    // Define ODF data paths
    std::string dataDirectory = "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_data/";
    std::vector< std::string > odfFiles = { "gralugf2012_097_0235smmmv1.odf" };

    // Laod raw ODF data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    for ( std::string odfFile : odfFiles )
    {
        rawOdfDataVector.push_back( std::make_shared< OdfRawFileContents >( dataDirectory + odfFile ) );
    }

    // Process ODF file data
    std::shared_ptr< ProcessedOdfFileContents< Time > > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents< Time > >( rawOdfDataVector, spacecraftName, true );
    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );

    // Create data structure that handles Observed Data in Tudat
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
            createCompressedDopplerCollection( observation_models::createOdfObservedObservationCollection< long double, Time >(
                    processedOdfFileContents, { dsn_n_way_averaged_doppler } ), 60.0, 10 );

    /****************************************************************************************
     ************************** CREATE OBSERVATION MODEL SETTINGS
     *****************************************************************************************/

    // Define observation model settings
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;

    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( firstOrderRelativisticLightTimeCorrectionSettings( { "Sun" } ) );
    std::vector< std::string > troposphericCorrectionFileNames =
            {"/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_data/grxlugf2012_092_2012_122.tro",
             "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_data/grxlugf2012_122_2012_153.tro"};
    std::vector< std::string > ionosphericCorrectionFileNames =
            {"/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_data/gralugf2012_092_2012_122.ion",
             "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_data/gralugf2012_122_2012_153.ion"};
    std::map< int, std::string > spacecraftNamePerSpacecraftId;
    spacecraftNamePerSpacecraftId[ 177 ] = "GRAIL-A";

    lightTimeCorrectionSettings.push_back( tabulatedTroposphericCorrectionSettings( troposphericCorrectionFileNames ) );
    lightTimeCorrectionSettings.push_back( tabulatedIonosphericCorrectionSettings( ionosphericCorrectionFileNames, spacecraftNamePerSpacecraftId ) );

    std::map < observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable =
            observedObservationCollection->getLinkEndsPerObservableType( );

    for ( auto it = linkEndsPerObservable.begin(); it != linkEndsPerObservable.end(); ++it )
    {
        for ( unsigned int i = 0; i < it->second.size( ); ++i )
        {
            if ( it->first == observation_models::dsn_n_way_averaged_doppler )
            {
                observationModelSettingsList.push_back(
                        observation_models::dsnNWayAveragedDopplerObservationSettings(
                                it->second.at( i ), lightTimeCorrectionSettings, constantTimeBias( 0.0, receiver ),
                                std::make_shared< LightTimeConvergenceCriteria >( true )  ) );
            }
        }
    }

    std::vector< std::shared_ptr< ObservationSimulatorBase< long double, Time > > > observationSimulators =
            createObservationSimulators< long double, Time >( observationModelSettingsList, bodies );

    /****************************************************************************************
     ************************** SIMULATE OBSERVATIONS AND COMPUTE RESIDUALS
     *****************************************************************************************/

    computeResidualsAndDependentVariables< long double, Time >( observedObservationCollection, observationSimulators, bodies );


    /****************************************************************************************
    ************************** FILTER OBSERVATIONS
    *****************************************************************************************/

    std::map< ObservableType, double > residualCutoffValuePerObservable;
    residualCutoffValuePerObservable[ dsn_n_way_averaged_doppler ] = 0.010;

    std::map< std::shared_ptr< ObservationCollectionParser >, std::shared_ptr< ObservationFilterBase > > filters;
    filters[ observationParser( dsn_n_way_averaged_doppler ) ] = observationFilter( residual_filtering, 0.010 );

    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredObservedObservationCollection =
            filterObservations( observedObservationCollection, filters );

    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > filteredObservationSimulationSettings =
            getObservationSimulationSettingsFromObservations( observedObservationCollection, bodies );

    /****************************************************************************************
    ************************** ARC STATISTICS
    *****************************************************************************************/

    Eigen::VectorXd startTimes;
    Eigen::VectorXd durations;
    Eigen::VectorXd meanValues;
    Eigen::VectorXd rmsValues;

    Eigen::VectorXd residualVector = observedObservationCollection->getConcatenatedResiduals( ).template cast< double >( );
    getResidualStatistics(
            observedObservationCollection,
            startTimes,
            durations,
            meanValues,
            rmsValues );

    Eigen::VectorXd numericalTimeBiasPartials = getNumericalObservationTimePartial< long double, Time >(
            filteredObservationSimulationSettings, observationSimulators, bodies, 5.0 );

    Eigen::VectorXd correctedResiduals;
    std::vector< double > timeBiases;
    estimateTimeBiasPerSet( observedObservationCollection, numericalTimeBiasPartials, timeBiases, correctedResiduals );

    input_output::writeMatrixToFile( startTimes, "grailTestStartTimes.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_StatisticsTest/");
    input_output::writeMatrixToFile( durations, "grailTestDurations.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_StatisticsTest/");
    input_output::writeMatrixToFile( meanValues, "grailMeanValues.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_StatisticsTest/");
    input_output::writeMatrixToFile( rmsValues, "grailRmsValues.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_StatisticsTest/");

    Eigen::VectorXd observationTimes = utilities::convertStlVectorToEigenVector(
            observedObservationCollection->getConcatenatedTimeVector( ) ).template cast< double >( );
    input_output::writeMatrixToFile( observationTimes, "grailTestTimes.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_StatisticsTest/");

    Eigen::VectorXd observationLinkEndsIds = utilities::convertStlVectorToEigenVector(
            observedObservationCollection->getConcatenatedLinkEndIds( ) ).template cast< double >( );
    input_output::writeMatrixToFile(observationLinkEndsIds , "grailTestLinkEnds.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_StatisticsTest/");

    input_output::writeMatrixToFile( residualVector, "grailUncorrectedResiduals.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_StatisticsTest/" );
    input_output::writeMatrixToFile( correctedResiduals, "grailTestCorrectedResiduals.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_StatisticsTest/");

    std::cout<<"Uncorrected mean: "<<linear_algebra::getVectorEntryMean( residualVector )<<std::endl;
    std::cout<<"Corrected mean: "<<linear_algebra::getVectorEntryMean( correctedResiduals )<<std::endl;

    std::cout<<"Uncorrected residual: "<<linear_algebra::getVectorEntryRootMeanSquare( residualVector )<<std::endl;
    std::cout<<"Corrected residual: "<<linear_algebra::getVectorEntryRootMeanSquare( correctedResiduals )<<std::endl;

//    {
//        std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface =
//                std::make_shared< propagators::SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >( nullptr, nullptr, 0, 1, std::vector< std::pair< int, int > >( ) );
//
//        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//        for ( auto it = linkEndsPerObservable.begin(); it != linkEndsPerObservable.end(); ++it )
//        {
//            for ( unsigned int i = 0; i < it->second.size( ); ++i )
//            {
//                parameterNames.push_back( timeObservationBias( it->second.at( i ), it->first ) );
//            }
//        }
//        parameterNames.push_back( groundStationPosition( "Earth", "DSS-24") );
//        std::shared_ptr< TranslationalStatePropagatorSettings< long double, Time > > propagatorSettings = nullptr;
//        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
//                createParametersToEstimate< long double, Time >( parameterNames, bodies, propagatorSettings );
//
//
//        std::shared_ptr< observation_models::ObservationManagerBase< long double, Time > > observationManagers = createObservationManagerBase< long double, Time >(
//                dsn_n_way_averaged_doppler,
//                observationModelSettingsList, bodies, parametersToEstimate,
//                stateTransitionMatrixInterface );
//    }


}