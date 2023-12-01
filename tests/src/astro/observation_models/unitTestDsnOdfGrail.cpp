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
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/grail_v07.tf" );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/gra_sclkscet_00013.tsc" );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/gra_sclkscet_00014.tsc" );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/grail_120301_120529_sci_v02.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/gra_rec_120402_120408.bc" );
    spice_interface::loadSpiceKernelInTudat( get_spice_kernels_path( ) + + "/moon_pa_de440_200625.bpc");

    // Create settings for default bodies
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon" };
    std::string globalFrameOrigin = "Earth";
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
        "Antenna", ( Eigen::Vector3d( ) << -0.082, 0.152, -0.810 ).finished( ) );
    std::cout<<"Number of reference points: "<<bodies.at( "GRAIL-A" )->getVehicleSystems( )->getBodyFixedReferencePoints( ).size( )<<std::endl;

    /****************************************************************************************
     ************************** LOAD ODF FILES
     *****************************************************************************************/

    // Define ODF data paths
    std::string dataDirectory = "/home/dominic/Tudat/Data/GRAIL_ODF/";
    std::vector< std::string > odfFiles = { "gralugf2012_097_0235smmmv1.odf" };

    // Laod raw ODF data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    for ( std::string odfFile : odfFiles )
    {
        rawOdfDataVector.push_back( std::make_shared< OdfRawFileContents >( dataDirectory + odfFile ) );
    }

    // Process ODF file data
    std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents =
        std::make_shared< ProcessedOdfFileContents >(
            rawOdfDataVector, spacecraftName, true );
    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );

    // Create data structure that handles Observed Data in Tudat
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedUncompressedObservationCollection =
        splitObservationSetsIntoArcs< long double, Time >( observation_models::createOdfObservedObservationCollection< long double, Time >(
            processedOdfFileContents, { dsn_n_way_averaged_doppler } ), 60.0, 10 );

    std::map< LinkEnds, std::vector< std::shared_ptr< observation_models::SingleObservationSet< long double, Time > > > > uncompressedObservationSets =
        observedUncompressedObservationCollection->getObservations( ).at( dsn_n_way_averaged_doppler );
    std::vector< std::shared_ptr< observation_models::SingleObservationSet< long double, Time > > > compressedObservationSets;
    for( auto it : uncompressedObservationSets )
    {
        for( unsigned int i = 0; i < it.second.size( ); i++ )
        {
            compressedObservationSets.push_back( compressDopplerData< long double, Time >( it.second.at( i ), 60 ) );
        }
    }

    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
        std::make_shared< observation_models::ObservationCollection< long double, Time > >( compressedObservationSets );


    /****************************************************************************************
     ************************** PRINT DATA SUMMARY
     *****************************************************************************************/


    std::map< int, observation_models::LinkEnds > linkEndIds = observedObservationCollection->getInverseLinkEndIdentifierMap( );
    for( auto it : linkEndIds )
    {
        std::cout<<it.first<<", ("<<it.second[ transmitter ].bodyName_<<", "<<it.second[ transmitter ].stationName_<<"); "
                 <<", ("<<it.second[ retransmitter ].bodyName_<<", "<<it.second[ retransmitter ].stationName_<<"); "
                 <<", ("<<it.second[ receiver ].bodyName_<<", "<<it.second[ receiver ].stationName_<<")"<<std::endl;

    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< double, double > > > > arcStartEndTimes;
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > arcStartEndIndices;

    std::pair< Time, Time > timeBounds = observedObservationCollection->getTimeBounds( );
    Time initialTime = timeBounds.first - 3600.0;
    Time finalTime = timeBounds.second + 3600.0;

    std::cout<<"Initial time: "<<basic_astrodynamics::getCalendarDateFromTime( initialTime ).isoString( false, 3 )<<std::endl;
    std::cout<<"Final time: "<<basic_astrodynamics::getCalendarDateFromTime( finalTime ).isoString( false, 3 )<<std::endl;

    /****************************************************************************************
     ************************** CREATE OBSERVATION MODEL SETTINGS
     *****************************************************************************************/

    // Define observation model settings
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;

    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( firstOrderRelativisticLightTimeCorrectionSettings( { "Sun" } ) );
    std::vector< std::string > troposphericCorrectionFileNames =
        {"/home/dominic/Tudat/Data/GRAIL_Ancilliary/grxlugf2012_092_2012_122.tro", "/home/dominic/Tudat/Data/GRAIL_Ancilliary/grxlugf2012_122_2012_153.tro"};
    std::vector< std::string > ionosphericCorrectionFileNames =
        {"/home/dominic/Tudat/Data/GRAIL_Ancilliary/gralugf2012_092_2012_122.ion", "/home/dominic/Tudat/Data/GRAIL_Ancilliary/gralugf2012_122_2012_153.ion"};
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
                        it->second.at( i ), lightTimeCorrectionSettings, constantAbsoluteBias( Eigen::Vector1d::Zero( ) ),
                        std::make_shared< LightTimeConvergenceCriteria >( true )  ) );
            }
        }
    }

    std::vector< std::shared_ptr< ObservationSimulatorBase< long double, Time > > > observationSimulators =
        createObservationSimulators< long double, Time >( observationModelSettingsList, bodies );


    /****************************************************************************************
     ************************** SIMULATE OBSERVATIONS AND COMPUTE RESIDUALS
     *****************************************************************************************/

    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > observationSimulationSettings =
        getObservationSimulationSettingsFromObservations( observedObservationCollection );
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > computedObservationCollection =
        simulateObservations( observationSimulationSettings, observationSimulators, bodies );

    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > residualObservationCollection =
        createResidualCollection( observedObservationCollection, computedObservationCollection );

    /****************************************************************************************
    ************************** FILTER OBSERVATIONS
    *****************************************************************************************/

    std::map< ObservableType, double > residualCutoffValuePerObservable;
    residualCutoffValuePerObservable[ dsn_n_way_averaged_doppler ] = 0.010;

    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredObservedObservationCollection =
        filterResidualOutliers( observedObservationCollection, residualObservationCollection, residualCutoffValuePerObservable );


    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > filteredObservationSimulationSettings =
        getObservationSimulationSettingsFromObservations( filteredObservedObservationCollection );
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredComputedObservationCollection =
        simulateObservations( filteredObservationSimulationSettings, observationSimulators, bodies );

    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredResidualObservationCollection =
        createResidualCollection( filteredObservedObservationCollection, filteredComputedObservationCollection );
    {
        Eigen::VectorXd residuals = filteredResidualObservationCollection->getObservationVector( ).template cast< double >( );
        input_output::writeMatrixToFile( residuals, "grailTestResiduals.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults_30s_bias/");

//        input_output::writeMatrixToFile( correctedResiduals, "grailTestCorrectedResiduals.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");

        Eigen::VectorXd observationTimes = utilities::convertStlVectorToEigenVector(
            filteredResidualObservationCollection->getConcatenatedTimeVector( ) ).template cast< double >( );
        input_output::writeMatrixToFile( observationTimes, "grailTestTimes.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults_30s_bias/");

        Eigen::VectorXd observationLinkEndsIds = utilities::convertStlVectorToEigenVector(
            filteredResidualObservationCollection->getConcatenatedLinkEndIds( ) ).template cast< double >( );
        input_output::writeMatrixToFile(observationLinkEndsIds , "grailTestLinkEnds.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults_30s_bias/");
    }

    {
        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfSpacecraft;
        accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
        accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared<AccelerationSettings>( radiation_pressure ));
        accelerationsOfSpacecraft[ "Earth" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
        accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared<SphericalHarmonicAccelerationSettings>( 256, 256 ));
//        accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared<AccelerationSettings>( radiation_pressure ));
        accelerationsOfSpacecraft[ "Moon" ].push_back( empiricalAcceleration( ));
        accelerationsOfSpacecraft[ "Mars" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
        accelerationsOfSpacecraft[ "Jupiter" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
        accelerationsOfSpacecraft[ "Saturn" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));

        accelerationMap[ "GRAIL-A" ] =accelerationsOfSpacecraft;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector<std::string> bodiesToEstimate = { "GRAIL-A" };
        std::vector<std::string> centralBodies = { "Moon" };

        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToEstimate, centralBodies );

        Eigen::VectorXd initialState = getInitialStatesOfBodies( bodiesToEstimate, centralBodies, bodies, initialTime );

        std::shared_ptr< TranslationalStatePropagatorSettings< long double, Time > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< long double, Time > >
                ( centralBodies, accelerationModelMap, bodiesToEstimate, initialState.template cast< long double >( ), Time( initialTime ),
                  numerical_integrators::rungeKuttaFixedStepSettings< Time >( 30.0,
                                                                              numerical_integrators::rungeKuttaFehlberg78 ),
                  std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );


        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< long double, Time >( propagatorSettings, bodies );

        std::vector<std::shared_ptr<estimatable_parameters::EstimatableParameterSettings> > additionalParameterNames;

        parameterNames.push_back( estimatable_parameters::radiationPressureCoefficient( "GRAIL-A" ));

        std::map<basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector<basic_astrodynamics::EmpiricalAccelerationFunctionalShapes> > empiricalComponentsToEstimate;
        empiricalComponentsToEstimate[ radial_empirical_acceleration_component ].push_back( constant_empirical );
        empiricalComponentsToEstimate[ radial_empirical_acceleration_component ].push_back( sine_empirical );
        empiricalComponentsToEstimate[ radial_empirical_acceleration_component ].push_back( cosine_empirical );

        empiricalComponentsToEstimate[ along_track_empirical_acceleration_component ].push_back( constant_empirical );
        empiricalComponentsToEstimate[ along_track_empirical_acceleration_component ].push_back( sine_empirical );
        empiricalComponentsToEstimate[ along_track_empirical_acceleration_component ].push_back( cosine_empirical );

        empiricalComponentsToEstimate[ across_track_empirical_acceleration_component ].push_back( constant_empirical );
        empiricalComponentsToEstimate[ across_track_empirical_acceleration_component ].push_back( sine_empirical );
        empiricalComponentsToEstimate[ across_track_empirical_acceleration_component ].push_back( cosine_empirical );

        parameterNames.push_back( std::make_shared<EmpiricalAccelerationEstimatableParameterSettings>(
            "GRAIL-A", "Moon", empiricalComponentsToEstimate ));
//        parameterNames.push_back( std::make_shared<EstimatableParameterSettings>( "GRAIL-A", reference_point_position, "Antenna" ) );

        std::map < ObservableType, std::vector< LinkEnds > > linkEndsPerObservable =
            filteredComputedObservationCollection->getLinkEndsPerObservableType( );
        for( auto it : linkEndsPerObservable )
        {
            for( unsigned int i = 0; i < it.second.size( ); i++ )
            {
//                if( it.second.at( i ).at( transmitter ) != it.second.at( i ).at( receiver ) )
                {
                    parameterNames.push_back( observationBias( it.second.at( i ), it.first ) );
                }
            }
        }
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
            createParametersToEstimate< long double, Time >( parameterNames, bodies, propagatorSettings );

        OrbitDeterminationManager< long double, Time > orbitDeterminationManager = OrbitDeterminationManager< long double, Time >(
            bodies, parametersToEstimate, observationModelSettingsList, propagatorSettings );

        Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
        int numberOfParameters = truthParameters.rows( );

//        Eigen::MatrixXd inverseAprioriCovariance = Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters );
//
//        for( unsigned int i = 0; i < 3; i++ )
//        {
//            inverseAprioriCovariance( numberOfParameters - 3 + i, numberOfParameters - 3 + i ) = 1.0 / 4.0;
//        }
        std::shared_ptr< EstimationInput< long double, Time > > estimationInput = std::make_shared< EstimationInput< long double, Time > >(
            filteredObservedObservationCollection );//, inverseAprioriCovariance );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 4 ) );
        estimationInput->defineEstimationSettings(
            0, 1, 0, 1, 1, 1 );
        std::shared_ptr< EstimationOutput< long double, Time > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

        input_output::writeMatrixToFile(estimationOutput->residuals_ , "grailPostFitResiduals.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults_30s_bias/");

        auto estimatedStateHistory =
            std::dynamic_pointer_cast< SingleArcVariationalSimulationResults< long double, Time > >( estimationOutput->getSimulationResults( ).back( ) )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution( );

        std::map< double, Eigen::VectorXd > finalStateHistory;
        std::map< double, Eigen::VectorXd > finalStateDifference;
        std::map< double, Eigen::VectorXd > finalStateDifferenceRsw;
        for( auto it : estimatedStateHistory )
        {
            finalStateHistory[ it.first ] = it.second.cast< double >( );
            Eigen::Matrix3d rotationMatrix = getInertialToRswSatelliteCenteredFrameRotationMatrix( it.second.cast< double >( ) );
            Eigen::VectorXd stateDifference  = it.second.cast< double >( ) - spice_interface::getBodyCartesianStateAtEpoch(
                "GRAIL-A", "Moon", globalFrameOrientation, "None", it.first );
            Eigen::VectorXd rswStateDifference = Eigen::Vector6d::Zero( );
            rswStateDifference.segment( 0, 3 ) = rotationMatrix * stateDifference.segment( 0, 3 );
            rswStateDifference.segment( 3, 3 ) = rotationMatrix * stateDifference.segment( 3, 3 );

            finalStateDifference[ it.first ] = stateDifference;
            finalStateDifferenceRsw[ it.first ] = rswStateDifference;
        }

        input_output::writeMatrixToFile(estimationOutput->getCorrelationMatrix( ) , "grailTestCorrelations.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults_30s_bias/");

        input_output::writeDataMapToTextFile( finalStateDifference,
                                              "stateDifference.dat",
                                              "/home/dominic/Tudat/Data/GRAIL_TestResults_30s_bias/",
                                              "", std::numeric_limits< double >::digits10,  std::numeric_limits< double >::digits10,  "," );

        input_output::writeDataMapToTextFile( finalStateDifferenceRsw,
                                              "stateDifferenceRsw.dat",
                                              "/home/dominic/Tudat/Data/GRAIL_TestResults_30s_bias/",
                                              "", std::numeric_limits< double >::digits10,  std::numeric_limits< double >::digits10,  "," );





    }
}