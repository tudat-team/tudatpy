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

using namespace tudat;

int main( )
{
    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/MRO_Spice/mro_psp8_ssd_mro110c.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/MRO_Spice/mro_psp9_ssd_mro110c.bsp" );

    // Create settings for default bodies
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Phobos", "Deimos" };
    std::string globalFrameOrigin = "SSB";
    std::string globalFrameOrientation = "J2000";
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, globalFrameOrigin, globalFrameOrientation );

    // Add high-accuracy Earth settings
    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
        basic_astrodynamics::iau_2006, globalFrameOrientation );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create vector of atmosphere dependent and independent variables
    std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables = {
        aerodynamics::specific_heat_ratio_dependent_atmosphere, aerodynamics::temperature_dependent_atmosphere,
        aerodynamics::density_dependent_atmosphere, aerodynamics::pressure_dependent_atmosphere,
        aerodynamics::gas_constant_dependent_atmosphere };
    std::vector< aerodynamics::AtmosphereIndependentVariables > independentVariables = {
        aerodynamics::longitude_dependent_atmosphere, aerodynamics::latitude_dependent_atmosphere,
        aerodynamics::altitude_dependent_atmosphere };
    // Create a tabulated atmosphere object.
    std::map< int, std::string > tabulatedAtmosphereFiles;
    tabulatedAtmosphereFiles[ 0 ] = paths::getAtmosphereTablesPath( ) + "/MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";
    tabulatedAtmosphereFiles[ 1 ] = paths::getAtmosphereTablesPath( ) + "/MCDMeanAtmosphereTimeAverage/temperature.dat";
    tabulatedAtmosphereFiles[ 2 ] = paths::getAtmosphereTablesPath( ) + "/MCDMeanAtmosphereTimeAverage/density.dat";
    tabulatedAtmosphereFiles[ 3 ] = paths::getAtmosphereTablesPath( ) + "/MCDMeanAtmosphereTimeAverage/pressure.dat";
    tabulatedAtmosphereFiles[ 4 ] = paths::getAtmosphereTablesPath( ) + "/MCDMeanAtmosphereTimeAverage/gasConstant.dat";

    bodySettings.at( "Mars" )->atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >(
        tabulatedAtmosphereFiles, independentVariables, dependentVariables );

//    bodySettings.at( "Mars" )->rotationModelSettings = getHighAccuracyMarsRotationModel(
//        globalFrameOrientation, "IAU_Mars" );
    // Add spacecraft settings
    std::string spacecraftName = "MRO";
    std::string spacecraftCentralBody = "Mars";
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings =
        std::make_shared< InterpolatedSpiceEphemerisSettings >(
            77300 * 3600.0 , 77333.0 * 3600.0, 10.0, spacecraftCentralBody, globalFrameOrientation );
//        std::make_shared< DirectSpiceEphemerisSettings >( spacecraftCentralBody, "J2000" );

    // Create spacecraft
    bodySettings.at( spacecraftName )->constantMass = 1200.0;

    // Create radiation pressure settings
    double referenceAreaRadiation = 30.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies = { "Mars" };
    bodySettings.at( spacecraftName )->radiationPressureSettings[ "Sun" ] =
        cannonBallRadiationPressureSettings( "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 30.0;
    bodySettings.at( spacecraftName )->aerodynamicCoefficientSettings = std::make_shared< ConstantAerodynamicCoefficientSettings >(
        referenceArea, ( Eigen::Vector3d( ) << 2.0, 0.0, 0.0 ).finished( ) );


    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    // Define ODF data paths
    std::string dataDirectory = "/home/dominic/Tudat/Data/pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/odf/";
    std::vector< std::string > odfFiles = { "mromagr2008_300_1630xmmmv1.odf" };//, "mromagr2008_301_1615xmmmv1.odf", "mromagr2008_302_1605xmmmv1.odf"};

    // Laod raw ODF data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    for ( std::string odfFile : odfFiles )
    {
        rawOdfDataVector.push_back( std::make_shared<OdfRawFileContents>( dataDirectory + odfFile ) );
    }

    // Process ODF file data
    std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents =
        std::make_shared< ProcessedOdfFileContents >(
            rawOdfDataVector, spacecraftName, true );
    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );

    // Create data structure that handles Observed Data in Tudat
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
        observation_models::createOdfObservedObservationCollection< long double, Time >(
            processedOdfFileContents, { dsn_n_way_averaged_doppler } );

    std::pair< Time, Time > timeBounds = observedObservationCollection->getTimeBounds( );
    Time initialTime = timeBounds.first - 3600.0;
    Time finalTime = timeBounds.second + 3600.0;

    std::cout<<initialTime<<" "<<finalTime<<std::endl;

    // Define observation model settings
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;



    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( firstOrderRelativisticLightTimeCorrectionSettings( { "Sun" } ) );
//    std::vector< std::string > troposphericCorrectionFileNames =
//        {"/home/dominic/Tudat/Data/MRO_Ancilliary/mromagr2008_275_2008_306.tro"};
//    std::vector< std::string > ionosphericCorrectionFileNames =
//        {"/home/dominic/Tudat/Data/MRO_Ancilliary/mromagr2008_275_2008_306.ion"};
//    std::map< int, std::string > spacecraftNamePerSpacecraftId;
//    spacecraftNamePerSpacecraftId[ 74 ] = "MRO";
//    lightTimeCorrectionSettings.push_back( tabulatedTroposphericCorrectionSettings( troposphericCorrectionFileNames ) );
//    lightTimeCorrectionSettings.push_back( tabulatedIonosphericCorrectionSettings( ionosphericCorrectionFileNames, spacecraftNamePerSpacecraftId ) );

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
                        it->second.at( i ), lightTimeCorrectionSettings, nullptr,
                        std::make_shared< LightTimeConvergenceCriteria >( true ) ) );
            }
        }
    }

    std::vector< std::shared_ptr< ObservationSimulatorBase< long double, Time > > > observationSimulators =
        createObservationSimulators< long double, Time >( observationModelSettingsList, bodies );

    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > observationSimulationSettings =
        getObservationSimulationSettingsFromObservations( observedObservationCollection );
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > computedObservationCollection =
        simulateObservations( observationSimulationSettings, observationSimulators, bodies );

    std::cout<<"Create residuals"<<std::endl;
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > residualObservationCollection =
        createResidualCollection( observedObservationCollection, computedObservationCollection );

    std::cout<<"Filter observations"<<std::endl;
    std::map< ObservableType, double > residualCutoffValuePerObservable;
    residualCutoffValuePerObservable[ dsn_n_way_averaged_doppler ] = 0.3;
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredObservedObservationCollection =
        filterResidualOutliers( observedObservationCollection, residualObservationCollection, residualCutoffValuePerObservable );


    std::cout<<"Computed filtered observations"<<std::endl;
    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > filteredObservationSimulationSettings =
        getObservationSimulationSettingsFromObservations( filteredObservedObservationCollection );
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredComputedObservationCollection =
        simulateObservations( filteredObservationSimulationSettings, observationSimulators, bodies );

    std::cout<<"Create filtered residuals"<<std::endl;
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredResidualObservationCollection =
        createResidualCollection( filteredObservedObservationCollection, filteredComputedObservationCollection );



    {
        Eigen::VectorXd residuals = filteredResidualObservationCollection->getObservationVector( ).template cast< double >( );
        input_output::writeMatrixToFile( residuals, "mroTestResiduals.dat", 16, "/home/dominic/Tudat/Data/MRO_TestResults/");

        Eigen::VectorXd observationTimes = utilities::convertStlVectorToEigenVector(
            filteredResidualObservationCollection->getConcatenatedTimeVector( ) ).template cast< double >( );
        input_output::writeMatrixToFile( observationTimes, "mroTestTimes.dat", 16, "/home/dominic/Tudat/Data/MRO_TestResults/");

        Eigen::VectorXd observationLinkEndsIds = utilities::convertStlVectorToEigenVector(
            filteredResidualObservationCollection->getConcatenatedLinkEndIds( ) ).template cast< double >( );
        input_output::writeMatrixToFile(observationLinkEndsIds , "mroTestLinkEnds.dat", 16, "/home/dominic/Tudat/Data/MRO_TestResults/");
    }

    {
        // Set accelerations on Vehicle that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Sun" ].push_back( pointMassGravityAcceleration( ) );
        accelerationsOfVehicle[ "Sun" ].push_back( cannonBallRadiationPressureAcceleration( ) );
        accelerationsOfVehicle[ "Mercury" ].push_back( pointMassGravityAcceleration( ) );
        accelerationsOfVehicle[ "Venus" ].push_back( pointMassGravityAcceleration( ) );
        accelerationsOfVehicle[ "Earth" ].push_back( pointMassGravityAcceleration( ) );
        accelerationsOfVehicle[ "Mars" ].push_back( sphericalHarmonicAcceleration( 32, 32 ) );
        accelerationsOfVehicle[ "Mars" ].push_back( relativisticAccelerationCorrection(  ) );
        accelerationsOfVehicle[ "Mars" ].push_back( aerodynamicAcceleration( ) );
        accelerationsOfVehicle[ "Phobos" ].push_back( pointMassGravityAcceleration( ) );
        accelerationsOfVehicle[ "Deimos" ].push_back( pointMassGravityAcceleration( ) );
        accelerationsOfVehicle[ "Jupiter" ].push_back( pointMassGravityAcceleration( ) );
        accelerationsOfVehicle[ "Saturn" ].push_back( pointMassGravityAcceleration( ) );
        accelerationMap[ spacecraftName ] = accelerationsOfVehicle;

        // Create acceleration models
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, { spacecraftName }, { spacecraftCentralBody } );

        std::shared_ptr< IntegratorSettings< Time > > integratorSettings = rungeKuttaFixedStepSettings< Time >(
            Time( 30.0 ), rungeKutta87DormandPrince, RungeKuttaCoefficients::OrderEstimateToIntegrate::higher );

        Eigen::Matrix< long double, 6, 1 > spacecraftInitialState = spice_interface::getBodyCartesianStateAtEpoch(
            spacecraftName, spacecraftCentralBody, globalFrameOrientation, "None", initialTime ).template cast< long double >( );

        std::shared_ptr< PropagationTerminationSettings > terminationSettings = propagationTimeTerminationSettings(
            finalTime );

        std::shared_ptr< TranslationalStatePropagatorSettings< long double, Time > > propagatorSettings = translationalStatePropagatorSettings<
            long double, Time >( { spacecraftCentralBody },
                                 accelerationModelMap,
                                 { spacecraftName },
                                 spacecraftInitialState,
                                 initialTime,
                                 integratorSettings,
                                 terminationSettings,
                                 cowell );


        // Select parameters to estimate
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames = getInitialStateParameterSettings< long double, Time >(
            propagatorSettings, bodies );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, radiation_pressure_coefficient ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, constant_drag_coefficient ) );

        // Create parameters
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
            createParametersToEstimate< long double, Time >( parameterNames, bodies );

        std::cout<<"Starting propagation "<<std::endl;
        // Create orbit determination object.
        OrbitDeterminationManager< long double, Time > orbitDeterminationManager =
            OrbitDeterminationManager< long double, Time >(
                bodies, parametersToEstimate,
                observationModelSettingsList, propagatorSettings, true );

        // Define estimation input
        std::shared_ptr< EstimationInput< long double, Time  > > estimationInput =
            std::make_shared< EstimationInput< long double, Time > >(
                filteredObservedObservationCollection );
//        estimationInput->saveStateHistoryForEachIteration_ = true;

        // Perform estimation
        std::shared_ptr< EstimationOutput< long double, Time > > estimationOutput = orbitDeterminationManager.estimateParameters(
            estimationInput );


        {
            input_output::writeMatrixToFile( estimationOutput->residuals_, "mroTestConvergedResiduals.dat", 16, "/home/dominic/Tudat/Data/MRO_TestResults/");
        }



    }


}