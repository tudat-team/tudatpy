/*    Copyright (c) 2010-2023, Delft University of Technology
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
#include "tudat/simulation/estimation.h"
#include "tudat/simulation/estimation_setup.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::input_output;
using namespace tudat::simulation_setup;
using namespace tudat;

void runSimulation(
        std::vector< std::string > odfFiles,
        std::vector< std::string > troposphericCorrectionFileNames,
        std::vector< std::string > ionosphericCorrectionFileNames,
        std::vector< std::string > weatherFileNames,
        std::string saveDirectory,
        std::string fileTag,
        Time initialEphemerisTime,
        Time finalEphemerisTime,
        bool useInterpolatedEphemerides,
        double epehemeridesTimeStep,
        std::vector< LightTimeCorrectionType > lightTimeCorrectionTypes,
        int sphericalHarmonicsOrder,
        double integrationTolerance,
        std::pair< double, double > integrationMinMaxStep,
        int estimationMaxIterations )
{

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = {
            "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Phobos", "Deimos",
            "Io", "Ganymede", "Callisto", "Europa", "Titan" };

    std::string baseFrameOrientation = "J2000";
    std::string baseFrameOrigin = "SSB";

    // Define light-time perturbing bodies
    std::vector< std::string > lightTimePerturbingBodies = bodiesToCreate;

    // Specify ephemeris time steps and buffers
    Time ephemerisTimeStepPlanets = Time( epehemeridesTimeStep );
    Time bufferPlanets = Time( 10.0 * ephemerisTimeStepPlanets );
    Time ephemerisTimeStepSpacecraft = Time( epehemeridesTimeStep );
    Time bufferSpacecraft = Time( 10.0 * ephemerisTimeStepSpacecraft );

    // Create bodies settings needed in simulation
    BodyListSettings bodySettings;
    if ( useInterpolatedEphemerides )
    {
        bodySettings = getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - bufferPlanets, finalEphemerisTime + bufferPlanets,
                    baseFrameOrigin, baseFrameOrientation, ephemerisTimeStepPlanets );
    }
    else
    {
        bodySettings = getDefaultBodySettings( bodiesToCreate, baseFrameOrigin, baseFrameOrientation );
    }

    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
            basic_astrodynamics::iau_2006, baseFrameOrientation );
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

    // Create spacecraft
    std::string spacecraftName = "MGS";
    bodySettings.addSettings( spacecraftName );
    if ( useInterpolatedEphemerides )
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialEphemerisTime - bufferSpacecraft, finalEphemerisTime + bufferSpacecraft,
                        ephemerisTimeStepSpacecraft, baseFrameOrigin, baseFrameOrientation,
                        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ), spacecraftName );
    }
    else
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< DirectSpiceEphemerisSettings >( baseFrameOrigin, baseFrameOrientation );
    }
    bodySettings.at( spacecraftName )->constantMass = 700.0;

    // Create radiation pressure settings
    double referenceAreaRadiation = 20.0;
    double radiationPressureCoefficient = 1.0;
    std::vector< std::string > occultingBodies = { "Mars" };
    bodySettings.at( spacecraftName )->radiationPressureSettings[ "Sun" ] =
            cannonBallRadiationPressureSettings( "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 17.5;
    bodySettings.at( spacecraftName )->aerodynamicCoefficientSettings = std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, ( Eigen::Vector3d( ) << 2.0, 0.0, 0.0 ).finished( ),
                true, true );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Sun" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Sun" ].push_back( cannonBallRadiationPressureAcceleration( ) );
//    accelerationsOfVehicle[ "Sun" ].push_back( relativisticAccelerationCorrection(  ) );
    accelerationsOfVehicle[ "Mercury" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Venus" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Earth" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Mars" ].push_back( sphericalHarmonicAcceleration( sphericalHarmonicsOrder, sphericalHarmonicsOrder ) );
    accelerationsOfVehicle[ "Mars" ].push_back( relativisticAccelerationCorrection(  ) );
    accelerationsOfVehicle[ "Mars" ].push_back( aerodynamicAcceleration( ) );
    accelerationsOfVehicle[ "Phobos" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Deimos" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Jupiter" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Io" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Ganymede" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Callisto" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Europa" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Saturn" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Titan" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Uranus" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Neptune" ].push_back( pointMassGravityAcceleration( ) );

    accelerationMap[ spacecraftName ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::string centralBody = "Mars";
    std::vector< std::string > centralBodies = { centralBody };
    bodiesToIntegrate.push_back( spacecraftName );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Create integrator settings
    double initialStep = 5.0;
    if ( integrationMinMaxStep.first == integrationMinMaxStep.second )
    {
        initialStep = integrationMinMaxStep.first;
    }
    std::shared_ptr< IntegratorSettings< Time > > integratorSettings = rungeKuttaVariableStepSettingsScalarTolerances< Time >(
            initialStep, rungeKutta87DormandPrince, integrationMinMaxStep.first,
            integrationMinMaxStep.second, integrationTolerance, integrationTolerance );

    // Set initial state from ephemerides
    Time initialPropagationTime = initialEphemerisTime;
    Eigen::Matrix< long double, 6, 1 > spacecraftInitialState =
            bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime ) -
            bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime );

    // Retrieve state history from SPICE
    std::vector< Time > sampledTimes;
    sampledTimes.push_back( initialPropagationTime );
    while( 1 )
    {
        Time newTime = sampledTimes.back( ) + 60.0;
        if ( newTime > finalEphemerisTime )
        {
            break;
        }
        else
        {
            sampledTimes.push_back( newTime );
        }
    }
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > spiceStateHistory;
    for ( auto it = sampledTimes.begin( ); it != sampledTimes.end( ); ++it )
    {
        spiceStateHistory[ (*it).getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( *it ) -
                    bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( *it );
    }

    // Create termination settings
    std::shared_ptr< PropagationTerminationSettings > terminationSettings = propagationTimeTerminationSettings(
            finalEphemerisTime );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< long double, Time > > propagatorSettings = translationalStatePropagatorSettings<
            long double, Time >( centralBodies,
                                 accelerationModelMap,
                                 bodiesToIntegrate,
                                 spacecraftInitialState,
                                 initialPropagationTime,
                                 integratorSettings,
                                 terminationSettings,
                                 gauss_modified_equinoctial );

    // Read and process ODF file data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    for ( std::string odfFile : odfFiles )
        rawOdfDataVector.push_back( std::make_shared< OdfRawFileContents >( odfFile ) );

    std::shared_ptr< ProcessedOdfFileContents< Time > > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents< Time > >( rawOdfDataVector, spacecraftName );

    // Create observed observation collection
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
            observation_models::createOdfObservedObservationCollection< long double, Time >(
                    processedOdfFileContents, { dsn_n_way_averaged_doppler } );

    // Set transmitting frequencies
    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );

    // Create computed observation collection
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;

    std::map< int, std::string > spacecraftNamePerSpacecraftId;
    spacecraftNamePerSpacecraftId[ 94 ] = "MGS";

    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;

    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), first_order_relativistic ) )
    {
        lightTimeCorrectionSettings.push_back( firstOrderRelativisticLightTimeCorrectionSettings( lightTimePerturbingBodies ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), tabulated_tropospheric ) )
    {
        lightTimeCorrectionSettings.push_back( tabulatedTroposphericCorrectionSettings( troposphericCorrectionFileNames ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), tabulated_ionospheric ) )
    {
        lightTimeCorrectionSettings.push_back( tabulatedIonosphericCorrectionSettings( ionosphericCorrectionFileNames, spacecraftNamePerSpacecraftId ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), saastamoinen_tropospheric ) )
    {
        input_output::setDsnWeatherDataInGroundStations(bodies, weatherFileNames  );

        lightTimeCorrectionSettings.push_back( saastamoinenTroposphericCorrectionSettings( ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), jakowski_vtec_ionospheric ) )
    {
        lightTimeCorrectionSettings.push_back( jakowskiIonosphericCorrectionSettings( ) );
    }

    std::map < observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable =
            observedObservationCollection->getLinkEndsPerObservableType( );
    std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
            std::make_shared< LightTimeConvergenceCriteria >( true );
    for ( auto it = linkEndsPerObservable.begin(); it != linkEndsPerObservable.end(); ++it )
    {
        for ( unsigned int i = 0; i < it->second.size(); ++i )
        {
//            observationModelSettingsList.push_back(
//                    std::make_shared< observation_models::ObservationModelSettings >(
//                            it->first, it->second.at( i ), lightTimeCorrectionSettings, nullptr, nullptr ) );
            if ( it->first == observation_models::dsn_n_way_averaged_doppler )
            {
                observationModelSettingsList.push_back(
                    std::make_shared< observation_models::DsnNWayAveragedDopplerObservationSettings >(
                            it->second.at( i ), lightTimeCorrectionSettings, nullptr,
                            lightTimeConvergenceCriteria ) );
//                observationModelSettingsList.push_back(
//                    std::make_shared< observation_models::NWayDifferencedRangeObservationSettings >(
//                            it->second.at( i ), lightTimeCorrectionSettings, nullptr,
//                            lightTimeConvergenceCriteria ) );
//                observationModelSettingsList.push_back(
//                    std::make_shared< observation_models::NWayRangeObservationSettings >(
//                            it->second.at( i ),
//                            nullptr, 3, nullptr,
//                            lightTimeConvergenceCriteria ) );
            }
        }
    }

    // Select parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< long double > >(
            spacecraftName, spacecraftInitialState, centralBody, baseFrameOrientation ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, radiation_pressure_coefficient ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, constant_drag_coefficient ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
            createParametersToEstimate< long double, Time >( parameterNames, bodies );

    // Create orbit determination object.
    OrbitDeterminationManager< long double, Time > orbitDeterminationManager =
            OrbitDeterminationManager< long double, Time >(
                bodies, parametersToEstimate,
                observationModelSettingsList, propagatorSettings, true );

    // Retrieve state history
//    std::shared_ptr< DynamicsSimulator< long double, Time > > dynamicsSimulator =
//            orbitDeterminationManager.getVariationalEquationsSolver( )->getDynamicsSimulatorBase( );
//    std::map < Time, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory =
//            std::dynamic_pointer_cast< SingleArcSimulationResults< long double, Time > >(
//                    dynamicsSimulator->getPropagationResults( ) )->getEquationsOfMotionNumericalSolution( );

    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory;
    for ( auto it = sampledTimes.begin( ); it != sampledTimes.end( ); ++it )
    {
        propagatedStateHistory[ (*it).getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( *it ) -
                    bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( *it );
    }

    writeDataMapToTextFile( propagatedStateHistory, "stateHistoryPropagatedPreFit_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
    writeDataMapToTextFile( spiceStateHistory, "stateHistorySpice_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );

    // Define estimation input
    std::shared_ptr< EstimationInput< long double, Time  > > estimationInput =
            std::make_shared< EstimationInput< long double, Time > >(
                    observedObservationCollection,
                    Eigen::MatrixXd::Zero( 0, 0 ),
                    std::make_shared< EstimationConvergenceChecker >( estimationMaxIterations ) );
    estimationInput->saveStateHistoryForEachIteration_ = true;

    // Perform estimation
    std::shared_ptr< EstimationOutput< long double, Time > > estimationOutput = orbitDeterminationManager.estimateParameters(
                estimationInput );

    // Retrieve post-fit state history
    std::shared_ptr< propagators::SimulationResults< long double, Time > > postFitSimulationResults =
            estimationOutput->getSimulationResults( ).at( estimationOutput->bestIteration_ );
    std::map < Time, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFit =
            std::dynamic_pointer_cast< SingleArcVariationalSimulationResults< long double, Time > >(
                    postFitSimulationResults )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution( );
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFitToWrite;
    for ( auto it = propagatedStateHistoryPostFit.begin( ); it != propagatedStateHistoryPostFit.end( ); ++it )
    {
        propagatedStateHistoryPostFitToWrite[ it->first.getSeconds< long double >() ] = it->second;
    }
    writeDataMapToTextFile( propagatedStateHistoryPostFitToWrite, "stateHistoryPropagatedPostFit_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );

    // Retrieve residuals and set them in matrix
    Eigen::MatrixXd residualHistory = estimationOutput->getResidualHistoryMatrix( );
    Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > parameterHistory = estimationOutput->getParameterHistoryMatrix( );

    Eigen::MatrixXd residualsWithTime;
    residualsWithTime.resize( residualHistory.rows( ), residualHistory.cols( ) + 1 );
    residualsWithTime.rightCols( residualHistory.cols( ) ) = residualHistory;

    for ( unsigned int i = 0; i < observedObservationCollection->getObservationVector( ).size( ); ++i )
    {
        residualsWithTime( i, 0 ) = static_cast< Time >( observedObservationCollection->getConcatenatedTimeVector( ).at( i )
                ).getSeconds< long double >();
    }

    std::ofstream file(saveDirectory + "residuals_" + fileTag + ".txt");
    file << std::setprecision( 17 ) << residualsWithTime;
    file.close();

    std::ofstream file3(saveDirectory + "parameters_" + fileTag + ".txt");
    file3 << std::setprecision( 21 ) << parameterHistory;
    file3.close();

    // Retrieve covariance matrix
    Eigen::MatrixXd normalizedCovarianceMatrix = estimationOutput->getNormalizedCovarianceMatrix( );
    Eigen::MatrixXd unnormalizedCovarianceMatrix = estimationOutput->getUnnormalizedCovarianceMatrix( );
    std::ofstream file4(saveDirectory + "covariance_" + fileTag + ".txt");
    file4 << std::setprecision( 17 ) << "Normalized covariance matrix: " << std::endl << normalizedCovarianceMatrix;
    file4 << std::endl << std::endl;
    file4 << "Unnormalized covariance matrix: " << std::endl << unnormalizedCovarianceMatrix;
    file4.close( );

    std::ofstream file2(saveDirectory + "observationsStartAndSize_" + fileTag + ".txt");
    std::map< ObservableType, std::map< int, std::vector< std::pair< int, int > > > > observationSetStartAndSize =
            observedObservationCollection->getObservationSetStartAndSizePerLinkEndIndex();
    for ( auto it = observationSetStartAndSize.begin(); it != observationSetStartAndSize.end(); ++it )
    {
        ObservableType observable = it->first;
        for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2 )
        {
            int linkEnd = it2->first;
            for ( unsigned int i = 0; i < it2->second.size(); ++i )
            {
                file2 << std::setprecision( 15 ) << observable << " " << linkEnd << " " << it2->second.at( i ).first <<
                    " " << it2->second.at( i ).second << std::endl;
            }
        }
    }
    file2.close();

}

BOOST_AUTO_TEST_SUITE( test_estimation_from_dsn_data )

BOOST_AUTO_TEST_CASE( testDsnNWayAveragedDopplerModel )
{
        // MGS kernels
//    spice_interface::loadStandardSpiceKernels( {
//        "/Users/pipas/Documents/mro-spice/de414.bsp",
//        "/Users/pipas/Documents/mgs-spice/mar063.bsp",
//        "/Users/pipas/Documents/mgs-spice/mgs_ext5_ipng_mgs95j.bsp",
//        "/Users/pipas/Documents/mgs-spice/mgs_ext6_ipng_mgs95j.bsp" } );
//    spice_interface::loadStandardSpiceKernels( {
//        "/Users/pipas/Documents/mro-spice/de414.bsp",
//        "/Users/pipas/Documents/mgs-spice/mar063.bsp",
//        "/Users/pipas/Documents/mgs-spice/mgs_ext22_ipng_mgs95j.bsp" } );
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgs-spice/mgs_ext22_ipng_mgs95j.bsp" );

    std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/estimation/";

    std::vector< std::string > odfFiles = { "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5332333a.odf" };
//    std::vector< std::string > odfFiles = { "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5327332a.odf" };

    std::vector< std::string > tropCorrectionFiles = { "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/tro/5305337a.tro" };
    std::vector< std::string > ionCorrectionFiles = {
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/ion/5305335g.ion",
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/ion/5335001a.ion" };
    std::vector< std::string > weatherFiles = {
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013321.wea",
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013324.wea",
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013326.wea" };

    // Select ephemeris time range (based on available data in loaded SPICE ephemeris)
    Time initialEphemerisTime = Time( 185976000 - 1.0 * 86400.0 ); // End of November 2005
    Time finalEphemerisTime = Time( 186580800 + 1.0 * 86400.0 ); // End of November 2005

    std::pair< double, double > integrationMinMaxStep = std::make_pair( 1e-16, 1e16 );

    int testCase = 4;

    if ( testCase == 0 )
    {
//        std::string fileTag = "5332333aOdf_interpState50";
        std::string fileTag = "5332333aOdf_interpState50";
        double ephemeridesTimeStep = 50.0;

        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorr",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 10, 1e-11,
                       integrationMinMaxStep, 0 );
    }
    else if ( testCase == 1 )
    {
//        std::string fileTag = "5332333aOdf_interpState50";
        std::string fileTag = "5332333aOdf_interpState50";
        double ephemeridesTimeStep = 50.0;

        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh10",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 10, 1e-11,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh20",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 20, 1e-11,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh30",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 30, 1e-11,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh40", \
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 40, 1e-11,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-11,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh60",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 60, 1e-11,
                       integrationMinMaxStep, 0 );
    }
    else if ( testCase == 2 )
    {
//        std::string fileTag = "5332333aOdf_interpState50";
        std::string fileTag = "5332333aOdf_interpState50";
        double ephemeridesTimeStep = 50.0;

        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50Tol5",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-5,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50Tol6",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-6,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50Tol7",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-7,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50Tol8",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-8,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50Tol9",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-9,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50Tol10",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-10,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50Tol11",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-11,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50Tol12",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-12,
                       integrationMinMaxStep, 0 );
    }
    else if ( testCase == 3)
    {
        std::string fileTag = "5332333aOdf_interpState50";
        double ephemeridesTimeStep = 50.0;
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh50Tol10",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 50, 1e-10,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh40Tol10",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 40, 1e-10,
                       integrationMinMaxStep, 0 );
        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_allCorrSh60Tol10",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 60, 1e-10,
                       integrationMinMaxStep, 0 );
    }
    else if ( testCase == 4)
    {
        std::string fileTag = "5332333aOdf_interpState50";
        double ephemeridesTimeStep = 50.0;

//        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
//                       saveDirectory, fileTag + "_allCorrSh120Tol10Moons", true, ephemeridesTimeStep,
//                       { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric }, 120, 1e-10, integrationMinMaxStep, 5 );

        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_noCorrSh120Tol10Moons",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { }, 120, 1e-10, integrationMinMaxStep, 5 );

//        runSimulation( odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
//                       saveDirectory, fileTag + "_noCorrSh120Tol10MoonsStep20",
//                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
//                       { }, 120, 1e16, std::make_pair( 20.0, 20.0 ), 1 );
    }


}

BOOST_AUTO_TEST_SUITE_END( )

}

}