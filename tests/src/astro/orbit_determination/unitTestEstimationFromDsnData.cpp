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
        std::string saveDirectory,
        std::string fileTag,
        bool useInterpolatedEphemerides,
        double epehemeridesTimeStep,
        std::vector< LightTimeCorrectionType > lightTimeCorrectionTypes )
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

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter" };

    std::string baseFrameOrientation = "J2000";
    std::string baseFrameOrigin = "SSB";

    // Define light-time perturbing bodies
    std::vector< std::string > lightTimePerturbingBodies = bodiesToCreate;

    // Specify initial time
    Time initialEphemerisTime = Time( 185976000 - 1.0 * 86400.0 ); // End of November 2005
    Time finalEphemerisTime = Time( 186580800 + 1.0 * 86400.0 ); // End of November 2005
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

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 10, 10 ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Venus" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Mercury" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
//    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );
//    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::aerodynamic ) );
    accelerationMap[ spacecraftName ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::string centralBody = "Mars";
    std::vector< std::string > centralBodies = { centralBody };
    bodiesToIntegrate.push_back( spacecraftName );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Create integrator settings
    double tolerance = 1e-11;
    double initialStep = 5.0;
    std::shared_ptr< IntegratorSettings< Time > > integratorSettings = rungeKuttaVariableStepSettingsScalarTolerances< Time >(
            initialStep, rungeKutta87DormandPrince, 1e-16, 1e16, tolerance, tolerance );

    // Set initial state from ephemerides
    Time initialPropagationTime = initialEphemerisTime;
    Eigen::Matrix< long double, 6, 1 > spacecraftInitialState =
            bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime ) -
            bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime );
    std::cout << std::setprecision(20) << "Initial state Spacecraft: " <<
        bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime ).transpose( ) << std::endl;
    std::cout << "Initial state Mars: " <<
        bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime ).transpose( ) << std::endl;

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
    std::map< Time, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > spiceStateHistory;
    for ( auto it = sampledTimes.begin( ); it != sampledTimes.end( ); ++it )
    {
        spiceStateHistory[ *it ] =
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
    std::shared_ptr< OdfRawFileContents > rawOdfFileContents =
            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5332333a.odf" );
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5327332a.odf" );

    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector = { rawOdfFileContents };
    std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents >(
                    rawOdfDataVector, bodies.getBody( "Earth" ), true, spacecraftName );

    // Create observed observation collection
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
            observation_models::createOdfObservedObservationCollection< long double, Time >(
                    processedOdfFileContents, { dsn_n_way_averaged_doppler } );

    // Set transmitting frequencies
    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );

    // Create computed observation collection
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;

    std::shared_ptr< input_output::CspRawFile > troposphericCspFile = std::make_shared< input_output::CspRawFile >(
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/tro/5305337a.tro"
            );
    std::shared_ptr< input_output::CspRawFile > ionosphericCspFile1 = std::make_shared< input_output::CspRawFile >(
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/ion/5305335g.ion"
            );
    std::shared_ptr< input_output::CspRawFile > ionosphericCspFile2 = std::make_shared< input_output::CspRawFile >(
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/ion/5335001a.ion"
            );

    std::map< int, std::string > spacecraftNamePerSpacecraftId;
    spacecraftNamePerSpacecraftId[ 94 ] = "MGS";

    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;

    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), first_order_relativistic ) )
    {
        lightTimeCorrectionSettings.push_back(
                std::make_shared< observation_models::FirstOrderRelativisticLightTimeCorrectionSettings >(
                    lightTimePerturbingBodies ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), tabulated_tropospheric ) )
    {
        lightTimeCorrectionSettings.push_back(
                std::make_shared< observation_models::TabulatedTroposphericCorrectionSettings >(
                        createTroposphericDryCorrectionAdjustment( { troposphericCspFile } ),
                        createTroposphericDryCorrectionAdjustment( { troposphericCspFile } ) ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), tabulated_ionospheric ) )
    {
        lightTimeCorrectionSettings.push_back(
                std::make_shared< observation_models::TabulatedIonosphericCorrectionSettings >(
                      createIonosphericCorrection( { ionosphericCspFile1, ionosphericCspFile2 },
                                                   spacecraftNamePerSpacecraftId ) ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), saastamoinen_tropospheric ) )
    {
        input_output::setDsnWeatherDataInGroundStations(
            bodies,
            std::vector< std::string >
                    { "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013321.wea",
                      "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013324.wea",
                      "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013326.wea" } );

        lightTimeCorrectionSettings.push_back( std::make_shared< SaastamoinenTroposphericCorrectionSettings >( ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), jakowski_vtec_ionospheric ) )
    {
        lightTimeCorrectionSettings.push_back( std::make_shared< JakowskiIonosphericCorrectionSettings >( ) );
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
//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, radiation_pressure_coefficient ) );
//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, constant_drag_coefficient ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
            createParametersToEstimate< long double, Time >( parameterNames, bodies );

    // Create orbit determination object.
    OrbitDeterminationManager< long double, Time > orbitDeterminationManager =
            OrbitDeterminationManager< long double, Time >(
                bodies, parametersToEstimate,
                observationModelSettingsList, propagatorSettings, true );

    // Retrieve state history
    std::shared_ptr< DynamicsSimulator< long double, Time > > dynamicsSimulator =
            orbitDeterminationManager.getVariationalEquationsSolver( )->getDynamicsSimulatorBase( );
//    std::map < Time, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory =
//            std::dynamic_pointer_cast< SingleArcSimulationResults< long double, Time > >(
//                    dynamicsSimulator->getPropagationResults( ) )->getEquationsOfMotionNumericalSolution( );

    std::map< Time, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory;
    for ( auto it = sampledTimes.begin( ); it != sampledTimes.end( ); ++it )
    {
        propagatedStateHistory[ *it ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( *it ) -
                    bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( *it );
    }

    writeDataMapToTextFile( propagatedStateHistory, "stateHistoryPropagated_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
    writeDataMapToTextFile( spiceStateHistory, "stateHistorySpice_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );

    // Define estimation input
    std::shared_ptr< EstimationInput< long double, Time  > > estimationInput =
            std::make_shared< EstimationInput< long double, Time > >( observedObservationCollection );

    // Perform estimation
    std::shared_ptr< EstimationOutput< long double, Time > > estimationOutput = orbitDeterminationManager.estimateParameters(
                estimationInput );

//    std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< long double, Time > > >
//            observationSimulators = observation_models::createObservationSimulators< long double, Time >(
//                    observationModelSettingsList, bodies );
//
//    std::vector< std::shared_ptr< ObservationSimulationSettings< Time > > > observationSimulationSettings =
//            observation_models::createOdfObservationSimulationSettingsList< long double, Time >(
//                    observedObservationCollection );
//
//    // Create observation collection
//    std::shared_ptr< observation_models::ObservationCollection< long double, Time > >
//            simulatedObservationCollection = simulation_setup::simulateObservations< long double, Time >(
//                    observationSimulationSettings, observationSimulators, bodies );
//
//    Eigen::Matrix< long double, Eigen::Dynamic, 1 > simulatedObservations =
//            simulatedObservationCollection->getObservationVector( );

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
    file << std::setprecision( 17 ) << residualHistory;
    file.close();

    std::ofstream file3(saveDirectory + "parameters_" + fileTag + ".txt");
    file3 << std::setprecision( 21 ) << parameterHistory;
    file3.close();

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
    std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/estimation/";

    int testCase = 0;

    if ( testCase == 0 )
    {
//        std::string fileTag = "5332333aOdf_interpState50";
        std::string fileTag = "5332333aOdf_interpState50";
        double ephemeridesTimeStep = 50.0;

        runSimulation( saveDirectory, fileTag + "_allCorr", true, ephemeridesTimeStep, { first_order_relativistic, tabulated_tropospheric, tabulated_ionospheric } );
    }


}

BOOST_AUTO_TEST_SUITE_END( )

}

}