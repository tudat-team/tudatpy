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
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"

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

void runEstimation(
        std::string saveDirectory,
        std::string fileTag,
        Time initialEphemerisTime,
        Time finalEphemerisTime,
        bool useInterpolatedEphemerides,
        double epehemeridesTimeStep,
        int sphericalHarmonicsOrder,
        double integrationTolerance,
        std::pair< double, double > integrationMinMaxStep,
        int estimationMaxIterations,
        double observationsSamplingTime )
{

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = {
            "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Phobos", "Deimos",
            "Io", "Ganymede", "Callisto", "Europa", "Titan" };

    std::string baseFrameOrientation = "J2000";
    std::string baseFrameOrigin = "SSB";

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
    std::vector< aerodynamics::AtmosphereDependentVariables > atmosphereDependentVariables = {
        aerodynamics::specific_heat_ratio_dependent_atmosphere, aerodynamics::temperature_dependent_atmosphere,
        aerodynamics::density_dependent_atmosphere, aerodynamics::pressure_dependent_atmosphere,
        aerodynamics::gas_constant_dependent_atmosphere };
    std::vector< aerodynamics::AtmosphereIndependentVariables > atmosphereIndependentVariables = {
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
            tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables );

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

    Time initialPropagationTime = initialEphemerisTime + 600.0;
    Time finalPropagationTime = finalEphemerisTime - 600.0;
    // Select observation times
    // Compute observed observations. NOTE: don't move this to after the creation of the OrbitDeterminationManager!
    std::vector< Time > observationTimes;
    std::vector< Eigen::Matrix< long double, Eigen::Dynamic, 1 > > observations;
    for ( Time t = initialPropagationTime; t < finalPropagationTime; t += observationsSamplingTime )
    {
        try
        {
            observations.push_back( bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ).segment( 0, 3 ) );
            observationTimes.push_back( t );
        }
        catch( ... )
        { }
    }

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
    Eigen::Matrix< long double, 6, 1 > spacecraftInitialState =
            bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime ) -
            bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime );

    // Retrieve state history from SPICE
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > spiceStateHistory;
    for ( Time t : observationTimes )
    {
        spiceStateHistory[ t.getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
                    bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    }

    // Create termination settings
    std::shared_ptr< PropagationTerminationSettings > terminationSettings = propagationTimeTerminationSettings(
            finalPropagationTime );

    // Select dependent variables
//    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables = {
//            keplerianStateDependentVariable( spacecraftName, "Mars" ) };

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

    // Create link ends
    LinkEnds linkEnds;
    linkEnds[ observed_body ] = spacecraftName;

    // Create observation model settings
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
    observationModelSettingsList.push_back( positionObservableSettings( linkEnds ) );

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
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory;
    for ( Time t : observationTimes )
    {
        propagatedStateHistory[ t.getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
                    bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    }

    writeDataMapToTextFile( propagatedStateHistory, "stateHistoryPropagatedPreFit_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
    writeDataMapToTextFile( spiceStateHistory, "stateHistorySpice_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );

    std::vector< std::shared_ptr< SingleObservationSet< long double, Time > > > observationSetList;
    observationSetList.push_back(
            std::make_shared< SingleObservationSet< long double, Time > >(
                    position_observable, linkEnds, observations, observationTimes, observed_body ) );
    std::shared_ptr< ObservationCollection< long double, Time > > observedObservationCollection =
            std::make_shared< ObservationCollection< long double, Time > >( observationSetList );

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
            estimationOutput->getBestIterationSimulationResults( );
    std::map < Time, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFitDynamic =
            std::dynamic_pointer_cast< SingleArcVariationalSimulationResults< long double, Time > >(
                    postFitSimulationResults )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution( );
    std::map < Time, Eigen::Matrix < long double, 6, 1 > > propagatedStateHistoryPostFit;
    for ( auto it = propagatedStateHistoryPostFitDynamic.begin( ); it != propagatedStateHistoryPostFitDynamic.end( ); ++it )
    {
        propagatedStateHistoryPostFit[ it->first ] = it->second;
    }
    std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 6, 1 > > > postFitStateInterpolator =
            propagators::createStateInterpolator< Time, long double >( propagatedStateHistoryPostFit );
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFitToWrite;
    for ( Time t : observationTimes )
    {
        propagatedStateHistoryPostFitToWrite[ t.getSeconds< long double >() ] = postFitStateInterpolator->interpolate( t );
    }
//    for ( auto it = propagatedStateHistoryPostFit.begin( ); it != propagatedStateHistoryPostFit.end( ); ++it )
//    {
//        propagatedStateHistoryPostFitToWrite[ it->first.getSeconds< long double >() ] = it->second;
//    }
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

    std::pair< double, double > integrationMinMaxStep = std::make_pair( 1e-16, 1e16 );

    int testCase = 0;

    if ( testCase == 0 )
    {
        spice_interface::loadStandardSpiceKernels( );
        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgs-spice/mgs_ext22_ipng_mgs95j.bsp" );

        std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/estimation_position/";
        std::string fileTag = "november2005";
        double ephemeridesTimeStep = 50.0;
        double observationsSamplingTime = 40.0;

        // Select ephemeris time range (based on available data in loaded SPICE ephemeris)
        Time initialEphemerisTime = Time( 185976000 - 1.0 * 86400.0 ); // 23 November 2005, 0h
        Time finalEphemerisTime = Time( 186580800 + 1.0 * 86400.0 ); // 30 November 2005, 0h

        runEstimation( saveDirectory, fileTag,
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       120, 1e-10, integrationMinMaxStep, 5,
                       observationsSamplingTime );
    }
    else if ( testCase == 1 )
    {
        spice_interface::loadStandardSpiceKernels( );
//        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgs-spice/mgs_ext18_ipng_mgs95j.bsp" );
//        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgs-spice/mgs_ext19_ipng_mgs95j.bsp" );
//        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgs-spice/mgs_ext20_ipng_mgs95j.bsp" );
//        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgs-spice/mgs_ext21_ipng_mgs95j.bsp" );
        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgs-spice/mgs_ext22_ipng_mgs95j.bsp" );

        std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/estimation_position/";
        std::string fileTag = "2005";
        double ephemeridesTimeStep = 50.0;
        double observationsSamplingTime = 8000.0;

        // Select ephemeris time range (based on available data in loaded SPICE ephemeris)
//        Time initialEphemerisTime = Time( 157809600 + 0.0 * 86400.0 ); // 1 January 2005, 0h
//        Time finalEphemerisTime = Time( 189345600 - 0.0 * 86400.0 ); // 1 January 2006, 0h
        Time initialEphemerisTime = Time( 184766400 + 1.0 * 86400.0 ); // 9 November 2005, 0h
        Time finalEphemerisTime = Time( 192024000 - 1.0 * 86400.0 ); // 1 February 2006, 0h

        runEstimation( saveDirectory, fileTag,
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       120, 1e-10, integrationMinMaxStep, 5,
                       observationsSamplingTime );
    }


}

BOOST_AUTO_TEST_SUITE_END( )

}

}