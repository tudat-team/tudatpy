/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <string>

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulatePseudoObservations.h"
#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::ephemerides;
using namespace tudat::estimatable_parameters;
using namespace tudat::numerical_integrators;
using namespace tudat::observation_models;
using namespace tudat::propagators;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_parallel_multi_arc_propagation )

// Replace the rotation models of the bodies by tabulated rotation models, so that no SPICE calls are made during propagation
void setTabulatedRotationModels( const SystemOfBodies& bodies,
                                 const std::vector< std::string >& bodyNames,
                                 const double initialTime,
                                 const double finalTime )
{
    for( const std::string& bodyName : bodyNames )
    {
        bodies.at( bodyName )
                ->setRotationalEphemeris( getTabulatedRotationalEphemeris< double, double >(
                        bodies.at( bodyName )->getRotationalEphemeris( ), initialTime, finalTime, 3600.0 ) );
    }
}

// Propagate the multi-arc dynamics of unitTestMultiArcDynamics.cpp (case 0), with environment models that require no SPICE calls.
// For numberOfThreads equal to 1, the acceleration models are created on the bodies, otherwise from the acceleration settings for
// each thread. Returns the propagated states per arc, and the resulting Moon ephemeris (evaluated at given times). For multiple
// threads, it is checked whether the Moon ephemeris of the bodies of each thread is (only) reset if resetThreadBodies is true.
std::pair< std::vector< std::map< double, Eigen::VectorXd > >, std::vector< Eigen::Vector6d > > propagateEarthMoonArcs(
        const unsigned int numberOfThreads,
        const bool resetThreadBodies = false )
{
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 2.0E7;
    double buffer = 5.0 * 3600.0;

    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( bodySettings.at( "Moon" )->ephemerisSettings )
            ->resetFrameOrigin( "Earth" );
    bodySettings.at( "Moon" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );
    bodySettings.at( "Earth" )->atmosphereSettings = nullptr;
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    setTabulatedRotationModels( bodies, { "Earth", "Moon" }, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Moon" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    std::vector< std::string > bodiesToIntegrate = { "Moon" };
    std::vector< std::string > centralBodies = { "SSB" };

    std::vector< double > integrationArcStarts, integrationArcEnds;
    double currentStartTime = initialEphemerisTime + 1.0E4;
    double currentEndTime = currentStartTime + 1.0E6;
    do
    {
        integrationArcStarts.push_back( currentStartTime );
        integrationArcEnds.push_back( currentEndTime );
        currentStartTime = currentEndTime - 1.0E4;
        currentEndTime = currentStartTime + 1.0E6;
    } while( currentEndTime < finalEphemerisTime - 1.0E4 );

    AccelerationMap accelerationModelMap;
    if( numberOfThreads == 1 )
    {
        accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );
    }

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < integrationArcStarts.size( ); i++ )
    {
        Eigen::VectorXd arcInitialState = spice_interface::getBodyCartesianStateAtEpoch(
                bodiesToIntegrate.at( 0 ), "Earth", "ECLIPJ2000", "NONE", integrationArcStarts.at( i ) );
        std::shared_ptr< IntegratorSettings<> > integratorSettings = std::make_shared< IntegratorSettings<> >( rungeKutta4, 120.0 );
        std::shared_ptr< PropagationTerminationSettings > terminationSettings =
                std::make_shared< PropagationTimeTerminationSettings >( integrationArcEnds.at( i ) );
        if( numberOfThreads == 1 )
        {
            arcPropagationSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                        accelerationModelMap,
                                                                                        bodiesToIntegrate,
                                                                                        arcInitialState,
                                                                                        integrationArcStarts.at( i ),
                                                                                        integratorSettings,
                                                                                        terminationSettings ) );
        }
        else
        {
            arcPropagationSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                        accelerationMap,
                                                                                        bodiesToIntegrate,
                                                                                        arcInitialState,
                                                                                        integrationArcStarts.at( i ),
                                                                                        integratorSettings,
                                                                                        terminationSettings ) );
        }
    }
    std::shared_ptr< MultiArcPropagatorSettings< double > > multiArcPropagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double > >( arcPropagationSettingsList );
    multiArcPropagatorSettings->getOutputSettings( )->setIntegratedResult( true );
    multiArcPropagatorSettings->setNumberOfThreads( numberOfThreads );

    MultiArcDynamicsSimulator<> dynamicsSimulator( bodies, multiArcPropagatorSettings, false );
    dynamicsSimulator.setResetIntegratedStatesInThreadBodies( resetThreadBodies );
    dynamicsSimulator.integrateEquationsOfMotion( multiArcPropagatorSettings->getInitialStates( ) );

    std::vector< Eigen::Vector6d > moonEphemerisStates;
    for( unsigned int i = 0; i < integrationArcStarts.size( ); i++ )
    {
        double testTime = ( integrationArcStarts.at( i ) + integrationArcEnds.at( i ) ) / 2.0;
        moonEphemerisStates.push_back( bodies.at( "Moon" )->getEphemeris( )->getCartesianState( testTime ) );
        for( SystemOfBodies threadBodies : dynamicsSimulator.getThreadBodies( ) )
        {
            bool isThreadEphemerisEqual =
                    ( threadBodies.at( "Moon" )->getEphemeris( )->getCartesianState( testTime ) == moonEphemerisStates.back( ) );
            BOOST_CHECK_EQUAL( isThreadEphemerisEqual, resetThreadBodies );
        }
    }

    return std::make_pair( dynamicsSimulator.getMultiArcPropagationResults( )->getConcatenatedEquationsOfMotionResults( ),
                           moonEphemerisStates );
}

//! Test whether multi-threaded multi-arc propagation gives results identical to single-threaded propagation
BOOST_AUTO_TEST_CASE( testParallelMultiArcPropagation )
{
    spice_interface::loadStandardSpiceKernels( );

    std::pair< std::vector< std::map< double, Eigen::VectorXd > >, std::vector< Eigen::Vector6d > > singleThreadResults =
            propagateEarthMoonArcs( 1 );
    BOOST_CHECK( singleThreadResults.first.size( ) > 3 );

    // Use fewer threads than arcs (several arcs per thread), and (nearly) as many threads as arcs. For 3 threads, the integrated
    // states are also set in the bodies of each thread
    for( unsigned int numberOfThreads : { 2, 3, 10 } )
    {
        std::pair< std::vector< std::map< double, Eigen::VectorXd > >, std::vector< Eigen::Vector6d > > multiThreadResults =
                propagateEarthMoonArcs( numberOfThreads, numberOfThreads == 3 );

        BOOST_CHECK_EQUAL( multiThreadResults.first.size( ), singleThreadResults.first.size( ) );
        for( unsigned int i = 0; i < singleThreadResults.first.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( multiThreadResults.first.at( i ).size( ), singleThreadResults.first.at( i ).size( ) );
            for( auto stateIterator : singleThreadResults.first.at( i ) )
            {
                BOOST_CHECK( multiThreadResults.first.at( i ).at( stateIterator.first ) == stateIterator.second );
            }
            BOOST_CHECK( multiThreadResults.second.at( i ) == singleThreadResults.second.at( i ) );
        }
    }
}

// Estimate arc-wise atmosphere parameters as in unitTestEstimationArcWiseAtmosphereParameters.cpp, using four arcs, and with
// environment models that require no SPICE calls during propagation. For numberOfThreads equal to 1, the acceleration models are
// created on the bodies, otherwise from the acceleration settings for each thread. Returns the estimated parameters and residuals.
std::pair< Eigen::VectorXd, Eigen::VectorXd > estimateArcWiseAtmosphereParameters( const unsigned int numberOfThreads )
{
    // Cassini T022 and T068 states (c/a-5min) from unitTestEstimationArcWiseAtmosphereParameters.cpp; the two additional arcs use
    // the same Titan-centered states one day later
    std::vector< double > arcStartTimes = {
        2.205723865899772644e+08, 2.205723865899772644e+08 + 86400.0, 3.275976261739959717e+08, 3.275976261739959717e+08 + 86400.0
    };
    std::vector< Eigen::Vector6d > arcInitialStates;
    arcInitialStates.push_back( Eigen::Vector6d{ 2.166092477145108860e+06,
                                                 -7.925218812290439382e+05,
                                                 3.553645666049534921e+06,
                                                 3.304561684998871897e+03,
                                                 2.420278071381293103e+03,
                                                 -4.244915762044645817e+03 } );
    arcInitialStates.push_back( arcInitialStates.at( 0 ) );
    arcInitialStates.push_back( Eigen::Vector6d{ -3.091501373934821691e+06,
                                                 9.134845736499809427e+05,
                                                 -2.880028589145255741e+06,
                                                 4.617009140001175183e+03,
                                                 3.555047142774454642e+03,
                                                 -4.711669772164187862e+02 } );
    arcInitialStates.push_back( arcInitialStates.at( 2 ) );

    double initialEnvironmentTime = arcStartTimes.front( ) - 10.0 * 3600.0;
    double finalEnvironmentTime = arcStartTimes.back( ) + 10.0 * 3600.0;
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Titan", "Saturn" }, initialEnvironmentTime, finalEnvironmentTime, "Titan", "ECLIPJ2000", 3600.0 );
    bodySettings.addSettings( "Cassini" );
    bodySettings.at( "Cassini" )->ephemerisSettings =
            std::make_shared< TabulatedEphemerisSettings >( std::map< double, Eigen::Vector6d >( ), "Titan", "ECLIPJ2000" );
    bodySettings.at( "Cassini" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    setTabulatedRotationModels( bodies, { "Titan", "Saturn" }, initialEnvironmentTime, finalEnvironmentTime );

    // Modify the bodies after their creation
    bodies.at( "Cassini" )->setConstantBodyMass( 2.0E3 );
    bodies.at( "Cassini" )
            ->setAerodynamicCoefficientInterface( createAerodynamicCoefficientInterface(
                    std::make_shared< ConstantAerodynamicCoefficientSettings >(
                            20.0, Eigen::Vector3d( 2.5, -0.1, 0.05 ), aerodynamics::negative_aerodynamic_frame_coefficients ),
                    "Cassini",
                    bodies ) );
    bodies.at( "Titan" )->setAtmosphereModel(
            createAtmosphereModel( std::make_shared< ExponentialAtmosphereSettings >( 7.58e04, 175, 3.0553447e-04 ), "Titan" ) );

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Cassini" ][ "Titan" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[ "Cassini" ][ "Titan" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Cassini" ][ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    std::vector< std::string > bodiesToEstimate = { "Cassini" };
    std::vector< std::string > centralBodies = { "Titan" };

    AccelerationMap accelerationModelMap;
    if( numberOfThreads == 1 )
    {
        accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies );
    }

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList = { singleAccelerationDependentVariable(
            aerodynamic, "Cassini", "Titan" ) };

    std::vector< double > arcEndTimes;
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    for( unsigned int j = 0; j < arcStartTimes.size( ); j++ )
    {
        arcEndTimes.push_back( arcStartTimes.at( j ) + 10 * 60 );
        if( numberOfThreads == 1 )
        {
            propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies,
                    accelerationModelMap,
                    bodiesToEstimate,
                    Eigen::VectorXd( arcInitialStates.at( j ) ),
                    arcStartTimes.at( j ),
                    rungeKuttaFixedStepSettings( 1.0, rungeKuttaFehlberg78 ),
                    propagationTimeTerminationSettings( arcEndTimes.at( j ) ),
                    cowell,
                    dependentVariablesList ) );
        }
        else
        {
            propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies,
                    accelerationMap,
                    bodiesToEstimate,
                    Eigen::VectorXd( arcInitialStates.at( j ) ),
                    arcStartTimes.at( j ),
                    rungeKuttaFixedStepSettings( 1.0, rungeKuttaFehlberg78 ),
                    propagationTimeTerminationSettings( arcEndTimes.at( j ) ),
                    cowell,
                    dependentVariablesList ) );
        }
    }
    std::shared_ptr< MultiArcPropagatorSettings< double > > propagatorSettings = std::make_shared< MultiArcPropagatorSettings< double > >(
            propagatorSettingsList, false, std::make_shared< MultiArcPropagatorProcessingSettings >( false, true ) );
    propagatorSettings->setNumberOfThreads( numberOfThreads );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double >( propagatorSettings, bodies, arcStartTimes );
    parameterNames.push_back( arcwiseExponentialAtmosphereBaseDensity( "Titan", arcStartTimes ) );
    parameterNames.push_back( arcwiseExponentialAtmosphereScaleHeight( "Titan", arcStartTimes ) );

    // Propagate the true dynamics, and simulate observations of the resulting Cassini states
    MultiArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings, true );
    std::vector< double > observationTimes;
    for( unsigned int k = 0; k < arcStartTimes.size( ); k++ )
    {
        double currentTime = arcStartTimes.at( k ) + 5.0;
        while( currentTime < arcEndTimes.at( k ) - 5.0 )
        {
            observationTimes.push_back( currentTime );
            currentTime += 1.0;
        }
    }
    std::pair< std::vector< std::shared_ptr< ObservationModelSettings > >, std::shared_ptr< ObservationCollection< double > > >
            observationCollectionAndModelSettings = simulatePseudoObservations( bodies, bodiesToEstimate, centralBodies, observationTimes );

    // Create the estimated parameters (from their settings), perturb the atmosphere parameters, and estimate
    OrbitDeterminationManager< double > orbitDeterminationManager(
            bodies, parameterNames, observationCollectionAndModelSettings.first, propagatorSettings );
    Eigen::VectorXd perturbation( 8 );
    perturbation << -3e-05, 6e-05, -3e-05, 6e-05, 1e03, -1.5e03, 1e03, -1.5e03;
    Eigen::VectorXd parameterVector = orbitDeterminationManager.getParametersToEstimate( )->getFullParameterValues< double >( );
    parameterVector.tail( perturbation.size( ) ) += perturbation;
    orbitDeterminationManager.resetParameterEstimate( parameterVector );
    std::shared_ptr< EstimationInput< double > > estimationInput =
            std::make_shared< EstimationInput< double > >( observationCollectionAndModelSettings.second );
    estimationInput->defineEstimationSettings( 0, 1, 0, 1, 1, 1 );
    std::shared_ptr< EstimationOutput<> > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    return std::make_pair( estimationOutput->parameterEstimate_, estimationOutput->residuals_ );
}

//! Test whether estimation with multi-threaded multi-arc propagation gives results identical to single-threaded propagation
BOOST_AUTO_TEST_CASE( testParallelMultiArcEstimation )
{
    spice_interface::loadStandardSpiceKernels( );

    std::pair< Eigen::VectorXd, Eigen::VectorXd > singleThreadResults = estimateArcWiseAtmosphereParameters( 1 );
    BOOST_CHECK( !singleThreadResults.first.hasNaN( ) );

    // Use fewer threads than arcs (several arcs per thread), and as many threads as arcs
    for( unsigned int numberOfThreads : { 2, 3, 4 } )
    {
        std::pair< Eigen::VectorXd, Eigen::VectorXd > multiThreadResults = estimateArcWiseAtmosphereParameters( numberOfThreads );
        BOOST_CHECK( multiThreadResults.first == singleThreadResults.first );
        BOOST_CHECK( multiThreadResults.second == singleThreadResults.second );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
