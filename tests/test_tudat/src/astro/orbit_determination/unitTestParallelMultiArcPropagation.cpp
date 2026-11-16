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

#include "tudat/basics/testMacros.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
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

// A small native-only acceptance case: ordinary settings, shared mu, unequal arc durations,
// and more arcs than worker environments. No SPICE, files, or Python callbacks are used.
struct ClosureArcRun {
    SystemOfBodies bodies;
    std::shared_ptr< MultiArcPropagatorSettings< double > > settings;
    std::shared_ptr< MultiArcVariationalEquationsSolver<> > solver;
};

ClosureArcRun makeClosureArcRun( const unsigned int workers, const bool customEnvironment = false )
{
    ClosureArcRun run;
    BodyListSettings bodySettings( "Central", "J2000" );
    bodySettings.addSettings( "Central" );
    bodySettings.at( "Central" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings.at( "Central" )->gravityFieldSettings = centralGravitySettings( 3.986004418E14 );
    bodySettings.addSettings( "Satellite" );
    bodySettings.at( "Satellite" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "Central", "J2000" );
    bodySettings.at( "Satellite" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    run.bodies = createSystemOfBodies( bodySettings );
    if( customEnvironment )
    {
        run.bodies.createEmptyBody( "Custom" );
        run.bodies.at( "Custom" )
                ->setEphemeris( std::make_shared< CustomEphemeris<> >(
                        []( double ) { return Eigen::Vector6d::Zero( ).eval( ); }, "SSB", "J2000" ) );
    }
    SelectedAccelerationMap accelerations;
    accelerations[ "Satellite" ][ "Central" ].push_back( pointMassGravityAcceleration( ) );
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcs;
    std::vector< double > starts;
    for( unsigned int arc = 0; arc < 4; ++arc )
    {
        const double start = 1000.0 * arc;
        starts.push_back( start );
        Eigen::Vector6d initialState;
        initialState << 7.0E6 + 100.0 * arc, 0.0, 0.0, 0.0, 7500.0, 10.0;
        auto settings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                std::vector< std::string >{ "Central" },
                AccelerationMap( ),
                std::vector< std::string >{ "Satellite" },
                initialState,
                start,
                rungeKuttaFixedStepSettings( 10.0, rungeKuttaFehlberg78 ),
                propagationTimeTerminationSettings( start + 100.0 + 20.0 * arc ),
                cowell,
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >{
                        relativePositionDependentVariable( "Satellite", "Central" ) } );
        settings->resetAccelerationModelsMap( accelerations, run.bodies );
        arcs.push_back( settings );
    }
    run.settings = std::make_shared< MultiArcPropagatorSettings< double > >( arcs, false );
    run.settings->setNumberOfThreads( workers );
    run.settings->getOutputSettings( )->setIntegratedResult( true );
    run.settings->getOutputSettings( )->setUpdateDependentVariableInterpolator( true );
    auto parameters = getInitialMultiArcParameterSettings< double >( run.settings, run.bodies, starts );
    parameters.push_back( std::make_shared< EstimatableParameterSettings >( "Central", gravitational_parameter ) );
    run.solver = std::make_shared< MultiArcVariationalEquationsSolver<> >( run.bodies, run.settings, parameters, false );
    const auto environments = run.solver->getDynamicsSimulator( )->getThreadBodies( );
    BOOST_CHECK_EQUAL( environments.size( ), workers > 1 ? workers : 0 );
    for( unsigned int worker = 0; worker < environments.size( ); ++worker )
    {
        const auto gravity = environments.at( worker ).at( "Central" )->getGravityFieldModel( );
        BOOST_CHECK( gravity != run.bodies.at( "Central" )->getGravityFieldModel( ) );
        for( unsigned int other = 0; other < worker; ++other )
        {
            BOOST_CHECK( gravity != environments.at( other ).at( "Central" )->getGravityFieldModel( ) );
        }
    }
    for( unsigned int arc = 0; arc < 4; ++arc )
    {
        BOOST_CHECK_EQUAL( run.solver->getDynamicsSimulator( )->getArcThreadIndex( arc ), workers > 1 ? arc % workers : 0 );
    }
    return run;
}

// Mixed absolute/relative comparison also checks entries whose reference is zero.
void checkClosureMatrix( const Eigen::MatrixXd& actual, const Eigen::MatrixXd& expected )
{
    BOOST_REQUIRE_EQUAL( actual.rows( ), expected.rows( ) );
    BOOST_REQUIRE_EQUAL( actual.cols( ), expected.cols( ) );
    BOOST_REQUIRE( actual.size( ) > 0 );
    BOOST_REQUIRE( actual.allFinite( ) && expected.allFinite( ) );
    const double normalizedError = ( ( actual - expected ).array( ).abs( ) / ( 1.0 + expected.array( ).abs( ) ) ).maxCoeff( );
    BOOST_CHECK_LE( normalizedError, 1.0E-11 );
}

void compareClosureRuns( ClosureArcRun& actual, ClosureArcRun& expected )
{
    BOOST_REQUIRE( actual.solver->getDynamicsSimulator( )->integrationCompletedSuccessfully( ) );
    BOOST_REQUIRE( expected.solver->getDynamicsSimulator( )->integrationCompletedSuccessfully( ) );
    const auto a = actual.solver->getDynamicsSimulator( )->getSingleArcDynamicsSimulators( );
    const auto b = expected.solver->getDynamicsSimulator( )->getSingleArcDynamicsSimulators( );
    const auto av = actual.solver->getNumericalVariationalEquationsSolution( );
    const auto bv = expected.solver->getNumericalVariationalEquationsSolution( );
    BOOST_REQUIRE_EQUAL( a.size( ), b.size( ) );
    for( unsigned int arc = 0; arc < a.size( ); ++arc )
    {
        const auto states = a.at( arc )->getEquationsOfMotionNumericalSolution( );
        BOOST_REQUIRE( !states.empty( ) );
        BOOST_REQUIRE_EQUAL( states.size( ), b.at( arc )->getEquationsOfMotionNumericalSolution( ).size( ) );
        for( const auto& entry : states )
        {
            checkClosureMatrix( entry.second, b.at( arc )->getEquationsOfMotionNumericalSolution( ).at( entry.first ) );
            checkClosureMatrix( a.at( arc )->getDependentVariableHistory( ).at( entry.first ),
                                b.at( arc )->getDependentVariableHistory( ).at( entry.first ) );
            for( unsigned int kind = 0; kind < 2; ++kind )
            {
                // Express the tiny d(state)/d(mu) column per fractional mu change so
                // an absolute tolerance cannot hide a missing or sign-flipped sensitivity.
                const double scale =
                        kind == 1 ? expected.bodies.at( "Central" )->getGravityFieldModel( )->getGravitationalParameter( ) : 1.0;
                checkClosureMatrix( scale * av.at( arc ).at( kind ).at( entry.first ), scale * bv.at( arc ).at( kind ).at( entry.first ) );
            }
        }
        BOOST_CHECK_EQUAL( a.at( arc )->getSingleArcPropagationResults( )->getTotalNumberOfFunctionEvaluations( ),
                           b.at( arc )->getSingleArcPropagationResults( )->getTotalNumberOfFunctionEvaluations( ) );
        BOOST_CHECK_EQUAL( a.at( arc )->getPropagationTerminationReason( )->getPropagationTerminationReason( ),
                           b.at( arc )->getPropagationTerminationReason( )->getPropagationTerminationReason( ) );
        const double epoch = states.rbegin( )->first;
        checkClosureMatrix( actual.bodies.at( "Satellite" )->getEphemeris( )->getCartesianState( epoch ),
                            expected.bodies.at( "Satellite" )->getEphemeris( )->getCartesianState( epoch ) );
        const double middle = 0.5 * ( states.begin( )->first + epoch );
        checkClosureMatrix( actual.solver->getDynamicsSimulator( )
                                    ->getMultiArcPropagationResults( )
                                    ->getMultiArcDependentVariablesInterface( )
                                    ->getDependentVariables( middle ),
                            expected.solver->getDynamicsSimulator( )
                                    ->getMultiArcPropagationResults( )
                                    ->getMultiArcDependentVariablesInterface( )
                                    ->getDependentVariables( middle ) );

        Eigen::MatrixXd matrix =
                actual.solver->getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( epoch );
        Eigen::MatrixXd referenceMatrix =
                expected.solver->getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( epoch );
        BOOST_REQUIRE_EQUAL( matrix.cols( ), 25 );  // 4 * 6 initial states, one shared mu column
        const double mu = expected.bodies.at( "Central" )->getGravityFieldModel( )->getGravitationalParameter( );
        matrix.col( 24 ) *= mu;
        referenceMatrix.col( 24 ) *= mu;
        checkClosureMatrix( matrix, referenceMatrix );
        checkClosureMatrix( actual.bodies.at( "Satellite" )->getEphemeris( )->getCartesianState( epoch ),
                            expected.bodies.at( "Satellite" )->getEphemeris( )->getCartesianState( epoch ) );
    }
}

BOOST_AUTO_TEST_CASE( testParallelSetupRejectsUnverifiedSharing )
{
    BOOST_CHECK_THROW( makeClosureArcRun( 2, true ), std::runtime_error );
    BOOST_CHECK_NO_THROW( makeClosureArcRun( 1, true ) );
}

BOOST_AUTO_TEST_CASE( testClosureArcParameterAndResultContract )
{
    auto reference = makeClosureArcRun( 1 );
    reference.solver->integrateVariationalAndDynamicalEquations( reference.settings->getInitialStateList( ), true );
    for( unsigned int workers : { 2, 3, 4 } )
    {
        auto parallel = makeClosureArcRun( workers );
        auto dynamics = parallel.solver->getDynamicsSimulator( );
        const auto results = dynamics->getMultiArcPropagationResults( );
        const auto variationalResults = parallel.solver->getVariationalPropagationResults( );
        parallel.solver->integrateVariationalAndDynamicalEquations( parallel.settings->getInitialStateList( ), true );
        const auto interface = parallel.solver->getStateTransitionMatrixInterface( );
        compareClosureRuns( parallel, reference );
        const auto environments = dynamics->getThreadBodies( );
        BOOST_REQUIRE_EQUAL( environments.size( ), workers );
        for( unsigned int arc = 0; arc < 4; ++arc )
        {
            BOOST_CHECK_EQUAL( dynamics->getArcThreadIndex( arc ), arc % workers );
        }
        for( unsigned int worker = 0; worker < workers; ++worker )
        {
            BOOST_CHECK( environments.at( worker ).at( "Central" )->getGravityFieldModel( ) !=
                         parallel.bodies.at( "Central" )->getGravityFieldModel( ) );
            for( unsigned int other = worker + 1; other < workers; ++other )
            {
                BOOST_CHECK( environments.at( worker ).at( "Central" )->getGravityFieldModel( ) !=
                             environments.at( other ).at( "Central" )->getGravityFieldModel( ) );
            }
        }
        // A fresh reference per worker count keeps this reset independent of previous iterations.
        auto changedReference = makeClosureArcRun( 1 );
        Eigen::VectorXd values = reference.solver->getParametersToEstimate( )->getFullParameterValues< double >( );
        values( 0 ) += 10.0;
        values( values.size( ) - 1 ) *= 1.00001;
        parallel.solver->resetParameterEstimate( values, true );
        changedReference.solver->resetParameterEstimate( values, true );
        compareClosureRuns( parallel, changedReference );
        const auto changedStates = dynamics->getSingleArcDynamicsSimulators( ).front( )->getEquationsOfMotionNumericalSolution( );
        BOOST_REQUIRE( !changedStates.empty( ) );
        BOOST_CHECK_SMALL( changedStates.begin( )->second( 0 ) - values( 0 ), 1.0E-8 );
        BOOST_CHECK( std::abs( changedStates.rbegin( )->second( 0 ) -
                               reference.solver->getDynamicsSimulator( )
                                       ->getSingleArcDynamicsSimulators( )
                                       .front( )
                                       ->getEquationsOfMotionNumericalSolution( )
                                       .rbegin( )
                                       ->second( 0 ) ) > 1.0 );

        BOOST_CHECK( dynamics->getMultiArcPropagationResults( ) == results );
        BOOST_CHECK( parallel.solver->getVariationalPropagationResults( ) == variationalResults );
        BOOST_CHECK( parallel.solver->getStateTransitionMatrixInterface( ) == interface );
        for( const auto& environment : environments )
        {
            BOOST_CHECK_EQUAL( environment.at( "Central" )->getGravityFieldModel( )->getGravitationalParameter( ),
                               values( values.size( ) - 1 ) );
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
    for( const std::string& bodyName : std::vector< std::string >{ "Titan", "Saturn" } )
    {
        // SPICE cannot be used during multi-threaded propagation
        bodySettings.at( bodyName )->rotationModelSettings =
                simpleRotationModelFromSpiceSettings( "ECLIPJ2000", "IAU_" + bodyName, "IAU_" + bodyName, initialEnvironmentTime );
    }
    bodySettings.at( "Titan" )->atmosphereSettings = std::make_shared< ExponentialAtmosphereSettings >( 7.58e04, 175, 3.0553447e-04 );
    bodySettings.addSettings( "Cassini" );
    bodySettings.at( "Cassini" )->ephemerisSettings =
            std::make_shared< TabulatedEphemerisSettings >( std::map< double, Eigen::Vector6d >( ), "Titan", "ECLIPJ2000" );
    bodySettings.at( "Cassini" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Cassini" )->constantMass = 2.0E3;
    bodySettings.at( "Cassini" )->aerodynamicCoefficientSettings = std::make_shared< ConstantAerodynamicCoefficientSettings >(
            20.0, Eigen::Vector3d( 2.5, -0.1, 0.05 ), aerodynamics::negative_aerodynamic_frame_coefficients );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Cassini" ][ "Titan" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[ "Cassini" ][ "Titan" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Cassini" ][ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    std::vector< std::string > bodiesToEstimate = { "Cassini" };
    std::vector< std::string > centralBodies = { "Titan" };

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList = { singleAccelerationDependentVariable(
            aerodynamic, "Cassini", "Titan" ) };

    std::vector< double > arcEndTimes;
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    for( unsigned int j = 0; j < arcStartTimes.size( ); j++ )
    {
        arcEndTimes.push_back( arcStartTimes.at( j ) + 10 * 60 );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies,
                        AccelerationMap( ),
                        bodiesToEstimate,
                        Eigen::VectorXd( arcInitialStates.at( j ) ),
                        arcStartTimes.at( j ),
                        rungeKuttaFixedStepSettings( 1.0, rungeKuttaFehlberg78 ),
                        propagationTimeTerminationSettings( arcEndTimes.at( j ) ),
                        cowell,
                        dependentVariablesList );
        translationalSettings->resetAccelerationModelsMap( accelerationMap, bodies );
        propagatorSettingsList.push_back( translationalSettings );
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

//! Test whether estimation with multi-threaded multi-arc propagation gives results equal (to numerical noise) to single-threaded propagation
BOOST_AUTO_TEST_CASE( testParallelMultiArcEstimation )
{
    spice_interface::loadStandardSpiceKernels( );

    std::pair< Eigen::VectorXd, Eigen::VectorXd > singleThreadResults = estimateArcWiseAtmosphereParameters( 1 );
    BOOST_CHECK( !singleThreadResults.first.hasNaN( ) );

    // Use fewer threads than arcs (several arcs per thread), and as many threads as arcs
    for( unsigned int numberOfThreads : { 2, 3, 4 } )
    {
        std::pair< Eigen::VectorXd, Eigen::VectorXd > multiThreadResults = estimateArcWiseAtmosphereParameters( numberOfThreads );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( multiThreadResults.first, singleThreadResults.first, 1.0E-8 );
        BOOST_CHECK_SMALL( ( multiThreadResults.second - singleThreadResults.second ).cwiseAbs( ).maxCoeff( ), 1.0E-6 );
    }
}

// Propagate Titan (around Saturn) and Cassini (around Titan) in two arcs, and return Cassini's ephemeris states at the middle of
// each arc. Cassini's ephemeris is w.r.t. Saturn, so setting the propagated states in its ephemeris requires Titan's state, which
// differs from Titan's ephemeris before propagation (initial state offset of 1 km)
std::vector< Eigen::Vector6d > propagateWithPropagatedCentralBody( const unsigned int numberOfThreads )
{
    std::vector< double > arcStartTimes = { 2.205723865899772644e+08, 2.205723865899772644e+08 + 86400.0 };
    double initialEnvironmentTime = arcStartTimes.front( ) - 86400.0;
    double finalEnvironmentTime = arcStartTimes.back( ) + 86400.0;
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Titan", "Saturn" }, initialEnvironmentTime, finalEnvironmentTime, "Saturn", "ECLIPJ2000", 600.0 );
    for( const std::string& bodyName : std::vector< std::string >{ "Titan", "Saturn" } )
    {
        bodySettings.at( bodyName )->rotationModelSettings =
                simpleRotationModelFromSpiceSettings( "ECLIPJ2000", "IAU_" + bodyName, "IAU_" + bodyName, initialEnvironmentTime );
    }
    bodySettings.at( "Titan" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.addSettings( "Cassini" );
    bodySettings.at( "Cassini" )->ephemerisSettings =
            std::make_shared< TabulatedEphemerisSettings >( std::map< double, Eigen::Vector6d >( ), "Saturn", "ECLIPJ2000" );
    bodySettings.at( "Cassini" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Titan" ][ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Cassini" ][ "Titan" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Cassini" ][ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    std::vector< std::string > bodiesToPropagate = { "Titan", "Cassini" };
    std::vector< std::string > centralBodies = { "Saturn", "Titan" };

    Eigen::Vector6d cassiniInitialState{ 2.166092477145108860e+06, -7.925218812290439382e+05, 3.553645666049534921e+06,
                                         3.304561684998871897e+03, 2.420278071381293103e+03,  -4.244915762044645817e+03 };
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    for( unsigned int j = 0; j < arcStartTimes.size( ); j++ )
    {
        Eigen::Vector6d titanInitialState =
                spice_interface::getBodyCartesianStateAtEpoch( "Titan", "Saturn", "ECLIPJ2000", "None", arcStartTimes.at( j ) );
        titanInitialState( 0 ) += 1.0E3;
        Eigen::VectorXd initialState( 12 );
        initialState << titanInitialState, cassiniInitialState;
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies,
                        AccelerationMap( ),
                        bodiesToPropagate,
                        initialState,
                        arcStartTimes.at( j ),
                        rungeKuttaFixedStepSettings( 10.0, rungeKuttaFehlberg78 ),
                        propagationTimeTerminationSettings( arcStartTimes.at( j ) + 3600.0 ) );
        translationalSettings->resetAccelerationModelsMap( accelerationMap, bodies );
        propagatorSettingsList.push_back( translationalSettings );
    }
    std::shared_ptr< MultiArcPropagatorSettings< double > > propagatorSettings = std::make_shared< MultiArcPropagatorSettings< double > >(
            propagatorSettingsList, false, std::make_shared< MultiArcPropagatorProcessingSettings >( false, true ) );
    propagatorSettings->setNumberOfThreads( numberOfThreads );
    MultiArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings, true );

    std::vector< Eigen::Vector6d > cassiniStates;
    for( unsigned int j = 0; j < arcStartTimes.size( ); j++ )
    {
        cassiniStates.push_back( bodies.at( "Cassini" )->getEphemeris( )->getCartesianState( arcStartTimes.at( j ) + 1800.0 ) );
    }
    return cassiniStates;
}

//! Test whether the propagated states set in the bodies are equal for multi-threaded and single-threaded propagation, when a central
//! body is itself propagated
BOOST_AUTO_TEST_CASE( testParallelPropagationWithPropagatedCentralBody )
{
    spice_interface::loadStandardSpiceKernels( );
    std::vector< Eigen::Vector6d > serialStates = propagateWithPropagatedCentralBody( 1 );
    std::vector< Eigen::Vector6d > parallelStates = propagateWithPropagatedCentralBody( 2 );
    for( unsigned int j = 0; j < serialStates.size( ); j++ )
    {
        double positionDifference = ( parallelStates.at( j ) - serialStates.at( j ) ).segment( 0, 3 ).norm( );
        BOOST_CHECK( !serialStates.at( j ).hasNaN( ) );
        BOOST_CHECK_SMALL( positionDifference, 1.0E-3 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
