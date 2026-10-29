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

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/desaturationDeltaV.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::estimatable_parameters;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::ephemerides;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_astrodynamics;
using namespace tudat::observation_models;

BOOST_AUTO_TEST_SUITE( test_reduced_sensitivity_matrix_columns )

enum ArcTestType { single_arc_test, non_sequential_single_arc_test, multi_arc_test, hybrid_arc_test };

// Reference epoch and duration of a single (multi-arc) propagation arc
const double initialTime = 1.0E7;
const double arcDuration = 43200.0;

//! Settings for a propagation of a Mars orbiter (and Mars, for hybrid arc) with arc-wise empirical accelerations
template< typename StateScalarType >
struct ReducedColumnsTestSetup {
    SystemOfBodies bodies;
    std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings;
    std::shared_ptr< EstimatableParameterSet< StateScalarType > > parametersToEstimate;
};

//! Function to create the settings for a propagation of a Mars orbiter (and Mars, for hybrid arc) with arc-wise empirical
//! accelerations, with (multi-arc) state arcs [t0+i*D, t0+(i+1)*D] (single-arc: [t0, t0+1.5*D]). Estimated parameters are the
//! initial states, Sun gravitational parameter, arc-wise empirical accelerations and an Earth station position.
template< typename StateScalarType >
ReducedColumnsTestSetup< StateScalarType > createTestSetup( const ArcTestType arcType,
                                                            const std::vector< double >& empiricalAccelerationArcStarts,
                                                            const bool reduceColumns )
{
    // Create bodies, and an Earth ground station (of which the position cannot influence the dynamics)
    ReducedColumnsTestSetup< StateScalarType > setup;
    setup.bodies = createSystemOfBodies(
            getDefaultBodySettings( { "Sun", "Mars", "Jupiter", "Earth" }, initialTime - 1.0E4, initialTime + 4.0 * arcDuration ) );
    SystemOfBodies& bodies = setup.bodies;
    createGroundStation( bodies.at( "Earth" ), "Station", ( Eigen::Vector3d( ) << 6378.0E3, 0.0, 0.0 ).finished( ) );
    bodies.createEmptyBody( "Orbiter" );
    bodies.at( "Orbiter" )->setConstantBodyMass( 1.0E3 );
    if( arcType == multi_arc_test || arcType == hybrid_arc_test )
    {
        bodies.at( "Orbiter" )
                ->setEphemeris( std::make_shared< MultiArcEphemeris >(
                        std::map< double, std::shared_ptr< Ephemeris > >( ), "Mars", "ECLIPJ2000" ) );
    }
    else
    {
        bodies.at( "Orbiter" )
                ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ),
                        "Mars",
                        "ECLIPJ2000" ) );
    }
    bodies.processBodyFrameDefinitions( );

    // Create orbiter propagation settings
    SelectedAccelerationMap orbiterAccelerationSettings;
    orbiterAccelerationSettings[ "Orbiter" ][ "Mars" ] = { pointMassGravityAcceleration( ), empiricalAcceleration( ) };
    orbiterAccelerationSettings[ "Orbiter" ][ "Sun" ] = { pointMassGravityAcceleration( ) };
    AccelerationMap orbiterAccelerations = createAccelerationModelsMap( bodies, orbiterAccelerationSettings, { "Orbiter" }, { "Mars" } );

    double marsGravitationalParameter = bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > orbiterArcSettings;
    for( unsigned int i = 0; i < 3; i++ )
    {
        double arcInitialTime = initialTime + i * arcDuration;
        std::shared_ptr< PropagationTerminationSettings > terminationSettings =
                propagationTimeTerminationSettings( arcInitialTime + arcDuration );
        if( arcType == single_arc_test )
        {
            terminationSettings = propagationTimeTerminationSettings( initialTime + 1.5 * arcDuration );
        }
        else if( arcType == non_sequential_single_arc_test )
        {
            arcInitialTime = initialTime + 0.75 * arcDuration;
            terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
                    propagationTimeTerminationSettings( initialTime + 1.5 * arcDuration ),
                    propagationTimeTerminationSettings( initialTime ) );
        }

        Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 6000.0E3, 0.05, 1.5, 4.1 - i, 0.4 + i, 2.4 + i ).finished( );
        orbiterArcSettings.push_back( translationalStatePropagatorSettings< StateScalarType >(
                { "Mars" },
                orbiterAccelerations,
                { "Orbiter" },
                convertKeplerianToCartesianElements( initialKeplerianState, marsGravitationalParameter )
                        .template cast< StateScalarType >( ),
                arcInitialTime,
                rungeKutta4Settings( 60.0 ),
                terminationSettings ) );
    }

    if( arcType == single_arc_test || arcType == non_sequential_single_arc_test )
    {
        setup.propagatorSettings = orbiterArcSettings.at( 0 );
    }
    else
    {
        std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > multiArcSettings =
                std::make_shared< MultiArcPropagatorSettings< StateScalarType > >( orbiterArcSettings );
        setup.propagatorSettings = multiArcSettings;
        if( arcType == hybrid_arc_test )
        {
            // Create single-arc propagation of Mars
            SelectedAccelerationMap marsAccelerationSettings;
            marsAccelerationSettings[ "Mars" ][ "Sun" ] = { pointMassGravityAcceleration( ) };
            marsAccelerationSettings[ "Mars" ][ "Jupiter" ] = { pointMassGravityAcceleration( ) };
            std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > marsSettings =
                    translationalStatePropagatorSettings< StateScalarType >(
                            { "SSB" },
                            createAccelerationModelsMap( bodies, marsAccelerationSettings, { "Mars" }, { "SSB" } ),
                            { "Mars" },
                            getInitialStatesOfBodies( { "Mars" }, { "SSB" }, bodies, initialTime ).template cast< StateScalarType >( ),
                            initialTime,
                            rungeKutta4Settings( 60.0 ),
                            propagationTimeTerminationSettings( initialTime + 3.0 * arcDuration ) );
            setup.propagatorSettings = std::make_shared< HybridArcPropagatorSettings< StateScalarType > >( marsSettings, multiArcSettings );
        }
    }
    setup.propagatorSettings->getOutputSettingsBase( )->setIntegratedResult( true );

    // Create parameters
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings = getInitialStateParameterSettings< StateScalarType >(
            setup.propagatorSettings, bodies, { initialTime, initialTime + arcDuration, initialTime + 2.0 * arcDuration } );
    parameterSettings.push_back( gravitationalParameter( "Sun" ) );
    parameterSettings.push_back( constantArcWiseEmpiricalAccelerationMagnitudes( "Orbiter", "Mars", empiricalAccelerationArcStarts ) );
    parameterSettings.push_back( groundStationPosition( "Earth", "Station" ) );
    setup.parametersToEstimate = createParametersToEstimate< StateScalarType >( parameterSettings, bodies, setup.propagatorSettings );
    setup.parametersToEstimate->setReduceSensitivityMatrixColumns( reduceColumns );
    return setup;
}

//! Output of a propagation of the variational equations, used to compare propagations with and without column reduction
struct VariationalEquationsTestOutput {
    // Propagated sensitivity matrix columns, for each variational equations object (empty if no reduction)
    std::vector< std::vector< int > > propagatedColumns;

    // Full state transition and sensitivity matrices at test epochs
    std::vector< Eigen::MatrixXd > fullMatrices;

    // Sensitivity matrix history (single arc only)
    std::map< double, Eigen::MatrixXd > sensitivityHistory;
};

//! Function to propagate the variational equations for given test settings
template< typename StateScalarType >
VariationalEquationsTestOutput propagateVariationalEquations( const ArcTestType arcType,
                                                              const std::vector< double >& empiricalAccelerationArcStarts,
                                                              const bool reduceColumns )
{
    ReducedColumnsTestSetup< StateScalarType > setup =
            createTestSetup< StateScalarType >( arcType, empiricalAccelerationArcStarts, reduceColumns );

    // Propagate variational equations, and retrieve propagated columns for each variational equations object
    VariationalEquationsTestOutput output;
    std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface;
    std::vector< std::shared_ptr< VariationalEquations > > variationalEquationsObjects;
    if( arcType == single_arc_test || arcType == non_sequential_single_arc_test )
    {
        SingleArcVariationalEquationsSolver< StateScalarType, double > solver(
                setup.bodies, setup.propagatorSettings, setup.parametersToEstimate );
        stateTransitionInterface = solver.getStateTransitionMatrixInterface( );
        variationalEquationsObjects.push_back( solver.getVariationalEquationsObject( ) );
        output.sensitivityHistory = solver.getSingleArcVariationalPropagationResults( )->getSensitivitySolution( );
    }
    else
    {
        std::shared_ptr< MultiArcVariationalEquationsSolver< StateScalarType, double > > multiArcSolver;
        if( arcType == multi_arc_test )
        {
            multiArcSolver = std::make_shared< MultiArcVariationalEquationsSolver< StateScalarType, double > >(
                    setup.bodies,
                    std::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType > >( setup.propagatorSettings ),
                    setup.parametersToEstimate,
                    true );
            stateTransitionInterface = multiArcSolver->getStateTransitionMatrixInterface( );
        }
        else
        {
            std::shared_ptr< HybridArcVariationalEquationsSolver< StateScalarType, double > > hybridArcSolver =
                    std::make_shared< HybridArcVariationalEquationsSolver< StateScalarType, double > >(
                            setup.bodies,
                            std::dynamic_pointer_cast< HybridArcPropagatorSettings< StateScalarType > >( setup.propagatorSettings ),
                            setup.parametersToEstimate,
                            true );
            stateTransitionInterface = hybridArcSolver->getStateTransitionMatrixInterface( );
            multiArcSolver = hybridArcSolver->getMultiArcSolver( );
            variationalEquationsObjects.push_back( hybridArcSolver->getSingleArcSolver( )->getVariationalEquationsObject( ) );
        }
        for( auto stateDerivative : multiArcSolver->getDynamicsStateDerivatives( ) )
        {
            variationalEquationsObjects.push_back( stateDerivative->getVariationalEquationsCalculator( ) );
        }
    }

    for( auto variationalEquations : variationalEquationsObjects )
    {
        output.propagatedColumns.push_back( variationalEquations->getReducedSensitivityMatrixColumns( ) );
    }

    // Retrieve full state transition and sensitivity matrices at epochs inside each state arc
    for( double testTime : { 0.3, 0.7, 1.3, 1.7, 2.3, 2.7 } )
    {
        if( ( arcType != single_arc_test && arcType != non_sequential_single_arc_test ) || testTime < 1.5 )
        {
            output.fullMatrices.push_back(
                    stateTransitionInterface->getFullCombinedStateTransitionAndSensitivityMatrix( initialTime + testTime * arcDuration ) );
        }
    }
    return output;
}

//! Function to estimate the parameters for given test settings, from simulated orbiter position observations
std::shared_ptr< EstimationOutput< double > > estimateParameters( const ArcTestType arcType,
                                                                  const std::vector< double >& empiricalAccelerationArcStarts,
                                                                  const bool reduceColumns )
{
    ReducedColumnsTestSetup< double > setup = createTestSetup< double >( arcType, empiricalAccelerationArcStarts, reduceColumns );

    LinkDefinition linkEnds;
    linkEnds[ observed_body ] = std::make_pair< std::string, std::string >( "Orbiter", "" );
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            setup.bodies,
            setup.parametersToEstimate,
            { std::make_shared< ObservationModelSettings >( position_observable, linkEnds ) },
            setup.propagatorSettings );

    // Simulate observations inside each state arc
    std::vector< double > observationTimes;
    for( double observationTime = initialTime + 300.0; observationTime < initialTime + 3.0 * arcDuration; observationTime += 600.0 )
    {
        if( ( arcType != single_arc_test || observationTime < initialTime + 1.5 * arcDuration ) &&
            std::fmod( observationTime - initialTime, arcDuration ) < arcDuration - 300.0 )
        {
            observationTimes.push_back( observationTime );
        }
    }
    std::shared_ptr< ObservationCollection<> > observations =
            simulateObservations< double, double >( { std::make_shared< TabulatedObservationSimulationSettings<> >(
                                                            position_observable, linkEnds, observationTimes, observed_body ) },
                                                    orbitDeterminationManager.getObservationSimulators( ),
                                                    setup.bodies );

    // Estimate parameters, with a priori covariance to regularize parameters that are not observed (e.g. station position)
    int numberOfParameters = setup.parametersToEstimate->getParameterSetSize( );
    std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( observations,
                                                                   Eigen::MatrixXd::Identity( numberOfParameters, numberOfParameters ),
                                                                   std::make_shared< EstimationConvergenceChecker >( 1 ) );
    estimationInput->defineEstimationSettings( true, true, true, false );
    return orbitDeterminationManager.estimateParameters( estimationInput );
}

//! Function to check that each column of a matrix is equal to that of a reference matrix, to within a relative tolerance
void checkMatrixColumnsClose( const Eigen::MatrixXd& matrix, const Eigen::MatrixXd& referenceMatrix, const double tolerance = 1.0E-12 )
{
    BOOST_CHECK_EQUAL( matrix.rows( ), referenceMatrix.rows( ) );
    BOOST_CHECK_EQUAL( matrix.cols( ), referenceMatrix.cols( ) );
    for( int i = 0; i < referenceMatrix.cols( ); i++ )
    {
        BOOST_CHECK_SMALL( ( matrix.col( i ) - referenceMatrix.col( i ) ).norm( ),
                           tolerance * referenceMatrix.col( i ).norm( ) + 1.0E-300 );
    }
}

BOOST_AUTO_TEST_CASE( testParameterActivityDetermination )
{
    // Check propagation time interval for backward propagation, and termination at unknown time
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            translationalStatePropagatorSettings< double >( { "Earth" },
                                                            AccelerationMap( ),
                                                            { "Vehicle" },
                                                            Eigen::Vector6d::Zero( ),
                                                            initialTime,
                                                            rungeKutta4Settings( 60.0 ),
                                                            propagationTimeTerminationSettings( initialTime - arcDuration ) );
    std::pair< double, double > timeInterval = getPropagationTimeInterval< double, double >( propagatorSettings );
    BOOST_CHECK_EQUAL( timeInterval.first, initialTime - arcDuration );
    BOOST_CHECK_EQUAL( timeInterval.second, initialTime );

    propagatorSettings->resetTerminationSettings( propagationCPUTimeTerminationSettings( 10.0 ) );
    BOOST_CHECK( std::isnan( getPropagationTimeInterval< double, double >( propagatorSettings ).second ) );

    // Check active entries of desaturation maneuvers, which act from 70 to 130 s and 970 to 1030 s (conservatively including the
    // maneuver start time)
    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > desaturationParameter = std::make_shared< DesaturationDeltaV >(
            std::make_shared< propulsion::MomentumWheelDesaturationThrustAcceleration >(
                    std::vector< double >( { 100.0, 1000.0 } ), std::vector< Eigen::Vector3d >( 2, Eigen::Vector3d::Zero( ) ), 60.0, 10.0 ),
            "Vehicle" );
    BOOST_CHECK( getParameterEntriesActiveInTimeInterval( desaturationParameter, std::make_pair( 0.0, 69.0 ) ) ==
                 std::vector< bool >( 6, false ) );
    BOOST_CHECK( getParameterEntriesActiveInTimeInterval( desaturationParameter, std::make_pair( 0.0, 70.0 ) ) ==
                 std::vector< bool >( { true, true, true, false, false, false } ) );
    BOOST_CHECK( getParameterEntriesActiveInTimeInterval( desaturationParameter, std::make_pair( 130.0, 1500.0 ) ) ==
                 std::vector< bool >( { false, false, false, true, true, true } ) );
    BOOST_CHECK( getParameterEntriesActiveInTimeInterval( desaturationParameter, std::make_pair( TUDAT_NAN, TUDAT_NAN ) ) ==
                 std::vector< bool >( 6, true ) );
}

// Empirical acceleration arc starts (in units of arc duration, w.r.t. initial time) for: arcs coinciding with state arcs,
// state arcs integer multiple of empirical acceleration arcs, state arcs not an integer multiple of empirical acceleration arcs
const std::vector< std::vector< double > > empiricalAccelerationArcStarts = { { 0.0, 1.0, 2.0 },
                                                                              { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 },
                                                                              { -1.0, 0.0, 0.7, 1.4, 2.1, 2.8 } };

//! Function to retrieve the empirical acceleration arc start times for a test case
std::vector< double > getEmpiricalAccelerationArcStarts( const unsigned int testCase )
{
    std::vector< double > arcStarts;
    for( double arcStart : empiricalAccelerationArcStarts.at( testCase ) )
    {
        arcStarts.push_back( initialTime + arcStart * arcDuration );
    }
    return arcStarts;
}

//! Function to retrieve the expected propagated sensitivity matrix columns for each variational equations object for a test case.
//! Sensitivity matrix columns are: Sun gravitational parameter (0), empirical accelerations (3 per arc), station position.
std::vector< std::vector< int > > getExpectedColumns( const ArcTestType arcType, const unsigned int testCase )
{
    // Active empirical acceleration arcs in single-arc propagation [0, 1.5], and in multi-arc propagation [0,1], [1,2], [2,3]. An arc
    // that starts at the final time of a propagation is active, since it is used in the final integrator stage.
    std::vector< std::vector< std::vector< int > > > activeArcsSingleArc = { { { 0, 1 } }, { { 0, 1, 2, 3 } }, { { 1, 2, 3 } } };
    std::vector< std::vector< std::vector< int > > > activeArcsMultiArc = { { { 0, 1 }, { 1, 2 }, { 2 } },
                                                                            { { 0, 1, 2 }, { 2, 3, 4 }, { 4, 5 } },
                                                                            { { 1, 2 }, { 2, 3 }, { 3, 4, 5 } } };
    std::vector< std::vector< int > > activeArcs = ( arcType == single_arc_test || arcType == non_sequential_single_arc_test )
            ? activeArcsSingleArc.at( testCase )
            : activeArcsMultiArc.at( testCase );

    // Single-arc Mars propagation (hybrid arc) only depends on Sun gravitational parameter
    std::vector< std::vector< int > > expectedColumns;
    if( arcType == hybrid_arc_test )
    {
        expectedColumns.push_back( { 0 } );
    }
    for( unsigned int i = 0; i < activeArcs.size( ); i++ )
    {
        expectedColumns.push_back( { 0 } );
        for( int arc : activeArcs.at( i ) )
        {
            for( int j = 0; j < 3; j++ )
            {
                expectedColumns.back( ).push_back( 1 + 3 * arc + j );
            }
        }
    }
    return expectedColumns;
}

template< typename StateScalarType >
void testReducedVariationalEquations( )
{
    for( ArcTestType arcType : { single_arc_test, non_sequential_single_arc_test, multi_arc_test, hybrid_arc_test } )
    {
        for( unsigned int testCase = 0; testCase < empiricalAccelerationArcStarts.size( ); testCase++ )
        {
            VariationalEquationsTestOutput reducedOutput =
                    propagateVariationalEquations< StateScalarType >( arcType, getEmpiricalAccelerationArcStarts( testCase ), true );
            VariationalEquationsTestOutput fullOutput =
                    propagateVariationalEquations< StateScalarType >( arcType, getEmpiricalAccelerationArcStarts( testCase ), false );

            // Check propagated columns (station position and inactive empirical acceleration arcs omitted)
            BOOST_CHECK( reducedOutput.propagatedColumns == getExpectedColumns( arcType, testCase ) );

            // Check that full matrices are equal to those without column reduction
            BOOST_CHECK_EQUAL( reducedOutput.fullMatrices.size( ), fullOutput.fullMatrices.size( ) );
            for( unsigned int i = 0; i < fullOutput.fullMatrices.size( ); i++ )
            {
                checkMatrixColumnsClose( reducedOutput.fullMatrices.at( i ), fullOutput.fullMatrices.at( i ) );
            }

            // Check full sensitivity matrix history
            BOOST_CHECK_EQUAL( reducedOutput.sensitivityHistory.size( ), fullOutput.sensitivityHistory.size( ) );
            for( auto it : fullOutput.sensitivityHistory )
            {
                checkMatrixColumnsClose( reducedOutput.sensitivityHistory.at( it.first ), it.second );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testReducedSensitivityMatrixColumns )
{
    spice_interface::loadStandardSpiceKernels( );
    testReducedVariationalEquations< double >( );
#if ( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
    testReducedVariationalEquations< long double >( );
#endif
}

BOOST_AUTO_TEST_CASE( testEstimationWithReducedSensitivityMatrixColumns )
{
    spice_interface::loadStandardSpiceKernels( );
    for( ArcTestType arcType : { single_arc_test, multi_arc_test, hybrid_arc_test } )
    {
        for( unsigned int testCase = 0; testCase < empiricalAccelerationArcStarts.size( ); testCase++ )
        {
            std::shared_ptr< EstimationOutput< double > > reducedOutput =
                    estimateParameters( arcType, getEmpiricalAccelerationArcStarts( testCase ), true );
            std::shared_ptr< EstimationOutput< double > > fullOutput =
                    estimateParameters( arcType, getEmpiricalAccelerationArcStarts( testCase ), false );

            // Check that observation partials, covariance and estimate are equal to those without column reduction
            checkMatrixColumnsClose( reducedOutput->getUnnormalizedDesignMatrix( ), fullOutput->getUnnormalizedDesignMatrix( ) );
            checkMatrixColumnsClose(
                    reducedOutput->getUnnormalizedCovarianceMatrix( ), fullOutput->getUnnormalizedCovarianceMatrix( ), 1.0E-8 );
            Eigen::VectorXd estimateDifference = ( reducedOutput->parameterEstimate_ - fullOutput->parameterEstimate_ )
                                                         .cwiseQuotient( fullOutput->getFormalErrorVector( ) );
            BOOST_CHECK_SMALL( estimateDifference.cwiseAbs( ).maxCoeff( ), 1.0E-6 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
