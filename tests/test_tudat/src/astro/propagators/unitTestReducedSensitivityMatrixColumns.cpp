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
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
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

BOOST_AUTO_TEST_SUITE( test_reduced_sensitivity_matrix_columns )

enum ArcTestType { single_arc_test, multi_arc_test, hybrid_arc_test };

// Reference epoch and duration of a single (multi-arc) propagation arc
const double initialTime = 1.0E7;
const double arcDuration = 43200.0;

//! Output of a propagation of the variational equations, used to compare propagations with and without column reduction
struct ReducedColumnsTestOutput {
    // Propagated sensitivity matrix columns, for each variational equations object
    std::vector< std::vector< int > > propagatedColumns;

    // Expected propagated sensitivity matrix columns, for each variational equations object
    std::vector< std::vector< int > > expectedColumns;

    // Full state transition and sensitivity matrices at test epochs
    std::vector< Eigen::MatrixXd > fullMatrices;

    // Sensitivity matrix history (single arc only)
    std::map< double, Eigen::MatrixXd > sensitivityHistory;
};

//! Function to create the expected propagated columns, from the indices of the (arc-wise) empirical acceleration arcs that are active
std::vector< int > getExpectedColumns( const int gravitationalParameterColumn,
                                       const int empiricalAccelerationStartColumn,
                                       const std::vector< int >& activeEmpiricalAccelerationArcs )
{
    std::vector< int > expectedColumns = { gravitationalParameterColumn };
    for( int arc : activeEmpiricalAccelerationArcs )
    {
        for( int i = 0; i < 3; i++ )
        {
            expectedColumns.push_back( empiricalAccelerationStartColumn + 3 * arc + i );
        }
    }
    std::sort( expectedColumns.begin( ), expectedColumns.end( ) );
    return expectedColumns;
}

//! Function to propagate the variational equations of a Mars orbiter (and Mars, for hybrid arc) with arc-wise empirical
//! accelerations, with (multi-arc) state arcs [t0+i*D, t0+(i+1)*D] (single-arc: [t0, t0+1.5*D])
ReducedColumnsTestOutput propagateWithArcWiseEmpiricalAccelerations(
        const ArcTestType arcType,
        const std::vector< double >& empiricalAccelerationArcStarts,
        const std::vector< std::vector< int > >& activeEmpiricalAccelerationArcsPerStateArc,
        const bool reduceColumns )
{
    // Create bodies, and an Earth ground station (of which the position cannot influence the dynamics)
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Sun", "Mars", "Jupiter", "Earth" }, initialTime - 1.0E4, initialTime + 4.0 * arcDuration );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    createGroundStation( bodies.at( "Earth" ), "Station", ( Eigen::Vector3d( ) << 6378.0E3, 0.0, 0.0 ).finished( ) );
    bodies.createEmptyBody( "Orbiter" );
    bodies.at( "Orbiter" )->setConstantBodyMass( 1.0E3 );
    bodies.at( "Orbiter" )
            ->setEphemeris(
                    std::make_shared< MultiArcEphemeris >( std::map< double, std::shared_ptr< Ephemeris > >( ), "Mars", "ECLIPJ2000" ) );
    bodies.processBodyFrameDefinitions( );

    // Create orbiter propagation settings
    SelectedAccelerationMap orbiterAccelerationSettings;
    orbiterAccelerationSettings[ "Orbiter" ][ "Mars" ] = { pointMassGravityAcceleration( ), empiricalAcceleration( ) };
    orbiterAccelerationSettings[ "Orbiter" ][ "Sun" ] = { pointMassGravityAcceleration( ) };
    AccelerationMap orbiterAccelerations = createAccelerationModelsMap( bodies, orbiterAccelerationSettings, { "Orbiter" }, { "Mars" } );

    double marsGravitationalParameter = bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    std::vector< double > stateArcStarts = { initialTime, initialTime + arcDuration, initialTime + 2.0 * arcDuration };
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > orbiterArcSettings;
    for( unsigned int i = 0; i < stateArcStarts.size( ); i++ )
    {
        Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 6000.0E3, 0.05, 1.5, 4.1 - i, 0.4 + i, 2.4 + i ).finished( );
        orbiterArcSettings.push_back( translationalStatePropagatorSettings< double >(
                { "Mars" },
                orbiterAccelerations,
                { "Orbiter" },
                convertKeplerianToCartesianElements( initialKeplerianState, marsGravitationalParameter ),
                stateArcStarts.at( i ),
                rungeKutta4Settings( 60.0 ),
                propagationTimeTerminationSettings( stateArcStarts.at( i ) + ( arcType == single_arc_test ? 1.5 : 1.0 ) * arcDuration ) ) );
    }

    std::shared_ptr< PropagatorSettings< double > > propagatorSettings;
    if( arcType == single_arc_test )
    {
        propagatorSettings = orbiterArcSettings.at( 0 );
    }
    else
    {
        std::shared_ptr< MultiArcPropagatorSettings< double > > multiArcSettings =
                std::make_shared< MultiArcPropagatorSettings< double > >( orbiterArcSettings );
        propagatorSettings = multiArcSettings;
        if( arcType == hybrid_arc_test )
        {
            // Create single-arc propagation of Mars
            SelectedAccelerationMap marsAccelerationSettings;
            marsAccelerationSettings[ "Mars" ][ "Sun" ] = { pointMassGravityAcceleration( ) };
            marsAccelerationSettings[ "Mars" ][ "Jupiter" ] = { pointMassGravityAcceleration( ) };
            std::shared_ptr< SingleArcPropagatorSettings< double > > marsSettings = translationalStatePropagatorSettings< double >(
                    { "SSB" },
                    createAccelerationModelsMap( bodies, marsAccelerationSettings, { "Mars" }, { "SSB" } ),
                    { "Mars" },
                    getInitialStatesOfBodies( { "Mars" }, { "SSB" }, bodies, initialTime ),
                    initialTime,
                    rungeKutta4Settings( 60.0 ),
                    propagationTimeTerminationSettings( initialTime + 3.0 * arcDuration ) );
            propagatorSettings = std::make_shared< HybridArcPropagatorSettings< double > >( marsSettings, multiArcSettings );
            propagatorSettings->getOutputSettingsBase( )->setIntegratedResult( true );
        }
    }

    // Create parameters: initial states, Sun gravitational parameter, arc-wise empirical accelerations and station position
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies, stateArcStarts );
    parameterSettings.push_back( gravitationalParameter( "Sun" ) );
    parameterSettings.push_back( constantArcWiseEmpiricalAccelerationMagnitudes( "Orbiter", "Mars", empiricalAccelerationArcStarts ) );
    parameterSettings.push_back( groundStationPosition( "Earth", "Station" ) );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterSettings, bodies, propagatorSettings );
    parametersToEstimate->setReduceSensitivityMatrixColumns( reduceColumns );

    int initialStateSize = parametersToEstimate->getInitialDynamicalStateParameterSize( );
    int gravitationalParameterColumn =
            parametersToEstimate->getIndicesForParameterType( std::make_pair( gravitational_parameter, std::make_pair( "Sun", "" ) ) )
                    .at( 0 )
                    .first -
            initialStateSize;
    int empiricalAccelerationStartColumn =
            parametersToEstimate
                    ->getIndicesForParameterType(
                            std::make_pair( arc_wise_empirical_acceleration_coefficients, std::make_pair( "Orbiter", "Mars" ) ) )
                    .at( 0 )
                    .first -
            initialStateSize;

    // Propagate variational equations, and retrieve propagated columns for each variational equations object
    ReducedColumnsTestOutput output;
    std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface;
    std::vector< std::shared_ptr< VariationalEquations > > variationalEquationsObjects;
    if( arcType == single_arc_test )
    {
        SingleArcVariationalEquationsSolver< double, double > solver( bodies, propagatorSettings, parametersToEstimate );
        stateTransitionInterface = solver.getStateTransitionMatrixInterface( );
        variationalEquationsObjects.push_back( solver.getVariationalEquationsObject( ) );
        output.sensitivityHistory = solver.getSingleArcVariationalPropagationResults( )->getSensitivitySolution( );
    }
    else
    {
        std::shared_ptr< MultiArcVariationalEquationsSolver< double, double > > multiArcSolver;
        std::shared_ptr< VariationalEquationsSolver< double, double > > solver;
        if( arcType == multi_arc_test )
        {
            multiArcSolver = std::make_shared< MultiArcVariationalEquationsSolver< double, double > >(
                    bodies,
                    std::dynamic_pointer_cast< MultiArcPropagatorSettings< double > >( propagatorSettings ),
                    parametersToEstimate,
                    true );
            solver = multiArcSolver;
        }
        else
        {
            std::shared_ptr< HybridArcVariationalEquationsSolver< double, double > > hybridArcSolver =
                    std::make_shared< HybridArcVariationalEquationsSolver< double, double > >(
                            bodies,
                            std::dynamic_pointer_cast< HybridArcPropagatorSettings< double > >( propagatorSettings ),
                            parametersToEstimate,
                            true );
            multiArcSolver = hybridArcSolver->getMultiArcSolver( );
            variationalEquationsObjects.push_back( hybridArcSolver->getSingleArcSolver( )->getVariationalEquationsObject( ) );
            output.expectedColumns.push_back( { gravitationalParameterColumn } );
            solver = hybridArcSolver;
        }
        stateTransitionInterface = solver->getStateTransitionMatrixInterface( );
        for( auto stateDerivative : multiArcSolver->getDynamicsStateDerivatives( ) )
        {
            variationalEquationsObjects.push_back( stateDerivative->getVariationalEquationsCalculator( ) );
        }
    }

    for( unsigned int i = 0; i < activeEmpiricalAccelerationArcsPerStateArc.size( ); i++ )
    {
        output.expectedColumns.push_back( getExpectedColumns(
                gravitationalParameterColumn, empiricalAccelerationStartColumn, activeEmpiricalAccelerationArcsPerStateArc.at( i ) ) );
    }
    for( auto variationalEquations : variationalEquationsObjects )
    {
        output.propagatedColumns.push_back( variationalEquations->getReducedSensitivityMatrixColumns( ) );
    }

    // Retrieve full state transition and sensitivity matrices at epochs inside each state arc
    for( double testTime : { 0.3, 0.7, 1.3, 1.7, 2.3, 2.7 } )
    {
        if( arcType != single_arc_test || testTime < 1.5 )
        {
            output.fullMatrices.push_back(
                    stateTransitionInterface->getFullCombinedStateTransitionAndSensitivityMatrix( initialTime + testTime * arcDuration ) );
        }
    }
    return output;
}

//! Function to check that each column of a matrix is equal to that of a reference matrix, to within a relative tolerance
void checkMatrixColumnsClose( const Eigen::MatrixXd& matrix, const Eigen::MatrixXd& referenceMatrix )
{
    BOOST_CHECK_EQUAL( matrix.rows( ), referenceMatrix.rows( ) );
    BOOST_CHECK_EQUAL( matrix.cols( ), referenceMatrix.cols( ) );
    for( int i = 0; i < referenceMatrix.cols( ); i++ )
    {
        BOOST_CHECK_SMALL( ( matrix.col( i ) - referenceMatrix.col( i ) ).norm( ), 1.0E-12 * referenceMatrix.col( i ).norm( ) + 1.0E-300 );
    }
}

BOOST_AUTO_TEST_CASE( testPropagationTimeInterval )
{
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            translationalStatePropagatorSettings< double >( { "Earth" },
                                                            AccelerationMap( ),
                                                            { "Vehicle" },
                                                            Eigen::Vector6d::Zero( ),
                                                            initialTime,
                                                            rungeKutta4Settings( 60.0 ),
                                                            propagationTimeTerminationSettings( initialTime - arcDuration ) );

    // Backward propagation
    std::pair< double, double > timeInterval = getPropagationTimeInterval< double, double >( propagatorSettings );
    BOOST_CHECK_EQUAL( timeInterval.first, initialTime - arcDuration );
    BOOST_CHECK_EQUAL( timeInterval.second, initialTime );

    // Hybrid termination, terminating on first time condition
    propagatorSettings->resetTerminationSettings(
            propagationHybridTerminationSettings( { propagationCPUTimeTerminationSettings( 10.0 ),
                                                    propagationTimeTerminationSettings( initialTime + 2.0 * arcDuration ),
                                                    propagationTimeTerminationSettings( initialTime + arcDuration ) },
                                                  true ) );
    timeInterval = getPropagationTimeInterval< double, double >( propagatorSettings );
    BOOST_CHECK_EQUAL( timeInterval.first, initialTime );
    BOOST_CHECK_EQUAL( timeInterval.second, initialTime + arcDuration );

    // Termination time not known before propagation
    propagatorSettings->resetTerminationSettings( propagationCPUTimeTerminationSettings( 10.0 ) );
    BOOST_CHECK( std::isnan( getPropagationTimeInterval< double, double >( propagatorSettings ).second ) );
}

BOOST_AUTO_TEST_CASE( testReducedSensitivityMatrixColumns )
{
    spice_interface::loadStandardSpiceKernels( );

    // Empirical acceleration arc starts (in units of arc duration, w.r.t. initial time) for: arcs coinciding with state arcs,
    // state arcs integer multiple of empirical acceleration arcs, state arcs not an integer multiple of empirical acceleration arcs
    std::vector< std::vector< double > > empiricalAccelerationArcStarts = { { 0.0, 1.0, 2.0 },
                                                                            { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 },
                                                                            { -1.0, 0.0, 0.7, 1.4, 2.1, 2.8 } };

    // Active empirical acceleration arcs in single-arc propagation [0, 1.5], and in multi-arc propagation [0,1], [1,2], [2,3]
    std::vector< std::vector< std::vector< int > > > activeArcsSingleArc = { { { 0, 1 } }, { { 0, 1, 2 } }, { { 1, 2, 3 } } };
    std::vector< std::vector< std::vector< int > > > activeArcsMultiArc = { { { 0 }, { 1 }, { 2 } },
                                                                            { { 0, 1 }, { 2, 3 }, { 4, 5 } },
                                                                            { { 1, 2 }, { 2, 3 }, { 3, 4, 5 } } };

    for( ArcTestType arcType : { single_arc_test, multi_arc_test, hybrid_arc_test } )
    {
        for( unsigned int testCase = 0; testCase < empiricalAccelerationArcStarts.size( ); testCase++ )
        {
            std::vector< double > arcStarts;
            for( double arcStart : empiricalAccelerationArcStarts.at( testCase ) )
            {
                arcStarts.push_back( initialTime + arcStart * arcDuration );
            }
            std::vector< std::vector< int > > activeArcs =
                    ( arcType == single_arc_test ) ? activeArcsSingleArc.at( testCase ) : activeArcsMultiArc.at( testCase );

            ReducedColumnsTestOutput reducedOutput = propagateWithArcWiseEmpiricalAccelerations( arcType, arcStarts, activeArcs, true );
            ReducedColumnsTestOutput fullOutput = propagateWithArcWiseEmpiricalAccelerations( arcType, arcStarts, activeArcs, false );

            // Check propagated columns (station position and inactive empirical acceleration arcs omitted)
            BOOST_CHECK_EQUAL( reducedOutput.propagatedColumns.size( ), reducedOutput.expectedColumns.size( ) );
            for( unsigned int i = 0; i < reducedOutput.propagatedColumns.size( ); i++ )
            {
                BOOST_CHECK( reducedOutput.propagatedColumns.at( i ) == reducedOutput.expectedColumns.at( i ) );
            }

            // Check that full matrices are equal to those without column reduction
            for( unsigned int i = 0; i < fullOutput.fullMatrices.size( ); i++ )
            {
                checkMatrixColumnsClose( reducedOutput.fullMatrices.at( i ), fullOutput.fullMatrices.at( i ) );
            }

            // Check full sensitivity matrix history, except at final epoch, where (without column reduction) the last integrator
            // stage evaluates an empirical acceleration arc that starts at the final epoch
            BOOST_CHECK_EQUAL( reducedOutput.sensitivityHistory.size( ), fullOutput.sensitivityHistory.size( ) );
            for( auto it : fullOutput.sensitivityHistory )
            {
                if( it.first < fullOutput.sensitivityHistory.rbegin( )->first )
                {
                    checkMatrixColumnsClose( reducedOutput.sensitivityHistory.at( it.first ), it.second );
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
