/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <map>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/propagators/stateTransitionMatrixInterface.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/simulation/estimation_setup/interArcContinuityConstraint.h"
#include "tudat/simulation/estimation_setup/interArcStateContinuityConstraintSettings.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/estimatableParameterSettings.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_inter_arc_continuity_constraint )

using namespace tudat;
using namespace tudat::estimatable_parameters;
using namespace tudat::numerical_integrators;
using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;

//! Own the propagated model and the interfaces consumed by the constraint assembler. The fixture also records
//! arc bounds and matrix dimensions so tests do not infer layout details from hard-coded sizes.
struct TwoArcFixture {
    //! Explicit-arc access to the propagated variational matrices.
    std::shared_ptr< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< double > > stmInterface;
    //! Numerical state histories used to evaluate the left and right boundary states.
    std::shared_ptr< MultiArcDynamicsSimulator< double, double > > simulator;
    //! Arc-wise initial states with the parameter ordering used by the variational matrices.
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate;
    std::vector< double > arcStartTimes;
    std::vector< double > arcEndTimes;
    int fullStateTransitionSize = 0;
    int fullSensitivitySize = 0;
};

//! Build a minimal two-arc Earth/Sun setup with adjacent arcs (no overlap) and return the full fixture
//! (variational solver, dynamics simulator, parameter set, STM interface). The variational equations are
//! integrated on construction so the interpolators are populated.
TwoArcFixture buildTwoArcFixture(
        const std::vector< Eigen::Matrix< double, 6, 1 > >& arcInitialStatePerturbations = std::vector< Eigen::Matrix< double, 6, 1 > >( ),
        const bool useOverlappingArcs = false )
{
    // SPICE supplies independent initial states for both arcs. This deliberately permits a non-zero boundary
    // discrepancy, which is required to exercise the cost and right-hand-side terms.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames = { "Earth", "Sun" };

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 4.0E7;
    const double buffer = 3.6E5;

    // The ephemeris interval extends beyond all integration epochs so interpolation coverage cannot influence
    // the constraint-specific assertions.
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    SystemOfBodies bodies = createSystemOfBodies< double, double >( bodySettings );

    // A single point-mass acceleration keeps the fixture inexpensive while still generating non-trivial state
    // transition matrices over each long arc.
    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Earth" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    std::vector< std::string > bodiesToIntegrate = { "Earth" };
    std::vector< std::string > centralBodies = { "SSB" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    const double arcDuration = 1.0E7;
    TwoArcFixture fixture;
    // Adjacent arcs exercise shared-boundary lookup; overlapping arcs exercise interpolation at an interior
    // connection epoch. Both modes retain the same duration and parameter layout.
    if( useOverlappingArcs )
    {
        const double overlapDuration = 1.0E6;
        fixture.arcStartTimes = { initialEphemerisTime + 1.0E5, initialEphemerisTime + 1.0E5 + arcDuration - overlapDuration };
        fixture.arcEndTimes = { fixture.arcStartTimes[ 0 ] + arcDuration, fixture.arcStartTimes[ 1 ] + arcDuration };
    }
    else
    {
        fixture.arcStartTimes = { initialEphemerisTime + 1.0E5, initialEphemerisTime + 1.0E5 + arcDuration };
        fixture.arcEndTimes = { fixture.arcStartTimes[ 0 ] + arcDuration, fixture.arcStartTimes[ 1 ] + arcDuration };
    }

    std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKutta4Settings< double >( 600.0 );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > propagatorSettingsList;
    for( unsigned int i = 0; i < fixture.arcStartTimes.size( ); ++i )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialState =
                getInitialStateOfBody< double, double >( "Earth", "SSB", bodies, fixture.arcStartTimes[ i ] );
        // Symmetric initial-state perturbations are supplied only by the finite-difference test. Applying them
        // before constructing the propagator makes that derivative independent of the STM accessor under test.
        if( !arcInitialStatePerturbations.empty( ) )
        {
            initialState += arcInitialStatePerturbations.at( i );
        }
        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                centralBodies,
                accelerationModelMap,
                bodiesToIntegrate,
                initialState,
                fixture.arcStartTimes[ i ],
                integratorSettings,
                propagationTimeTerminationSettings( fixture.arcEndTimes[ i ] ) ) );
    }
    std::shared_ptr< MultiArcPropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double, double > >( propagatorSettingsList, false );

    // Estimate each arc's translational initial state, then propagate state and variational equations together.
    // This produces the exact manager interfaces used by production constraint assembly.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double, double >( propagatorSettings, bodies, fixture.arcStartTimes );

    fixture.parametersToEstimate = createParametersToEstimate< double >( parameterNames, bodies );

    auto variationalSolver = simulation_setup::createVariationalEquationsSolver< double, double >(
            bodies, propagatorSettings, fixture.parametersToEstimate, true );

    fixture.stmInterface = std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< double > >(
            variationalSolver->getStateTransitionMatrixInterface( ) );
    fixture.simulator =
            std::dynamic_pointer_cast< MultiArcDynamicsSimulator< double, double > >( variationalSolver->getDynamicsSimulatorBase( ) );

    // A failed cast means the fixture did not create a pure multi-arc system and all subsequent checks would be
    // testing a different estimator configuration.
    BOOST_REQUIRE( fixture.stmInterface != nullptr );
    BOOST_REQUIRE( fixture.simulator != nullptr );

    fixture.fullStateTransitionSize = fixture.stmInterface->getFullStateTransitionMatrixSize( );
    fixture.fullSensitivitySize = fixture.stmInterface->getFullSensitivityMatrixSize( );

    return fixture;
}

//! Assemble one full-state boundary discrepancy for finite-difference tests. Keeping this operation in one
//! helper ensures that the perturbed and nominal propagations use identical constraint settings and normalization.
Eigen::VectorXd getSinglePairDiscrepancy( const TwoArcFixture& fixture )
{
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    // Unit weights prevent the finite-difference signal from being masked while leaving the discrepancy itself
    // independent of the chosen penalty strength.
    auto settings = simulation_setup::fullStateContinuity( { "Earth" }, { { "Earth", { fixture.arcStartTimes[ 1 ] } } }, 1.0, 1.0, 1.0 );
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );
    // The assembler is the public path that defines discrepancy ordering, so returning its diagnostic avoids
    // duplicating the production left/right sign convention in the finite-difference test.
    auto contribution = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                                  fixture.parametersToEstimate,
                                                                                  fixture.simulator,
                                                                                  fixture.stmInterface,
                                                                                  columnNormalizationFactors,
                                                                                  totalParameterSize );
    return contribution.perPairDiscrepancies.at( 0 );
}

//! Build a two-body counterpart of buildTwoArcFixture. Earth and Mars occupy distinct state-row blocks, making
//! it possible to detect accidental use of the first body's state or variational rows when Mars is constrained.
TwoArcFixture buildTwoBodyTwoArcFixture( )
{
    // Each body's independent SPICE state makes row-selection errors visible in both the discrepancy and design
    // matrix, instead of producing two numerically similar blocks.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames = { "Earth", "Mars", "Sun" };

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 4.0E7;
    const double buffer = 3.6E5;

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Mars" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    SystemOfBodies bodies = createSystemOfBodies< double, double >( bodySettings );

    // Both bodies use the same simple force-model class but separate acceleration instances and state blocks.
    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Earth" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Mars" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    std::vector< std::string > bodiesToIntegrate = { "Earth", "Mars" };
    std::vector< std::string > centralBodies = { "SSB", "SSB" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    const double arcDuration = 1.0E7;
    TwoArcFixture fixture;
    // Overlap permits distinct, interior connection epochs for the two bodies in one settings object.
    const double overlapDuration = 1.0E6;
    fixture.arcStartTimes = { initialEphemerisTime + 1.0E5, initialEphemerisTime + 1.0E5 + arcDuration - overlapDuration };
    fixture.arcEndTimes = { fixture.arcStartTimes[ 0 ] + arcDuration, fixture.arcStartTimes[ 1 ] + arcDuration };

    std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKutta4Settings< double >( 600.0 );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > propagatorSettingsList;
    for( unsigned int i = 0; i < fixture.arcStartTimes.size( ); ++i )
    {
        // Concatenating Earth before Mars creates a known non-zero row offset for the Mars regression test.
        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialState( 12 );
        initialState.segment( 0, 6 ) = getInitialStateOfBody< double, double >( "Earth", "SSB", bodies, fixture.arcStartTimes[ i ] );
        initialState.segment( 6, 6 ) = getInitialStateOfBody< double, double >( "Mars", "SSB", bodies, fixture.arcStartTimes[ i ] );
        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                centralBodies,
                accelerationModelMap,
                bodiesToIntegrate,
                initialState,
                fixture.arcStartTimes[ i ],
                integratorSettings,
                propagationTimeTerminationSettings( fixture.arcEndTimes[ i ] ) ) );
    }
    std::shared_ptr< MultiArcPropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double, double > >( propagatorSettingsList, false );

    // Generate one arc-wise initial-state parameter per body and a multi-arc variational solution with matching
    // row metadata. The assembler must use that metadata rather than assuming six rows in total.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double, double >( propagatorSettings, bodies, fixture.arcStartTimes );
    fixture.parametersToEstimate = createParametersToEstimate< double >( parameterNames, bodies );

    auto variationalSolver = simulation_setup::createVariationalEquationsSolver< double, double >(
            bodies, propagatorSettings, fixture.parametersToEstimate, true );

    fixture.stmInterface = std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< double > >(
            variationalSolver->getStateTransitionMatrixInterface( ) );
    fixture.simulator =
            std::dynamic_pointer_cast< MultiArcDynamicsSimulator< double, double > >( variationalSolver->getDynamicsSimulatorBase( ) );

    // These casts establish that the fixture exercises the production pure-multi-arc code path.
    BOOST_REQUIRE( fixture.stmInterface != nullptr );
    BOOST_REQUIRE( fixture.simulator != nullptr );

    fixture.fullStateTransitionSize = fixture.stmInterface->getFullStateTransitionMatrixSize( );
    fixture.fullSensitivitySize = fixture.stmInterface->getFullSensitivityMatrixSize( );

    return fixture;
}

//! Verify that exact endpoints use stored endpoint states while nearby interior epochs are genuinely
//! interpolated. The large absolute epoch specifically guards against an absolute-time tolerance snapping a
//! valid interior request to a boundary.
BOOST_AUTO_TEST_CASE( test_EvaluateArcStateAtTime_NearBoundaryEpochsAreInterpolated )
{
    const double arcInitialTime = 1.0E9;
    const double arcFinalTime = arcInitialTime + physical_constants::JULIAN_YEAR;

    // A linear six-component state has an exact Lagrange-interpolation result, so any error observed below is
    // caused by endpoint selection/windowing rather than approximation error in the reference function.
    auto stateAtTime = [ arcInitialTime ]( const double time ) {
        Eigen::VectorXd state( 6 );
        const double elapsedTime = time - arcInitialTime;
        for( int i = 0; i < state.size( ); ++i )
        {
            state( i ) = 10.0 * static_cast< double >( i + 1 ) + 1.0E-3 * static_cast< double >( i + 1 ) * elapsedTime;
        }
        return state;
    };

    // Samples near both ends provide enough support for the adaptive interpolation window without filling the
    // entire year-long interval.
    std::map< double, Eigen::VectorXd > arcSolution;
    for( int i = 0; i <= 8; ++i )
    {
        const double time = arcInitialTime + static_cast< double >( i );
        arcSolution[ time ] = stateAtTime( time );
    }
    for( int i = 8; i >= 0; --i )
    {
        const double time = arcFinalTime - static_cast< double >( i );
        arcSolution[ time ] = stateAtTime( time );
    }

    // Exact endpoint requests must bypass interpolation and return the corresponding stored samples.
    const Eigen::VectorXd initialState =
            simulation_setup::evaluateArcStateAtTime< double, double >( arcSolution, arcInitialTime, arcInitialTime, arcFinalTime, 0, 6 );
    const Eigen::VectorXd finalState =
            simulation_setup::evaluateArcStateAtTime< double, double >( arcSolution, arcFinalTime, arcInitialTime, arcFinalTime, 0, 6 );
    BOOST_CHECK_SMALL( ( initialState - stateAtTime( arcInitialTime ) ).norm( ), std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( ( finalState - stateAtTime( arcFinalTime ) ).norm( ), std::numeric_limits< double >::epsilon( ) );

    // A point just inside the initial boundary must reproduce the analytical state and must not be snapped back
    // to the initial sample merely because the absolute epoch is large.
    const double nearStartEpoch = arcInitialTime + 1.25;
    const Eigen::VectorXd nearStartState =
            simulation_setup::evaluateArcStateAtTime< double, double >( arcSolution, nearStartEpoch, arcInitialTime, arcFinalTime, 0, 6 );
    BOOST_CHECK_SMALL( ( nearStartState - stateAtTime( nearStartEpoch ) ).norm( ), 1.0E-6 );
    BOOST_CHECK_GT( ( nearStartState - stateAtTime( arcInitialTime ) ).norm( ), 1.0E-3 );

    // Repeat the same distinction at the final boundary to cover both branches of the endpoint handling.
    const double nearFinalEpoch = arcFinalTime - 1.25;
    const Eigen::VectorXd nearFinalState =
            simulation_setup::evaluateArcStateAtTime< double, double >( arcSolution, nearFinalEpoch, arcInitialTime, arcFinalTime, 0, 6 );
    BOOST_CHECK_SMALL( ( nearFinalState - stateAtTime( nearFinalEpoch ) ).norm( ), 1.0E-6 );
    BOOST_CHECK_GT( ( nearFinalState - stateAtTime( arcFinalTime ) ).norm( ), 1.0E-3 );
}

//! Verify explicit-arc STM selection and its input contract. At a shared boundary the accessor must distinguish
//! the propagated left arc from the initial right arc, pad each into the correct parameter block, and reject
//! invalid arc indices or epochs rather than silently falling back to time-based arc selection.
BOOST_AUTO_TEST_CASE( test_StmForArc_SelectionAndValidation )
{
    auto fixture = buildTwoArcFixture( );

    // The test addresses the second of exactly two arcs, so establish that fixture invariant before indexing it.
    BOOST_REQUIRE_EQUAL( fixture.arcStartTimes.size( ), 2u );
    const double connectionEpoch = fixture.arcStartTimes[ 1 ];

    // Verify the structural property for arc 1: identity-Phi, zero-S at the arc start.
    Eigen::MatrixXd fullArc1 = fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, connectionEpoch );

    // The padded matrix keeps one propagated-state row block and all estimated-parameter columns.
    BOOST_CHECK_EQUAL( fullArc1.rows( ), 6 );
    BOOST_CHECK_EQUAL( fullArc1.cols( ), fixture.fullStateTransitionSize + fixture.fullSensitivitySize );

    // Layout: column block [0, 6) holds arc 0's 6x6 STM slot, [6, 12) holds arc 1's. Sensitivity columns follow.
    const Eigen::Matrix< double, 6, 6 > arc0Block = fullArc1.block< 6, 6 >( 0, 0 );
    const Eigen::Matrix< double, 6, 6 > arc1Block = fullArc1.block< 6, 6 >( 0, 6 );
    const Eigen::Matrix< double, 6, 6 > identity6 = Eigen::Matrix< double, 6, 6 >::Identity( );
    // At its own initial epoch, arc 1 responds identically to its initial state and not at all to arc 0's state.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( arc1Block, identity6, 1.0E-12 );
    for( int row = 0; row < 6; ++row )
    {
        for( int col = 0; col < 6; ++col )
        {
            BOOST_CHECK_SMALL( std::fabs( arc0Block( row, col ) ), 1.0E-12 );
        }
    }

    // Sensitivity block (columns fullStateTransitionSize..end) must be exactly zero at the arc start.
    if( fixture.fullSensitivitySize > 0 )
    {
        Eigen::MatrixXd sensitivityBlock = fullArc1.block( 0, fixture.fullStateTransitionSize, 6, fixture.fullSensitivitySize );
        for( int row = 0; row < sensitivityBlock.rows( ); ++row )
        {
            for( int col = 0; col < sensitivityBlock.cols( ); ++col )
            {
                BOOST_CHECK_SMALL( std::fabs( sensitivityBlock( row, col ) ), 1.0E-12 );
            }
        }
    }

    // Contrast: the per-arc-index accessor for arc 0 at the same boundary epoch (which is arc 0's end time) returns the
    // propagated Phi in arc 0's 6-block (not identity, since arcDuration > 0) and zeros in arc 1's 6-block.
    // This is the matrix needed for the left-arc variational block in the inter-arc continuity assembly; the time-keyed lookup cannot
    // retrieve it because at the shared boundary the hunt scheme picks arc 1 as the "current" arc.
    Eigen::MatrixXd fullArc0 = fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 0, connectionEpoch );
    const Eigen::Matrix< double, 6, 6 > arc0BlockFromLeft = fullArc0.block< 6, 6 >( 0, 0 );
    const Eigen::Matrix< double, 6, 6 > arc1BlockFromLeft = fullArc0.block< 6, 6 >( 0, 6 );
    // Arc 0 has propagated for a full arc duration, proving the accessor selected its final variational state
    // rather than the identity belonging to arc 1 at the same numerical epoch.
    BOOST_CHECK_GT( ( arc0BlockFromLeft - identity6 ).norm( ), 1.0E-3 );
    for( int row = 0; row < 6; ++row )
    {
        for( int col = 0; col < 6; ++col )
        {
            BOOST_CHECK_SMALL( std::fabs( arc1BlockFromLeft( row, col ) ), 1.0E-12 );
        }
    }

    // Negative and upper-bound indices cover both invalid sides of the explicit arc selector.
    BOOST_CHECK_THROW( fixture.stmInterface->getCombinedStateTransitionAndSensitivityMatrixForArc( -1, fixture.arcStartTimes[ 0 ] ),
                       std::runtime_error );
    BOOST_CHECK_THROW( fixture.stmInterface->getCombinedStateTransitionAndSensitivityMatrixForArc( 2, fixture.arcStartTimes[ 0 ] ),
                       std::runtime_error );
    // An epoch after every arc is invalid even when the selected arc index itself is valid.
    BOOST_CHECK_THROW( fixture.stmInterface->getCombinedStateTransitionAndSensitivityMatrixForArc( 0, fixture.arcEndTimes[ 1 ] + 1.0 ),
                       std::runtime_error );
    // The full padded overload must enforce the same interval contract for an epoch before the selected arc.
    BOOST_CHECK_THROW(
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, fixture.arcStartTimes[ 0 ] - 1.0 ),
            std::runtime_error );
}

//! Exact assembly check: the design matrix is the right-arc variational block minus the left-arc variational block,
//! its columns are normalized before the normal equations are formed, the residual sign is negative state discrepancy,
//! dense rank-deficient PSD weights are handled without whitening/Cholesky assumptions, and accepted
//! round-off-level asymmetry is removed consistently from the rank, cost, gradient, and normal matrix.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_DensePsdExactNormalEquations )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    const double connectionEpoch = fixture.arcStartTimes[ 1 ];

    // An outer product exercises a dense PSD matrix whose constrained dimension is one. The skew perturbation is
    // deliberately below validation tolerance but large enough to reveal use of the raw matrix in any one term.
    Eigen::Matrix< double, 6, 1 > rankOneWeightVector;
    rankOneWeightVector << 2.0, -1.0, 0.5, 0.25, -0.75, 1.5;
    Eigen::Matrix< double, 6, 6 > denseRankOneWeightMatrix = rankOneWeightVector * rankOneWeightVector.transpose( );
    Eigen::Matrix< double, 6, 6 > nearSymmetricRankOneWeightMatrix = denseRankOneWeightMatrix;
    nearSymmetricRankOneWeightMatrix( 0, 1 ) += 1.0E-13;
    nearSymmetricRankOneWeightMatrix( 1, 0 ) -= 1.0E-13;
    const double constraintScalingFactor = 2.0;
    auto settings = simulation_setup::generalContinuity(
            { "Earth" },
            { { "Earth", { connectionEpoch } } },
            std::map< std::string, std::vector< Eigen::MatrixXd > >{ { "Earth", { nearSymmetricRankOneWeightMatrix } } },
            constraintScalingFactor );

    // Distinct factors per column make an omitted or reversed normalization operation observable.
    Eigen::VectorXd columnNormalizationFactors( totalParameterSize );
    for( int i = 0; i < totalParameterSize; ++i )
    {
        columnNormalizationFactors( i ) = 0.7 + 0.13 * static_cast< double >( i + 1 );
    }

    auto contribution = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                                  fixture.parametersToEstimate,
                                                                                  fixture.simulator,
                                                                                  fixture.stmInterface,
                                                                                  columnNormalizationFactors,
                                                                                  totalParameterSize );

    // The assembled quantities must occupy the estimator's full parameter space and retain one diagnostic for the
    // single configured boundary. This structural check is kept here because the exact algebra below exercises the
    // same output more strongly than a separate position-only smoke test would.
    BOOST_REQUIRE_EQUAL( contribution.additionalNormalMatrix.rows( ), totalParameterSize );
    BOOST_REQUIRE_EQUAL( contribution.additionalNormalMatrix.cols( ), totalParameterSize );
    BOOST_REQUIRE_EQUAL( contribution.additionalRightHandSide.size( ), totalParameterSize );
    BOOST_REQUIRE_EQUAL( contribution.perPairDiscrepancies.size( ), 1u );

    // Reconstruct the expected design matrix directly from the two explicit-arc accessors. This independently
    // checks the right-minus-left sign and the application of estimator column normalization.
    Eigen::MatrixXd leftArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 0, connectionEpoch );
    Eigen::MatrixXd rightArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, connectionEpoch );
    Eigen::MatrixXd continuityDesignMatrix = rightArcVariationalMatrix - leftArcVariationalMatrix;
    for( int col = 0; col < totalParameterSize; ++col )
    {
        continuityDesignMatrix.col( col ) /= columnNormalizationFactors( col );
    }

    // The symmetric rank-one matrix is the effective weight: the accepted skew part must disappear everywhere.
    const Eigen::VectorXd& stateDiscrepancy = contribution.perPairDiscrepancies.at( 0 );
    BOOST_CHECK_EQUAL( stateDiscrepancy.rows( ), 6 );
    BOOST_CHECK( std::isfinite( stateDiscrepancy.norm( ) ) );
    const Eigen::Matrix< double, 6, 6 > scaledConstraintWeight =
            denseRankOneWeightMatrix / constraintScalingFactor;  // rank(weight matrix)=1 -> total constrained dimension=1
    const Eigen::MatrixXd expectedNormalMatrixContribution =
            continuityDesignMatrix.transpose( ) * scaledConstraintWeight * continuityDesignMatrix;
    const Eigen::VectorXd expectedRightHandSideContribution =
            -continuityDesignMatrix.transpose( ) * ( scaledConstraintWeight * stateDiscrepancy );
    const double expectedCost = 0.5 * stateDiscrepancy.transpose( ) * scaledConstraintWeight * stateDiscrepancy;

    // Compare all three products of the quadratic model separately so a sign, factor-of-two, or symmetrization
    // error cannot be hidden by agreement in another output.
    BOOST_CHECK_LT( ( contribution.additionalNormalMatrix - expectedNormalMatrixContribution ).norm( ) /
                            std::max( expectedNormalMatrixContribution.norm( ), 1.0E-30 ),
                    1.0E-12 );
    BOOST_CHECK_LT( ( contribution.additionalRightHandSide - expectedRightHandSideContribution ).norm( ) /
                            std::max( expectedRightHandSideContribution.norm( ), 1.0E-30 ),
                    1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( contribution.totalConstraintCost, expectedCost, 1.0E-12 );

    // The exact matrix equality already checks symmetry algebraically; the eigenspectrum additionally confirms
    // that a dense rank-deficient PSD weight produces a finite positive-semidefinite information contribution.
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( contribution.additionalNormalMatrix );
    const double largestEigenvalue = std::max( solver.eigenvalues( ).maxCoeff( ), 1.0 );
    BOOST_CHECK_GE( solver.eigenvalues( ).minCoeff( ), -1.0E-9 * largestEigenvalue );
    BOOST_CHECK_GE( contribution.totalConstraintCost, 0.0 );
}

//! Orbit14-style mid-gap connection epoch. The adjacent arcs overlap around the connection epoch so both sides are
//! evaluated strictly inside their propagated intervals; neither STM block is the arc-start identity.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_MidGapConnectionEpoch )
{
    auto fixture = buildTwoArcFixture( {}, true );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    const double connectionEpoch = 0.5 * ( fixture.arcStartTimes[ 1 ] + fixture.arcEndTimes[ 0 ] );
    // Establish that the chosen epoch is strictly interior to both overlapping arcs; endpoint shortcuts would
    // otherwise make this an ineffective interpolation test.
    BOOST_CHECK_GT( connectionEpoch, fixture.arcStartTimes[ 1 ] );
    BOOST_CHECK_LT( connectionEpoch, fixture.arcEndTimes[ 0 ] );

    auto settings = simulation_setup::fullStateContinuity( { "Earth" }, { { "Earth", { connectionEpoch } } }, 1.0, 1.0, 1.0 );
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );
    auto contribution = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                                  fixture.parametersToEstimate,
                                                                                  fixture.simulator,
                                                                                  fixture.stmInterface,
                                                                                  columnNormalizationFactors,
                                                                                  totalParameterSize );

    // Interior evaluation must still produce one finite discrepancy and a valid non-negative quadratic cost.
    BOOST_REQUIRE_EQUAL( contribution.perPairDiscrepancies.size( ), 1u );
    BOOST_CHECK( std::isfinite( contribution.perPairDiscrepancies.at( 0 ).norm( ) ) );
    BOOST_CHECK_GE( contribution.totalConstraintCost, 0.0 );

    Eigen::MatrixXd leftArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 0, connectionEpoch );
    Eigen::MatrixXd rightArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, connectionEpoch );
    const Eigen::Matrix< double, 6, 6 > identity6 = Eigen::Matrix< double, 6, 6 >::Identity( );

    // Neither side may be an initial-state identity at the mid-gap epoch; both variational matrices must have
    // been interpolated from propagated solutions.
    BOOST_CHECK_GT( ( leftArcVariationalMatrix.block< 6, 6 >( 0, 0 ) - identity6 ).norm( ), 1.0E-3 );
    BOOST_CHECK_GT( ( rightArcVariationalMatrix.block< 6, 6 >( 0, 6 ) - identity6 ).norm( ), 1.0E-3 );
}

//! Finite-difference check of the continuity design matrix. This repropagates the two-arc problem with one
//! perturbed arc initial-state component at a time, so it is independent of the STM accessor used by assembly.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_FiniteDifferenceInitialStatePartials )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    const double connectionEpoch = fixture.arcStartTimes[ 1 ];

    // The perturbation loop below assumes two six-component arc-initial-state parameter blocks and no globals.
    BOOST_REQUIRE_EQUAL( totalParameterSize, 12 );

    Eigen::MatrixXd leftArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 0, connectionEpoch );
    Eigen::MatrixXd rightArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, connectionEpoch );
    const Eigen::MatrixXd analyticalContinuityDesignMatrix = rightArcVariationalMatrix - leftArcVariationalMatrix;

    // Build an entirely propagation-based central-difference Jacobian for comparison with the variational result.
    Eigen::MatrixXd finiteDifferenceContinuityDesignMatrix( 6, totalParameterSize );
    finiteDifferenceContinuityDesignMatrix.setZero( );

    for( int column = 0; column < totalParameterSize; ++column )
    {
        // Position and velocity use perturbations appropriate to their units, keeping truncation and numerical
        // propagation errors small relative to the derivative signal.
        const double perturbation = ( column % 6 < 3 ) ? 10.0 : 1.0E-3;
        std::vector< Eigen::Matrix< double, 6, 1 > > positivePerturbations( 2, Eigen::Matrix< double, 6, 1 >::Zero( ) );
        std::vector< Eigen::Matrix< double, 6, 1 > > negativePerturbations( 2, Eigen::Matrix< double, 6, 1 >::Zero( ) );
        positivePerturbations.at( static_cast< unsigned int >( column / 6 ) )( column % 6 ) = perturbation;
        negativePerturbations.at( static_cast< unsigned int >( column / 6 ) )( column % 6 ) = -perturbation;

        auto positiveFixture = buildTwoArcFixture( positivePerturbations );
        auto negativeFixture = buildTwoArcFixture( negativePerturbations );
        const Eigen::VectorXd positiveDiscrepancy = getSinglePairDiscrepancy( positiveFixture );
        const Eigen::VectorXd negativeDiscrepancy = getSinglePairDiscrepancy( negativeFixture );
        // Both perturbed runs must preserve the translational discrepancy layout before they are differenced.
        BOOST_CHECK_EQUAL( positiveDiscrepancy.rows( ), 6 );
        BOOST_CHECK_EQUAL( negativeDiscrepancy.rows( ), 6 );
        finiteDifferenceContinuityDesignMatrix.col( column ) = ( positiveDiscrepancy - negativeDiscrepancy ) / ( 2.0 * perturbation );
    }

    // Scale each comparison by the analytical column magnitude so weak and strong sensitivities receive a
    // meaningful common relative tolerance.
    for( int column = 0; column < totalParameterSize; ++column )
    {
        const double columnScale = std::max( analyticalContinuityDesignMatrix.col( column ).norm( ), 1.0 );
        BOOST_CHECK_LT( ( finiteDifferenceContinuityDesignMatrix.col( column ) - analyticalContinuityDesignMatrix.col( column ) ).norm( ) /
                                columnScale,
                        1.0E-4 );
    }
}

//! Component masks must affect only the requested discrepancy components, with constrained dimension equal to the
//! rank of the constraint weight matrix.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_ComponentMasksAndRelativeRank )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    const double connectionEpoch = fixture.arcStartTimes[ 1 ];
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );

    // Use unequal position and velocity weights so incorrect component selection cannot accidentally cancel.
    auto positionSettings = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { connectionEpoch } } }, 2.0, 5.0 );
    auto velocitySettings = simulation_setup::velocityOnlyContinuity( { "Earth" }, { { "Earth", { connectionEpoch } } }, 4.0, 5.0 );
    auto fullSettings = simulation_setup::fullStateContinuity( { "Earth" }, { { "Earth", { connectionEpoch } } }, 2.0, 4.0, 5.0 );

    auto positionContribution = assembleInterArcContinuityContribution< double, double >( { positionSettings },
                                                                                          fixture.parametersToEstimate,
                                                                                          fixture.simulator,
                                                                                          fixture.stmInterface,
                                                                                          columnNormalizationFactors,
                                                                                          totalParameterSize );
    auto velocityContribution = assembleInterArcContinuityContribution< double, double >( { velocitySettings },
                                                                                          fixture.parametersToEstimate,
                                                                                          fixture.simulator,
                                                                                          fixture.stmInterface,
                                                                                          columnNormalizationFactors,
                                                                                          totalParameterSize );
    auto fullContribution = assembleInterArcContinuityContribution< double, double >( { fullSettings },
                                                                                      fixture.parametersToEstimate,
                                                                                      fixture.simulator,
                                                                                      fixture.stmInterface,
                                                                                      columnNormalizationFactors,
                                                                                      totalParameterSize );

    // Compute each expected cost from the same discrepancy. The denominators encode ranks three, three, and six
    // respectively, directly testing constrained-dimension accounting for each mask.
    const Eigen::VectorXd& stateDiscrepancy = fullContribution.perPairDiscrepancies.at( 0 );
    BOOST_CHECK_EQUAL( stateDiscrepancy.rows( ), 6 );
    const double expectedPositionCost = 0.5 * ( 2.0 / ( 5.0 * 3.0 ) ) * stateDiscrepancy.head( 3 ).squaredNorm( );
    const double expectedVelocityCost = 0.5 * ( 4.0 / ( 5.0 * 3.0 ) ) * stateDiscrepancy.tail( 3 ).squaredNorm( );
    const double expectedFullCost = 0.5 * ( 1.0 / ( 5.0 * 6.0 ) ) *
            ( 2.0 * stateDiscrepancy.head( 3 ).squaredNorm( ) + 4.0 * stateDiscrepancy.tail( 3 ).squaredNorm( ) );

    // Each preset must penalize exactly its selected state components with its own component weights.
    BOOST_CHECK_CLOSE_FRACTION( positionContribution.totalConstraintCost, expectedPositionCost, 1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( velocityContribution.totalConstraintCost, expectedVelocityCost, 1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( fullContribution.totalConstraintCost, expectedFullCost, 1.0E-12 );

    // Very small positive weights must still be counted by relative rank; this used to fail with an absolute
    // 1e-12 floor in rankOf6x6PsdMatrix.
    auto tinySettings = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { connectionEpoch } } }, 1.0E-30, 1.0 );
    auto tinyContribution = assembleInterArcContinuityContribution< double, double >( { tinySettings },
                                                                                      fixture.parametersToEstimate,
                                                                                      fixture.simulator,
                                                                                      fixture.stmInterface,
                                                                                      columnNormalizationFactors,
                                                                                      totalParameterSize );
    // Scale-aware rank detection must retain the three active directions even when every eigenvalue is tiny.
    BOOST_CHECK_GT( tinyContribution.totalConstraintCost, 0.0 );
    BOOST_CHECK( std::isfinite( tinyContribution.totalConstraintCost ) );
}

//! Multi-body settings must use each body's own state/STM rows and connection epoch. The shared settings path is
//! compared with independent Earth and Mars contributions, while the Mars contribution is also reconstructed from
//! production layout metadata to catch accidental use of the first body's rows.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_MultiBodyRowsAndEpochs )
{
    auto fixture = buildTwoBodyTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );

    const double earthConnectionEpoch = fixture.arcStartTimes[ 1 ];
    const double marsConnectionEpoch = 0.5 * ( fixture.arcStartTimes[ 1 ] + fixture.arcEndTimes[ 0 ] );
    // The epochs must be genuinely distinct and both valid for their arcs; otherwise the body-to-epoch mapping
    // could be wrong without changing the result.
    BOOST_REQUIRE_NE( earthConnectionEpoch, marsConnectionEpoch );
    BOOST_CHECK_GT( marsConnectionEpoch, fixture.arcStartTimes[ 1 ] );
    BOOST_CHECK_LT( marsConnectionEpoch, fixture.arcEndTimes[ 0 ] );

    // Assemble both bodies once through a shared settings object, which is the API path under test.
    auto multiBodySettings = simulation_setup::positionOnlyContinuity(
            { "Earth", "Mars" },
            std::map< std::string, std::vector< double > >{ { "Earth", { earthConnectionEpoch } }, { "Mars", { marsConnectionEpoch } } },
            1.0,
            2.0 );
    auto multiBodyContribution = assembleInterArcContinuityContribution< double, double >( { multiBodySettings },
                                                                                           fixture.parametersToEstimate,
                                                                                           fixture.simulator,
                                                                                           fixture.stmInterface,
                                                                                           columnNormalizationFactors,
                                                                                           totalParameterSize );

    // Independent one-body settings form a reference with no shared map lookup or multi-body iteration. Keeping
    // their contributions separate also permits a direct inspection of the Mars row selection.
    auto earthSettings = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { earthConnectionEpoch } } }, 1.0, 2.0 );
    auto marsSettings = simulation_setup::positionOnlyContinuity( { "Mars" }, { { "Mars", { marsConnectionEpoch } } }, 1.0, 2.0 );
    auto earthContribution = assembleInterArcContinuityContribution< double, double >( { earthSettings },
                                                                                       fixture.parametersToEstimate,
                                                                                       fixture.simulator,
                                                                                       fixture.stmInterface,
                                                                                       columnNormalizationFactors,
                                                                                       totalParameterSize );
    auto marsContribution = assembleInterArcContinuityContribution< double, double >( { marsSettings },
                                                                                      fixture.parametersToEstimate,
                                                                                      fixture.simulator,
                                                                                      fixture.stmInterface,
                                                                                      columnNormalizationFactors,
                                                                                      totalParameterSize );

    // Read the actual multi-arc layout instead of assuming Mars follows Earth by six rows. The non-zero offsets make
    // this assertion sensitive to an implementation that always slices the first propagated body's variational rows.
    const auto layout = fixture.stmInterface->getArcWiseAndFullSolutionInitialStateIndices( );
    const int marsFullRowsLeft = layout.at( 0 ).at( "Mars" ).second.first.first;
    const int marsFullRowsRight = layout.at( 1 ).at( "Mars" ).second.first.first;
    BOOST_REQUIRE_NE( marsFullRowsLeft, 0 );
    BOOST_REQUIRE_NE( marsFullRowsRight, 0 );

    Eigen::MatrixXd leftArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 0, marsConnectionEpoch );
    Eigen::MatrixXd rightArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, marsConnectionEpoch );
    Eigen::MatrixXd expectedMarsDesignMatrix = rightArcVariationalMatrix.block( marsFullRowsRight, 0, 6, totalParameterSize ) -
            leftArcVariationalMatrix.block( marsFullRowsLeft, 0, 6, totalParameterSize );
    Eigen::Matrix< double, 6, 6 > scaledMarsWeight = Eigen::Matrix< double, 6, 6 >::Zero( );
    scaledMarsWeight.block< 3, 3 >( 0, 0 ) = Eigen::Matrix3d::Identity( ) / ( 2.0 * 3.0 );
    const Eigen::MatrixXd expectedMarsNormalMatrix = expectedMarsDesignMatrix.transpose( ) * scaledMarsWeight * expectedMarsDesignMatrix;

    // Agreement with the independently reconstructed Mars-only normal matrix proves that both explicit arc
    // variational matrices were sliced at the Mars rows.
    BOOST_CHECK_LT( ( marsContribution.additionalNormalMatrix - expectedMarsNormalMatrix ).norm( ) /
                            std::max( expectedMarsNormalMatrix.norm( ), 1.0E-30 ),
                    1.0E-12 );

    // Both paths must preserve one six-component discrepancy per body in body-list order.
    BOOST_REQUIRE_EQUAL( multiBodyContribution.perPairDiscrepancies.size( ), 2u );
    BOOST_CHECK_EQUAL( multiBodyContribution.perPairDiscrepancies.at( 0 ).rows( ), 6 );
    BOOST_CHECK_EQUAL( multiBodyContribution.perPairDiscrepancies.at( 1 ).rows( ), 6 );
    BOOST_CHECK_LT( ( multiBodyContribution.perPairDiscrepancies.at( 0 ) - earthContribution.perPairDiscrepancies.at( 0 ) ).norm( ),
                    1.0E-12 );
    BOOST_CHECK_LT( ( multiBodyContribution.perPairDiscrepancies.at( 1 ) - marsContribution.perPairDiscrepancies.at( 0 ) ).norm( ),
                    1.0E-12 );
    // Each one-body reference is normalized by rank three, whereas the shared setting is normalized by the global
    // rank six. The factor one half therefore converts their sum to the expected multi-body quadratic contribution.
    const Eigen::MatrixXd independentNormalMatrix =
            0.5 * ( earthContribution.additionalNormalMatrix + marsContribution.additionalNormalMatrix );
    const Eigen::VectorXd independentRightHandSide =
            0.5 * ( earthContribution.additionalRightHandSide + marsContribution.additionalRightHandSide );
    BOOST_CHECK_LT( ( multiBodyContribution.additionalNormalMatrix - independentNormalMatrix ).norm( ) /
                            std::max( independentNormalMatrix.norm( ), 1.0E-30 ),
                    1.0E-12 );
    BOOST_CHECK_LT( ( multiBodyContribution.additionalRightHandSide - independentRightHandSide ).norm( ) /
                            std::max( independentRightHandSide.norm( ), 1.0E-30 ),
                    1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( multiBodyContribution.totalConstraintCost,
                                0.5 * ( earthContribution.totalConstraintCost + marsContribution.totalConstraintCost ),
                                1.0E-12 );
}

//! Verify both estimator-coordinate scaling and global constrained-dimension scaling. Parameter normalization may
//! change the coordinates of the normal equations but not their physical update, while duplicating a constraint
//! term doubles the global rank and number of terms without changing their accumulated information.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_NormalisationAndGlobalDimension )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;

    auto settings = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { fixture.arcStartTimes[ 1 ] } } }, 1.0, 1.0 );

    // A uniform change of parameter normalization isolates the coordinate transformation from any change in
    // relative column scales.
    Eigen::VectorXd unitColumnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );
    Eigen::VectorXd scaledColumnNormalizationFactors = Eigen::VectorXd::Constant( totalParameterSize, 3.5 );

    auto unit = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                          fixture.parametersToEstimate,
                                                                          fixture.simulator,
                                                                          fixture.stmInterface,
                                                                          unitColumnNormalizationFactors,
                                                                          totalParameterSize );
    auto scaled = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                            fixture.parametersToEstimate,
                                                                            fixture.simulator,
                                                                            fixture.stmInterface,
                                                                            scaledColumnNormalizationFactors,
                                                                            totalParameterSize );
    // A separate assembly changes only the number of observations, testing the Orbit14 relative scaling
    // independently of parameter normalization.
    const int numberOfObservations = 7;
    auto observationScaled = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                                       fixture.parametersToEstimate,
                                                                                       fixture.simulator,
                                                                                       fixture.stmInterface,
                                                                                       unitColumnNormalizationFactors,
                                                                                       totalParameterSize,
                                                                                       numberOfObservations );

    // Parameter normalization must not change the physical cost. Observation count, by contrast, must scale the
    // normal matrix, gradient term, and cost linearly by the same factor.
    BOOST_CHECK_CLOSE_FRACTION( unit.totalConstraintCost, scaled.totalConstraintCost, 1.0E-12 );
    BOOST_CHECK_LT( ( observationScaled.additionalNormalMatrix - numberOfObservations * unit.additionalNormalMatrix ).norm( ) /
                            std::max( observationScaled.additionalNormalMatrix.norm( ), 1.0E-30 ),
                    1.0E-12 );
    BOOST_CHECK_LT( ( observationScaled.additionalRightHandSide - numberOfObservations * unit.additionalRightHandSide ).norm( ) /
                            std::max( observationScaled.additionalRightHandSide.norm( ), 1.0E-30 ),
                    1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( observationScaled.totalConstraintCost, numberOfObservations * unit.totalConstraintCost, 1.0E-12 );

    // Uniformly scaling every column-normalization factor by k scales the normal matrix contribution by 1/k^2
    // and the right-hand-side contribution by 1/k. The normalized parameter update must therefore be divided
    // by k to recover the same physical update as the unit-normalization case.
    const double tinyPriorScale = 1.0E-3;
    Eigen::MatrixXd priorNormalMatrix = tinyPriorScale * Eigen::MatrixXd::Identity( totalParameterSize, totalParameterSize );

    Eigen::VectorXd unitParameterUpdate = ( unit.additionalNormalMatrix + priorNormalMatrix ).ldlt( ).solve( unit.additionalRightHandSide );
    Eigen::VectorXd scaledParameterUpdate =
            ( scaled.additionalNormalMatrix + priorNormalMatrix / ( 3.5 * 3.5 ) ).ldlt( ).solve( scaled.additionalRightHandSide );
    Eigen::VectorXd scaledPhysicalParameterUpdate = scaledParameterUpdate / 3.5;
    BOOST_CHECK_LT( ( unitParameterUpdate - scaledPhysicalParameterUpdate ).norm( ) / std::max( unitParameterUpdate.norm( ), 1.0E-30 ),
                    1.0E-5 );

    // Adding the same rank-three pair twice doubles the global constrained dimension as well as the number of
    // terms. Consequently, each term receives half the weight and the aggregate quadratic model is unchanged.
    auto duplicateSettings =
            simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { fixture.arcStartTimes[ 1 ] } } }, 1.0, 1.0 );
    auto duplicated = assembleInterArcContinuityContribution< double, double >( { settings, duplicateSettings },
                                                                                fixture.parametersToEstimate,
                                                                                fixture.simulator,
                                                                                fixture.stmInterface,
                                                                                unitColumnNormalizationFactors,
                                                                                totalParameterSize );
    BOOST_CHECK_LT( ( duplicated.additionalNormalMatrix - unit.additionalNormalMatrix ).norm( ),
                    1.0E-9 * std::max( unit.additionalNormalMatrix.norm( ), 1.0 ) );
    BOOST_CHECK_CLOSE_FRACTION( duplicated.totalConstraintCost, unit.totalConstraintCost, 1.0E-9 );
    // Diagnostics intentionally preserve both configured pairs even though their numerical values are identical.
    BOOST_CHECK_EQUAL( duplicated.perPairDiscrepancies.size( ), 2u );
}

//! Exercise continuity regularization through the complete iterative orbit-determination manager. The test first
//! establishes the constraint-free no-op behavior, then verifies per-iteration diagnostics, selected-iteration
//! diagnostics, reduction of the constrained jump, and preservation of an acceptable observation fit.
BOOST_AUTO_TEST_CASE( test_OdLoop_WithInterArcContinuity_EndToEnd )
{
    // A one-way Earth-to-Mars range link provides observations sensitive to both arc initial states without
    // introducing extra estimated dynamical parameters.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames = { "Earth", "Mars", "Sun" };
    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 4.0E7;
    const double buffer = 3.6E5;

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    SystemOfBodies bodies = createSystemOfBodies< double, double >( bodySettings );

    std::pair< std::string, std::string > marsStation( "Mars", "MarsStation" );
    createGroundStation( bodies.at( "Mars" ),
                         "MarsStation",
                         ( Eigen::Vector3d( ) << 100.0, 0.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );

    // Keep the propagated dynamics identical to the smaller fixtures so this test focuses on manager integration.
    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Earth" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    std::vector< std::string > bodiesToIntegrate = { "Earth" };
    std::vector< std::string > centralBodies = { "SSB" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    const double arcDuration = 1.0E7;
    const std::vector< double > arcStartTimes = { initialEphemerisTime + 1.0E5, initialEphemerisTime + 1.0E5 + arcDuration };
    const std::vector< double > arcEndTimes = { arcStartTimes[ 0 ] + arcDuration, arcStartTimes[ 1 ] + arcDuration };

    auto integratorSettings = rungeKutta4Settings< double >( 600.0 );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > propagatorSettingsList;
    for( unsigned int i = 0; i < arcStartTimes.size( ); ++i )
    {
        Eigen::VectorXd initialState = getInitialStateOfBody< double, double >( "Earth", "SSB", bodies, arcStartTimes[ i ] );
        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                centralBodies,
                accelerationModelMap,
                bodiesToIntegrate,
                initialState,
                arcStartTimes[ i ],
                integratorSettings,
                propagationTimeTerminationSettings( arcEndTimes[ i ] ) ) );
    }
    auto propagatorSettings = std::make_shared< MultiArcPropagatorSettings< double, double > >( propagatorSettingsList, false );

    // Estimate both arc initial states; this is the pure multi-arc parameter layout supported by the feature.
    auto parameterNames = getInitialMultiArcParameterSettings< double, double >( propagatorSettings, bodies, arcStartTimes );
    auto parametersToEstimate = createParametersToEstimate< double >( parameterNames, bodies );

    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth", "" );
    linkEnds[ receiver ] = marsStation;

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList = { std::make_shared< ObservationModelSettings >(
            one_way_range, linkEnds ) };

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettingsList, propagatorSettings );

    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );

    // Keep observations away from arc endpoints and distribute them uniformly across both arcs, giving the
    // baseline estimator enough information to converge before regularization is introduced.
    const int observationsPerArc = 200;
    std::vector< double > observationTimes;
    for( unsigned int arc = 0; arc < arcStartTimes.size( ); ++arc )
    {
        const double dt = ( arcEndTimes[ arc ] - arcStartTimes[ arc ] - 2.0 * 12000.0 ) / static_cast< double >( observationsPerArc - 1 );
        double t = arcStartTimes[ arc ] + 12000.0;
        for( int i = 0; i < observationsPerArc; ++i )
        {
            observationTimes.push_back( t );
            t += dt;
        }
    }
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementInput = {
        std::make_shared< TabulatedObservationSimulationSettings< double > >( one_way_range, linkEnds, observationTimes, receiver )
    };
    auto observations =
            simulateObservations< double, double >( measurementInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Apply unit-aware perturbations to every arc so iterative correction and best-iteration selection are
    // exercised rather than returning the simulation truth immediately.
    Eigen::VectorXd initialEstimate = truthParameters;
    for( unsigned int arc = 0; arc < arcStartTimes.size( ); ++arc )
    {
        initialEstimate( 6 * arc + 0 ) += 1.0;
        initialEstimate( 6 * arc + 1 ) += 1.0;
        initialEstimate( 6 * arc + 2 ) += 1.0;
        initialEstimate( 6 * arc + 3 ) += 1.0E-5;
        initialEstimate( 6 * arc + 4 ) += 1.0E-5;
        initialEstimate( 6 * arc + 5 ) += 1.0E-5;
    }
    parametersToEstimate->resetParameterValues( initialEstimate );

    // The baseline run establishes backwards compatibility: no configured priors means no continuity work and no
    // continuity diagnostics in the output.
    auto estimationInputNoConstraint = std::make_shared< EstimationInput< double, double > >( observations );
    auto outputNoConstraint = orbitDeterminationManager.estimateParameters( estimationInputNoConstraint );
    BOOST_REQUIRE( outputNoConstraint != nullptr );
    BOOST_CHECK( outputNoConstraint->getInterArcContinuityCostHistory( ).empty( ) );
    BOOST_CHECK( outputNoConstraint->getInterArcContinuityDiscrepancyHistory( ).empty( ) );

    // Reset to the identical initial estimate, then enable a weak position prior. The attached constraint is the
    // only behavioral difference between the two manager runs.
    parametersToEstimate->resetParameterValues( initialEstimate );
    auto estimationInputWithConstraint = std::make_shared< EstimationInput< double, double > >( observations );
    auto constraint = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { arcStartTimes[ 1 ] } } }, 1.0, 1.0E14 );
    estimationInputWithConstraint->setInterArcContinuityConstraints( { constraint } );
    auto outputWithConstraint = orbitDeterminationManager.estimateParameters( estimationInputWithConstraint );

    // The manager must record one cost and one discrepancy collection for every iteration that evaluated the
    // constraint, with exactly one boundary discrepancy per collection.
    BOOST_REQUIRE( outputWithConstraint != nullptr );
    BOOST_REQUIRE_GT( outputWithConstraint->getInterArcContinuityCostHistory( ).size( ), 0u );
    BOOST_CHECK_EQUAL( outputWithConstraint->getInterArcContinuityCostHistory( ).size( ),
                       outputWithConstraint->getInterArcContinuityDiscrepancyHistory( ).size( ) );
    for( const auto& perIter : outputWithConstraint->getInterArcContinuityDiscrepancyHistory( ) )
    {
        BOOST_CHECK_EQUAL( perIter.size( ), 1u );
    }
    for( double cost : outputWithConstraint->getInterArcContinuityCostHistory( ) )
    {
        BOOST_CHECK( std::isfinite( cost ) );
        BOOST_CHECK_GE( cost, 0.0 );
    }

    // The scalar diagnostic exposed on the final output must come from the same iteration selected by the combined
    // observation-plus-continuity objective.
    const auto& history = outputWithConstraint->getInterArcContinuityDiscrepancyHistory( );
    const int bestIterationIndex =
            std::max( 0, std::min( outputWithConstraint->bestIteration_, static_cast< int >( history.size( ) ) - 1 ) );
    BOOST_CHECK_EQUAL( outputWithConstraint->getInterArcContinuityCost( ),
                       outputWithConstraint->getInterArcContinuityCostHistory( ).at( bestIterationIndex ) );

    const auto& bestIterationDiscrepancies = history.at( bestIterationIndex );
    const auto& outputDiscrepancies = outputWithConstraint->getInterArcContinuityDiscrepancies( );
    // Selected output vectors must be exact copies of that history entry, not values from the final iteration.
    BOOST_REQUIRE_EQUAL( outputDiscrepancies.size( ), bestIterationDiscrepancies.size( ) );
    for( unsigned int i = 0; i < outputDiscrepancies.size( ); ++i )
    {
        BOOST_CHECK_SMALL( ( outputDiscrepancies.at( i ) - bestIterationDiscrepancies.at( i ) ).norm( ), 1.0E-14 );
    }

    // The prior must have a useful, non-catastrophic effect: the selected iteration reduces the constrained
    // position jump from its initial value while retaining a sub-metre observation residual RMS.
    const Eigen::VectorXd& bestIterationStateDiscrepancy = history.at( bestIterationIndex ).at( 0 );
    BOOST_CHECK_EQUAL( bestIterationStateDiscrepancy.rows( ), 6 );
    BOOST_CHECK_LT( bestIterationStateDiscrepancy.head( 3 ).norm( ), history.front( ).at( 0 ).head( 3 ).norm( ) );
    BOOST_CHECK_LT( linear_algebra::getVectorEntryRootMeanSquare( outputWithConstraint->residuals_ ), 1.0 );
    BOOST_CHECK( std::isfinite( bestIterationStateDiscrepancy.norm( ) ) );
}

//! Invalid configurations must fail with actionable diagnostics. This test covers both an epoch outside the
//! selected arcs and an ambiguous parameter set containing the same body's multi-arc initial state twice.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_InvalidConfigurationDiagnostics )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );

    // Place the epoch beyond both arcs so whichever side is validated first must reject it.
    const double outOfRangeConnectionEpoch = fixture.arcEndTimes[ 1 ] + 1.0E6;
    auto badSettings = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { outOfRangeConnectionEpoch } } }, 1.0, 1.0 );

    try
    {
        assembleInterArcContinuityContribution< double, double >( { badSettings },
                                                                  fixture.parametersToEstimate,
                                                                  fixture.simulator,
                                                                  fixture.stmInterface,
                                                                  columnNormalizationFactors,
                                                                  totalParameterSize );
        // Reaching this statement would mean the assembler attempted extrapolation or silently selected an arc.
        BOOST_FAIL( "Expected runtime_error from out-of-range connection epoch was not thrown." );
    }
    catch( const std::runtime_error& error )
    {
        // Check semantic parts of the diagnostic rather than its formatting, ensuring users receive enough
        // context to correct the configuration while allowing harmless wording changes.
        const std::string what = error.what( );
        BOOST_CHECK( what.find( "Inter-arc continuity connection epoch" ) != std::string::npos );
        BOOST_CHECK( what.find( "Earth" ) != std::string::npos );
        BOOST_CHECK( what.find( "outside the propagated interval" ) != std::string::npos );
        BOOST_CHECK( what.find( "Extend the arc propagation interval" ) != std::string::npos );
    }

    // Next, begin from the one valid Earth parameter and deliberately insert the same identifier twice. Sharing
    // the pointer is sufficient because the ambiguity is in parameter-set lookup, not parameter values.
    auto multiArcStateParameters = fixture.parametersToEstimate->getEstimatedMultiArcInitialStateParameters( );
    BOOST_REQUIRE_EQUAL( multiArcStateParameters.size( ), 1u );

    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< double, Eigen::Dynamic, 1 > > > >
            duplicatedInitialStateParameters = { multiArcStateParameters.at( 0 ), multiArcStateParameters.at( 0 ) };
    auto duplicatedParametersToEstimate = std::make_shared< EstimatableParameterSet< double > >(
            std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
            std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ),
            duplicatedInitialStateParameters );

    auto validSettings = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { fixture.arcStartTimes[ 1 ] } } }, 1.0, 1.0 );
    try
    {
        assembleInterArcContinuityContribution< double, double >( { validSettings },
                                                                  duplicatedParametersToEstimate,
                                                                  fixture.simulator,
                                                                  fixture.stmInterface,
                                                                  columnNormalizationFactors,
                                                                  totalParameterSize );
        // Silent success would make the selected parameter depend on container ordering.
        BOOST_FAIL( "Expected runtime_error from duplicate body initial-state parameters was not thrown." );
    }
    catch( const std::runtime_error& error )
    {
        // The diagnostic must identify both the affected body and the exact multiplicity/type requirement.
        const std::string what = error.what( );
        BOOST_CHECK( what.find( "Earth" ) != std::string::npos );
        BOOST_CHECK( what.find( "matches 2 arc-wise initial state parameters" ) != std::string::npos );
        BOOST_CHECK( what.find( "translational" ) != std::string::npos );
    }
}

//! Compare otherwise identical covariance analyses with and without an attached continuity prior. The added
//! information must be positive semidefinite and non-zero, and the covariance output must expose the discrepancy
//! and cost evaluated at its linearization point.
BOOST_AUTO_TEST_CASE( test_CovarianceAnalysis_WithInterArcContinuity_Tightens )
{
    // Recreate the production manager setup here because this test exercises CovarianceAnalysisInput/Output, not
    // only the lower-level assembler used by the reusable fixtures.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames = { "Earth", "Mars", "Sun" };
    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 4.0E7;
    const double buffer = 3.6E5;

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    SystemOfBodies bodies = createSystemOfBodies< double, double >( bodySettings );

    std::pair< std::string, std::string > marsStation( "Mars", "MarsStation" );
    createGroundStation( bodies.at( "Mars" ),
                         "MarsStation",
                         ( Eigen::Vector3d( ) << 100.0, 0.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Earth" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    std::vector< std::string > bodiesToIntegrate = { "Earth" };
    std::vector< std::string > centralBodies = { "SSB" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    const double arcDuration = 1.0E7;
    const std::vector< double > arcStartTimes = { initialEphemerisTime + 1.0E5, initialEphemerisTime + 1.0E5 + arcDuration };
    const std::vector< double > arcEndTimes = { arcStartTimes[ 0 ] + arcDuration, arcStartTimes[ 1 ] + arcDuration };

    auto integratorSettings = rungeKutta4Settings< double >( 600.0 );
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > propagatorSettingsList;
    for( unsigned int i = 0; i < arcStartTimes.size( ); ++i )
    {
        Eigen::VectorXd initialState = getInitialStateOfBody< double, double >( "Earth", "SSB", bodies, arcStartTimes[ i ] );
        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                centralBodies,
                accelerationModelMap,
                bodiesToIntegrate,
                initialState,
                arcStartTimes[ i ],
                integratorSettings,
                propagationTimeTerminationSettings( arcEndTimes[ i ] ) ) );
    }
    auto propagatorSettings = std::make_shared< MultiArcPropagatorSettings< double, double > >( propagatorSettingsList, false );

    auto parameterNames = getInitialMultiArcParameterSettings< double, double >( propagatorSettings, bodies, arcStartTimes );
    auto parametersToEstimate = createParametersToEstimate< double >( parameterNames, bodies );

    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth", "" );
    linkEnds[ receiver ] = marsStation;
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList = { std::make_shared< ObservationModelSettings >(
            one_way_range, linkEnds ) };

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettingsList, propagatorSettings );

    // Uniform observation coverage in both arcs makes the baseline normal matrix well formed, so any information
    // increase can be attributed to the attached prior.
    const int observationsPerArc = 100;
    std::vector< double > observationTimes;
    for( unsigned int arc = 0; arc < arcStartTimes.size( ); ++arc )
    {
        const double dt = ( arcEndTimes[ arc ] - arcStartTimes[ arc ] - 2.0 * 12000.0 ) / static_cast< double >( observationsPerArc - 1 );
        double t = arcStartTimes[ arc ] + 12000.0;
        for( int i = 0; i < observationsPerArc; ++i )
        {
            observationTimes.push_back( t );
            t += dt;
        }
    }
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementInput = {
        std::make_shared< TabulatedObservationSimulationSettings< double > >( one_way_range, linkEnds, observationTimes, receiver )
    };
    auto observations =
            simulateObservations< double, double >( measurementInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Compute reference and constrained covariances from the same propagated model and observation collection.
    auto unconstrainedInput = std::make_shared< CovarianceAnalysisInput< double, double > >( observations );
    auto unconstrainedOutput = orbitDeterminationManager.computeCovariance( unconstrainedInput );

    auto constrainedInput = std::make_shared< CovarianceAnalysisInput< double, double > >( observations );
    auto constraint = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { arcStartTimes[ 1 ] } } }, 1.0, 1.0E-15 );
    constrainedInput->setInterArcContinuityConstraints( { constraint } );
    auto constrainedOutput = orbitDeterminationManager.computeCovariance( constrainedInput );

    // Both analyses must complete before their information matrices and diagnostics can be compared.
    BOOST_REQUIRE( unconstrainedOutput != nullptr );
    BOOST_REQUIRE( constrainedOutput != nullptr );

    // Compare normalized inverse covariance (= normal matrix). Constrained must >= unconstrained in the PSD
    // sense: the difference (constrained - unconstrained) is PSD with smallest eigenvalue >= 0.
    Eigen::MatrixXd normalMatrixDifference =
            constrainedOutput->inverseNormalizedCovarianceMatrix_ - unconstrainedOutput->inverseNormalizedCovarianceMatrix_;
    Eigen::MatrixXd symmetricNormalMatrixDifference = 0.5 * ( normalMatrixDifference + normalMatrixDifference.transpose( ) );
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( symmetricNormalMatrixDifference );
    const double largestEigenvalue = std::max( solver.eigenvalues( ).maxCoeff( ), 1.0 );
    BOOST_CHECK_GE( solver.eigenvalues( ).minCoeff( ), -1.0E-9 * largestEigenvalue );
    // A strictly positive eigenvalue rules out an implementation that accepts the settings but adds a zero matrix.
    BOOST_CHECK_GT( solver.eigenvalues( ).maxCoeff( ), 0.0 );
    // Covariance-specific diagnostics must be populated at the same linearization point and remain finite.
    BOOST_CHECK_GT( constrainedOutput->getInterArcContinuityCost( ), 0.0 );
    BOOST_REQUIRE_EQUAL( constrainedOutput->getInterArcContinuityDiscrepancies( ).size( ), 1u );
    BOOST_CHECK( std::isfinite( constrainedOutput->getInterArcContinuityDiscrepancies( ).at( 0 ).norm( ) ) );
}

//! Verify the complete least-squares extension in one shared dense problem: empty soft additions preserve the old
//! call exactly, non-empty additions match equivalent a-priori information, malformed dimensions are rejected,
//! and soft additions compose correctly with pre-existing hard equality constraints.
BOOST_AUTO_TEST_CASE( test_LeastSquaresSoftAdditionsCompatibilityAndComposition )
{
    // Use a small dense, non-symmetric design problem so every parameter and normal-matrix entry participates.
    Eigen::MatrixXd designMatrix( 5, 3 );
    designMatrix << 1.0, 0.5, -0.2, 0.3, 1.2, 0.1, -0.7, 0.4, 1.0, 0.2, -0.3, 0.9, 1.1, 0.8, -0.5;
    Eigen::VectorXd residuals( 5 );
    residuals << 0.1, -0.2, 0.05, -0.05, 0.15;
    Eigen::VectorXd weights = Eigen::VectorXd::Constant( 5, 1.0 );
    Eigen::MatrixXd inverseApriori = 0.01 * Eigen::MatrixXd::Identity( 3, 3 );

    // Compare the pre-existing call signature with the extended signature whose optional matrices are empty.
    auto baseline =
            tudat::linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrix, residuals, weights, inverseApriori );
    auto withEmptyAdditions = tudat::linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                                                                    residuals,
                                                                                                    weights,
                                                                                                    inverseApriori,
                                                                                                    1.0E8,
                                                                                                    Eigen::MatrixXd( 0, 0 ),
                                                                                                    Eigen::VectorXd( 0 ),
                                                                                                    Eigen::MatrixXd( 0, 0 ),
                                                                                                    Eigen::VectorXd( 0 ),
                                                                                                    Eigen::MatrixXd( 0, 0 ),
                                                                                                    Eigen::VectorXd( 0 ) );

    // Shape equality is required before checking bit-level equality of the correction and complete augmented
    // normal matrix.
    BOOST_REQUIRE_EQUAL( baseline.first.size( ), withEmptyAdditions.first.size( ) );
    BOOST_REQUIRE_EQUAL( baseline.second.rows( ), withEmptyAdditions.second.rows( ) );
    BOOST_REQUIRE_EQUAL( baseline.second.cols( ), withEmptyAdditions.second.cols( ) );
    for( int i = 0; i < baseline.first.size( ); ++i )
    {
        BOOST_CHECK_EQUAL( baseline.first( i ), withEmptyAdditions.first( i ) );
    }
    for( int row = 0; row < baseline.second.rows( ); ++row )
    {
        for( int col = 0; col < baseline.second.cols( ); ++col )
        {
            BOOST_CHECK_EQUAL( baseline.second( row, col ), withEmptyAdditions.second( row, col ) );
        }
    }

    // Round-trip cross-check: supplying a non-empty addition should change the result in a known way.
    // Adding alpha*I to the normal matrix is equivalent to scaling the a-priori covariance inverse.
    const double alpha = 0.5;
    Eigen::MatrixXd alphaIdentity = alpha * Eigen::MatrixXd::Identity( 3, 3 );
    auto withAddition = tudat::linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                                                              residuals,
                                                                                              weights,
                                                                                              inverseApriori,
                                                                                              1.0E8,
                                                                                              Eigen::MatrixXd( 0, 0 ),
                                                                                              Eigen::VectorXd( 0 ),
                                                                                              Eigen::MatrixXd( 0, 0 ),
                                                                                              Eigen::VectorXd( 0 ),
                                                                                              alphaIdentity,
                                                                                              Eigen::VectorXd::Zero( 3 ) );
    auto withScaledPrior = tudat::linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix(
            designMatrix, residuals, weights, inverseApriori + alphaIdentity );
    for( int row = 0; row < withAddition.second.rows( ); ++row )
    {
        for( int col = 0; col < withAddition.second.cols( ); ++col )
        {
            BOOST_CHECK_CLOSE_FRACTION( withAddition.second( row, col ), withScaledPrior.second( row, col ), 1.0E-12 );
        }
    }

    // A non-empty soft normal matrix with the wrong parameter dimension must be rejected.
    BOOST_CHECK_THROW( tudat::linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                                                             residuals,
                                                                                             weights,
                                                                                             inverseApriori,
                                                                                             1.0E8,
                                                                                             Eigen::MatrixXd( 0, 0 ),
                                                                                             Eigen::VectorXd( 0 ),
                                                                                             Eigen::MatrixXd( 0, 0 ),
                                                                                             Eigen::VectorXd( 0 ),
                                                                                             Eigen::MatrixXd::Identity( 4, 4 ),
                                                                                             Eigen::VectorXd( 0 ) ),
                       std::runtime_error );
    // The soft right-hand side is validated independently, including when its paired matrix is empty.
    BOOST_CHECK_THROW( tudat::linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                                                             residuals,
                                                                                             weights,
                                                                                             inverseApriori,
                                                                                             1.0E8,
                                                                                             Eigen::MatrixXd( 0, 0 ),
                                                                                             Eigen::VectorXd( 0 ),
                                                                                             Eigen::MatrixXd( 0, 0 ),
                                                                                             Eigen::VectorXd( 0 ),
                                                                                             Eigen::MatrixXd( 0, 0 ),
                                                                                             Eigen::VectorXd::Ones( 4 ) ),
                       std::runtime_error );

    // Add a hard constraint x[0] = 0.5 and a soft quadratic term to the same dense problem. Combining these checks
    // here avoids repeating the common design matrix, residual, weights, and a-priori setup in a separate test.
    Eigen::MatrixXd constraintMultiplier( 1, 3 );
    constraintMultiplier << 1.0, 0.0, 0.0;
    Eigen::VectorXd constraintRhs = Eigen::VectorXd::Constant( 1, 0.5 );
    Eigen::MatrixXd soft = 0.2 * Eigen::MatrixXd::Identity( 3, 3 );
    Eigen::VectorXd softRhs( 3 );
    softRhs << 0.1, -0.05, 0.07;

    auto withBoth = tudat::linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                                                          residuals,
                                                                                          weights,
                                                                                          inverseApriori,
                                                                                          1.0E8,
                                                                                          constraintMultiplier,
                                                                                          constraintRhs,
                                                                                          Eigen::MatrixXd( 0, 0 ),
                                                                                          Eigen::VectorXd( 0 ),
                                                                                          soft,
                                                                                          softRhs );
    auto withHardOnly = tudat::linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix(
            designMatrix, residuals, weights, inverseApriori, 1.0E8, constraintMultiplier, constraintRhs );

    // One Lagrange multiplier extends the three-parameter solution to four entries, while the parameter correction
    // must satisfy the hard equality independently of the simultaneous soft term.
    BOOST_REQUIRE_EQUAL( withBoth.first.size( ), 4 );
    const Eigen::VectorXd parameterCorrection = withBoth.first.head( 3 );
    BOOST_CHECK_CLOSE_FRACTION( ( constraintMultiplier * parameterCorrection )( 0 ), constraintRhs( 0 ), 1.0E-10 );
    // The soft term must visibly affect the parameter solution and augment exactly the parameter block of the
    // hard-constrained normal matrix, leaving the multiplier rows and columns untouched.
    BOOST_CHECK_GT( ( withBoth.first.head( 3 ) - withHardOnly.first.head( 3 ) ).norm( ), 1.0E-6 );
    const Eigen::MatrixXd normalDifference = withBoth.second.topLeftCorner( 3, 3 ) - withHardOnly.second.topLeftCorner( 3, 3 );
    BOOST_CHECK_LT( ( normalDifference - soft ).norm( ), 1.0E-10 );
}

//! CovarianceAnalysisOutput's pre-existing final positional argument is the propagation-exception flag. New
//! continuity diagnostics must remain appended after it to avoid silently converting that bool to a cost.
BOOST_AUTO_TEST_CASE( test_CovarianceOutputConstructorRetainsPositionalCompatibility )
{
    // Supply the legacy positional arguments only, with the final boolean set true. Defaulted continuity fields
    // must be initialized without consuming or reinterpreting that boolean.
    auto output = CovarianceAnalysisOutput< double, double >( Eigen::MatrixXd::Identity( 1, 1 ),
                                                              Eigen::VectorXd::Ones( 1 ),
                                                              Eigen::VectorXd::Ones( 1 ),
                                                              Eigen::MatrixXd::Identity( 1, 1 ),
                                                              Eigen::MatrixXd::Zero( 0, 0 ),
                                                              Eigen::VectorXd::Zero( 0 ),
                                                              Eigen::MatrixXd::Zero( 0, 0 ),
                                                              Eigen::MatrixXd::Zero( 0, 0 ),
                                                              true );

    // Preserve the legacy flag and verify both newly appended diagnostics receive their no-constraint defaults.
    BOOST_CHECK( output.exceptionDuringPropagation_ );
    BOOST_CHECK_EQUAL( output.getInterArcContinuityCost( ), 0.0 );
    BOOST_CHECK( output.getInterArcContinuityDiscrepancies( ).empty( ) );
}

//! Verify the matrices produced by each settings factory, single-matrix broadcasting, and every validation rule
//! that protects assembly from invalid scaling, epochs, matrix properties, pair definitions, and duplicate bodies.
BOOST_AUTO_TEST_CASE( test_InterArcStateContinuityConstraintSettings_PresetsAndValidation )
{
    const std::vector< double > epochs = { 100.0, 200.0 };

    // Position-only settings must populate all three Cartesian position entries and leave velocity unpenalized.
    auto positionOnly = simulation_setup::positionOnlyContinuity( { "Sat" }, { { "Sat", epochs } }, 2.5, 1.0 );
    const auto& positionWeightMatrix = positionOnly->weightMatrixForBodyAndPair( "Sat", 0 );
    BOOST_CHECK_CLOSE_FRACTION( positionWeightMatrix( 0, 0 ), 2.5, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( positionWeightMatrix( 1, 1 ), 2.5, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( positionWeightMatrix( 2, 2 ), 2.5, 1.0E-15 );
    BOOST_CHECK_EQUAL( positionWeightMatrix( 3, 3 ), 0.0 );
    BOOST_CHECK_EQUAL( positionWeightMatrix( 4, 4 ), 0.0 );
    BOOST_CHECK_EQUAL( positionWeightMatrix( 5, 5 ), 0.0 );

    // Velocity-only settings must invert the active block without leaking a weight into position.
    auto velocityOnly = simulation_setup::velocityOnlyContinuity( { "Sat" }, { { "Sat", epochs } }, 0.1 );
    const auto& velocityWeightMatrix = velocityOnly->weightMatrixForBodyAndPair( "Sat", 0 );
    BOOST_CHECK_EQUAL( velocityWeightMatrix( 0, 0 ), 0.0 );
    BOOST_CHECK_CLOSE_FRACTION( velocityWeightMatrix( 3, 3 ), 0.1, 1.0E-15 );

    // Full-state settings preserve independent position and velocity weights in their respective diagonal blocks.
    auto fullState = simulation_setup::fullStateContinuity( { "Sat" }, { { "Sat", epochs } }, 1.5, 0.7 );
    const auto& fullStateWeightMatrix = fullState->weightMatrixForBodyAndPair( "Sat", 0 );
    BOOST_CHECK_CLOSE_FRACTION( fullStateWeightMatrix( 0, 0 ), 1.5, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( fullStateWeightMatrix( 3, 3 ), 0.7, 1.0E-15 );

    // The preset stores one matrix but exposes both requested boundary pairs through broadcast semantics.
    BOOST_CHECK_EQUAL( positionOnly->numberOfPairsForBody( "Sat" ), 2u );

    // Zero, negative, and non-finite scaling factors cannot define a finite positive penalty denominator.
    BOOST_CHECK_THROW( simulation_setup::positionOnlyContinuity( { "Sat" }, { { "Sat", epochs } }, 1.0, 0.0 ), std::runtime_error );
    BOOST_CHECK_THROW( simulation_setup::positionOnlyContinuity( { "Sat" }, { { "Sat", epochs } }, 1.0, -1.0 ), std::runtime_error );
    BOOST_CHECK_THROW(
            simulation_setup::positionOnlyContinuity( { "Sat" }, { { "Sat", epochs } }, 1.0, std::numeric_limits< double >::infinity( ) ),
            std::runtime_error );

    // Epochs and weights must be finite to prevent invalid comparisons/eigendecompositions during assembly.
    BOOST_CHECK_THROW( simulation_setup::positionOnlyContinuity(
                               { "Sat" }, { { "Sat", { 100.0, std::numeric_limits< double >::quiet_NaN( ) } } }, 1.0, 1.0 ),
                       std::runtime_error );
    Eigen::Matrix< double, 6, 6 > nonFiniteWeight = Eigen::Matrix< double, 6, 6 >::Identity( );
    nonFiniteWeight( 0, 0 ) = std::numeric_limits< double >::infinity( );
    BOOST_CHECK_THROW( simulation_setup::generalContinuity(
                               { "Sat" },
                               { { "Sat", epochs } },
                               std::map< std::string, std::vector< Eigen::MatrixXd > >{ { "Sat", { nonFiniteWeight } } } ),
                       std::runtime_error );

    // Asymmetry above numerical tolerance is rejected because the matrix no longer defines the documented
    // self-adjoint quadratic weight.
    Eigen::Matrix< double, 6, 6 > asymmetric = Eigen::Matrix< double, 6, 6 >::Zero( );
    asymmetric( 0, 1 ) = 1.0;
    BOOST_CHECK_THROW(
            simulation_setup::generalContinuity( { "Sat" },
                                                 { { "Sat", epochs } },
                                                 std::map< std::string, std::vector< Eigen::MatrixXd > >{ { "Sat", { asymmetric } } } ),
            std::runtime_error );

    // A negative eigenvalue would make the continuity cost and information matrix indefinite.
    Eigen::Matrix< double, 6, 6 > indefinite = Eigen::Matrix< double, 6, 6 >::Zero( );
    indefinite( 0, 0 ) = -1.0;
    BOOST_CHECK_THROW(
            simulation_setup::generalContinuity( { "Sat" },
                                                 { { "Sat", epochs } },
                                                 std::map< std::string, std::vector< Eigen::MatrixXd > >{ { "Sat", { indefinite } } } ),
            std::runtime_error );

    // The feature is translational-only: generic weights must still be 6x6 position/velocity weights.
    Eigen::MatrixXd wrongSizeWeight = Eigen::Matrix3d::Identity( );
    BOOST_CHECK_THROW( simulation_setup::generalContinuity(
                               { "Sat" },
                               { { "Sat", epochs } },
                               std::map< std::string, std::vector< Eigen::MatrixXd > >{ { "Sat", { wrongSizeWeight } } } ),
                       std::runtime_error );

    // Explicit arc-pair metadata must provide exactly one pair for every connection epoch.
    BOOST_CHECK_THROW(
            InterArcStateContinuityConstraintSettings( { "Sat" },
                                                       { { "Sat", epochs } },
                                                       std::map< std::string, std::vector< Eigen::MatrixXd > >{
                                                               { "Sat",
                                                                 { tudat::simulation_setup::createCartesianStateWeightMatrix(
                                                                         Eigen::Vector3d::Ones( ), Eigen::Vector3d::Ones( ) ) } } },
                                                       1.0,
                                                       { { "Sat", { { 0, 1 } } } } ),
            std::runtime_error );

    // The implemented model connects adjacent arcs only; a gap in either supplied pair is invalid.
    BOOST_CHECK_THROW(
            InterArcStateContinuityConstraintSettings( { "Sat" },
                                                       { { "Sat", epochs } },
                                                       std::map< std::string, std::vector< Eigen::MatrixXd > >{
                                                               { "Sat",
                                                                 { tudat::simulation_setup::createCartesianStateWeightMatrix(
                                                                         Eigen::Vector3d::Ones( ), Eigen::Vector3d::Ones( ) ) } } },
                                                       1.0,
                                                       { { "Sat", { { 0, 2 }, { 1, 3 } } } } ),
            std::runtime_error );

    // Per-pair matrices must either broadcast from one entry or match the connection count exactly.
    BOOST_CHECK_THROW( InterArcStateContinuityConstraintSettings(
                               { "Sat" },
                               { { "Sat", epochs } },
                               std::map< std::string, std::vector< Eigen::MatrixXd > >{
                                       { "Sat",
                                         std::vector< Eigen::MatrixXd >( 3,
                                                                         tudat::simulation_setup::createCartesianStateWeightMatrix(
                                                                                 Eigen::Vector3d::Ones( ), Eigen::Vector3d::Ones( ) ) ) } },
                               1.0 ),
                       std::runtime_error );

    // Duplicate body names would apply the same settings twice while making output ordering ambiguous.
    BOOST_CHECK_THROW(
            InterArcStateContinuityConstraintSettings( { "Sat", "Sat" },
                                                       { { "Sat", epochs } },
                                                       std::map< std::string, std::vector< Eigen::MatrixXd > >{
                                                               { "Sat",
                                                                 { tudat::simulation_setup::createCartesianStateWeightMatrix(
                                                                         Eigen::Vector3d::Ones( ), Eigen::Vector3d::Ones( ) ) } } },
                                                       1.0 ),
            std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
