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

#define BOOST_TEST_MAIN

#include <limits>
#include <map>
#include <string>

#include <boost/test/included/unit_test.hpp>

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
TwoArcFixture buildTwoArcFixture( )
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
    // Adjacent arcs share one numerical boundary, which is the ambiguous lookup case under test.
    fixture.arcStartTimes = { initialEphemerisTime + 1.0E5, initialEphemerisTime + 1.0E5 + arcDuration };
    fixture.arcEndTimes = { fixture.arcStartTimes[ 0 ] + arcDuration, fixture.arcStartTimes[ 1 ] + arcDuration };

    std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKutta4Settings< double >( 600.0 );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > propagatorSettingsList;
    for( unsigned int i = 0; i < fixture.arcStartTimes.size( ); ++i )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialState =
                getInitialStateOfBody< double, double >( "Earth", "SSB", bodies, fixture.arcStartTimes[ i ] );
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

//! Verify boundary evaluation for propagated states and variational matrices: nearby state epochs must be
//! interpolated, while explicit arc selection must distinguish the two STMs at a shared numerical epoch.
BOOST_AUTO_TEST_CASE( test_ArcStateAndStmBoundaryEvaluation )
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

    // A propagated two-arc fixture checks the analogous shared-boundary behavior of the explicit STM accessor.
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
    BOOST_CHECK_SMALL( arc0Block.norm( ), 1.0E-12 );

    // Sensitivity block (columns fullStateTransitionSize..end) must be exactly zero at the arc start.
    if( fixture.fullSensitivitySize > 0 )
    {
        Eigen::MatrixXd sensitivityBlock = fullArc1.block( 0, fixture.fullStateTransitionSize, 6, fixture.fullSensitivitySize );
        BOOST_CHECK_SMALL( sensitivityBlock.norm( ), 1.0E-12 );
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
    BOOST_CHECK_SMALL( arc1BlockFromLeft.norm( ), 1.0E-12 );

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

//! Verify the complete quadratic model, including the discrepancy sign, normalized design matrix, dense positive
//! semi-definite weight,
//! cost, gradient, normal matrix, observation-count scaling, and scale-aware rank handling.
BOOST_AUTO_TEST_CASE( test_InterArcContinuityQuadraticModelAndScaling )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    const double connectionEpoch = fixture.arcStartTimes[ 1 ];

    // An outer product exercises a dense positive semi-definite matrix whose constrained dimension is one. The skew
    // perturbation is below validation tolerance but large enough to reveal use of the raw matrix in any one term.
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
    // that a dense rank-deficient positive semi-definite weight produces a finite positive semi-definite
    // information contribution.
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( contribution.additionalNormalMatrix );
    const double largestEigenvalue = std::max( solver.eigenvalues( ).maxCoeff( ), 1.0 );
    BOOST_CHECK_GE( solver.eigenvalues( ).minCoeff( ), -1.0E-9 * largestEigenvalue );
    BOOST_CHECK_GE( contribution.totalConstraintCost, 0.0 );

    // Observation averaging scales every term in the quadratic model by the same number of observations.
    const int numberOfObservations = 7;
    auto observationScaled = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                                       fixture.parametersToEstimate,
                                                                                       fixture.simulator,
                                                                                       fixture.stmInterface,
                                                                                       columnNormalizationFactors,
                                                                                       totalParameterSize,
                                                                                       numberOfObservations );
    BOOST_CHECK_SMALL( ( observationScaled.additionalNormalMatrix - numberOfObservations * contribution.additionalNormalMatrix ).norm( ),
                       1.0E-9 * std::max( observationScaled.additionalNormalMatrix.norm( ), 1.0 ) );
    BOOST_CHECK_SMALL( ( observationScaled.additionalRightHandSide - numberOfObservations * contribution.additionalRightHandSide ).norm( ),
                       1.0E-9 * std::max( observationScaled.additionalRightHandSide.norm( ), 1.0 ) );
    BOOST_CHECK_CLOSE_FRACTION( observationScaled.totalConstraintCost, numberOfObservations * contribution.totalConstraintCost, 1.0E-12 );

    // Relative rank detection must retain the three active position directions even at a very small absolute scale.
    auto tinySettings = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { connectionEpoch } } }, 1.0E-30, 1.0 );
    auto tinyContribution = assembleInterArcContinuityContribution< double, double >( { tinySettings },
                                                                                      fixture.parametersToEstimate,
                                                                                      fixture.simulator,
                                                                                      fixture.stmInterface,
                                                                                      Eigen::VectorXd::Ones( totalParameterSize ),
                                                                                      totalParameterSize );
    BOOST_CHECK_GT( tinyContribution.totalConstraintCost, 0.0 );
    BOOST_CHECK( std::isfinite( tinyContribution.totalConstraintCost ) );
}

//! Verify overlapping-arc evaluation and multi-body layout. Earth and Mars use different connection epochs, and
//! Mars' contribution is reconstructed from layout metadata to catch accidental use of the first body's rows.
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

    // A Mars-only reference permits direct inspection of its non-zero row offset at the interior overlap epoch.
    auto marsSettings = simulation_setup::positionOnlyContinuity( { "Mars" }, { { "Mars", { marsConnectionEpoch } } }, 1.0, 2.0 );
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

    // The shared settings path must return one finite six-component discrepancy per body in body-list order.
    BOOST_REQUIRE_EQUAL( multiBodyContribution.perPairDiscrepancies.size( ), 2u );
    BOOST_CHECK_EQUAL( multiBodyContribution.perPairDiscrepancies.at( 0 ).rows( ), 6 );
    BOOST_CHECK_EQUAL( multiBodyContribution.perPairDiscrepancies.at( 1 ).rows( ), 6 );
    BOOST_CHECK( std::isfinite( multiBodyContribution.perPairDiscrepancies.at( 0 ).norm( ) ) );
    BOOST_CHECK_LT( ( multiBodyContribution.perPairDiscrepancies.at( 1 ) - marsContribution.perPairDiscrepancies.at( 0 ) ).norm( ),
                    1.0E-12 );
}

//! Exercise covariance analysis and iterative estimation through one shared manager and observation setup.
//! Verify that the continuity prior adds positive semi-definite information, populates both output types, and
//! improves the arc jump.
BOOST_AUTO_TEST_CASE( test_EstimationAndCovariance_WithInterArcContinuity )
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

    // Compare covariance information at one linearization point before perturbing the parameters for estimation.
    auto unconstrainedCovarianceOutput =
            orbitDeterminationManager.computeCovariance( std::make_shared< CovarianceAnalysisInput< double, double > >( observations ) );
    auto constrainedCovarianceInput = std::make_shared< CovarianceAnalysisInput< double, double > >( observations );
    auto covarianceConstraint =
            simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { arcStartTimes[ 1 ] } } }, 1.0, 1.0E-15 );
    constrainedCovarianceInput->setInterArcContinuityConstraints( { covarianceConstraint } );
    auto constrainedCovarianceOutput = orbitDeterminationManager.computeCovariance( constrainedCovarianceInput );
    BOOST_REQUIRE( unconstrainedCovarianceOutput != nullptr );
    BOOST_REQUIRE( constrainedCovarianceOutput != nullptr );

    // Added continuity information must be non-zero and positive semi-definite up to round-off.
    Eigen::MatrixXd informationDifference = constrainedCovarianceOutput->inverseNormalizedCovarianceMatrix_ -
            unconstrainedCovarianceOutput->inverseNormalizedCovarianceMatrix_;
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > informationSolver( 0.5 *
                                                                        ( informationDifference + informationDifference.transpose( ) ) );
    const double largestInformationEigenvalue = std::max( informationSolver.eigenvalues( ).maxCoeff( ), 1.0 );
    BOOST_CHECK_GE( informationSolver.eigenvalues( ).minCoeff( ), -1.0E-9 * largestInformationEigenvalue );
    BOOST_CHECK_GT( informationSolver.eigenvalues( ).maxCoeff( ), 0.0 );
    BOOST_CHECK_GT( constrainedCovarianceOutput->getInterArcContinuityCost( ), 0.0 );
    BOOST_REQUIRE_EQUAL( constrainedCovarianceOutput->getInterArcContinuityDiscrepancies( ).size( ), 1u );

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

    // Enable a weak position prior for the iterative estimation from the perturbed initial states.
    auto estimationInputWithConstraint = std::make_shared< EstimationInput< double, double > >( observations );
    auto estimationConstraint =
            simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { arcStartTimes[ 1 ] } } }, 1.0, 1.0E14 );
    estimationInputWithConstraint->setInterArcContinuityConstraints( { estimationConstraint } );
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

    // Preserve the legacy positional constructor contract: the final pre-feature argument remains the propagation
    // exception flag, while continuity diagnostics retain their defaults.
    auto legacyOutput = CovarianceAnalysisOutput< double, double >( Eigen::MatrixXd::Identity( 1, 1 ),
                                                                    Eigen::VectorXd::Ones( 1 ),
                                                                    Eigen::VectorXd::Ones( 1 ),
                                                                    Eigen::MatrixXd::Identity( 1, 1 ),
                                                                    Eigen::MatrixXd::Zero( 0, 0 ),
                                                                    Eigen::VectorXd::Zero( 0 ),
                                                                    Eigen::MatrixXd::Zero( 0, 0 ),
                                                                    Eigen::MatrixXd::Zero( 0, 0 ),
                                                                    true );
    BOOST_CHECK( legacyOutput.exceptionDuringPropagation_ );
    BOOST_CHECK_EQUAL( legacyOutput.getInterArcContinuityCost( ), 0.0 );
}

//! Invalid configurations must fail before assembly. Cover an epoch outside the selected arcs and an ambiguous
//! parameter set containing the same body's multi-arc initial state twice.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_InvalidConfigurationDiagnostics )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );

    // Place the epoch beyond both arcs so whichever side is validated first must reject it.
    const double outOfRangeConnectionEpoch = fixture.arcEndTimes[ 1 ] + 1.0E6;
    auto badSettings = simulation_setup::positionOnlyContinuity( { "Earth" }, { { "Earth", { outOfRangeConnectionEpoch } } }, 1.0, 1.0 );

    BOOST_CHECK_THROW( ( assembleInterArcContinuityContribution< double, double >( { badSettings },
                                                                                   fixture.parametersToEstimate,
                                                                                   fixture.simulator,
                                                                                   fixture.stmInterface,
                                                                                   columnNormalizationFactors,
                                                                                   totalParameterSize ) ),
                       std::runtime_error );

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
    BOOST_CHECK_THROW( ( assembleInterArcContinuityContribution< double, double >( { validSettings },
                                                                                   duplicatedParametersToEstimate,
                                                                                   fixture.simulator,
                                                                                   fixture.stmInterface,
                                                                                   columnNormalizationFactors,
                                                                                   totalParameterSize ) ),
                       std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
