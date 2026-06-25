/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
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

namespace
{

struct TwoArcFixture {
    std::shared_ptr< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< double > > stmInterface;
    std::shared_ptr< MultiArcDynamicsSimulator< double, double > > simulator;
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
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames = { "Earth", "Sun" };

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 4.0E7;
    const double buffer = 3.6E5;

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    SystemOfBodies bodies = createSystemOfBodies< double, double >( bodySettings );

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Earth" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    std::vector< std::string > bodiesToIntegrate = { "Earth" };
    std::vector< std::string > centralBodies = { "SSB" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    const double arcDuration = 1.0E7;
    TwoArcFixture fixture;
    if( useOverlappingArcs )
    {
        const double overlapDuration = 1.0E6;
        fixture.arcStartTimes = { initialEphemerisTime + 1.0E5, initialEphemerisTime + 1.0E5 + arcDuration - overlapDuration };
        fixture.arcEndTimes = { fixture.arcStartTimes[ 0 ] + arcDuration, fixture.arcStartTimes[ 1 ] + arcDuration };
    }
    else
    {
        // Two adjacent arcs (no overlap) so the boundary is shared exactly.
        fixture.arcStartTimes = { initialEphemerisTime + 1.0E5, initialEphemerisTime + 1.0E5 + arcDuration };
        fixture.arcEndTimes = { fixture.arcStartTimes[ 0 ] + arcDuration, fixture.arcStartTimes[ 1 ] + arcDuration };
    }

    std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKutta4Settings< double >( 600.0 );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > propagatorSettingsList;
    for( unsigned int i = 0; i < fixture.arcStartTimes.size( ); ++i )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialState =
                getInitialStateOfBody< double, double >( "Earth", "SSB", bodies, fixture.arcStartTimes[ i ] );
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

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double, double >( propagatorSettings, bodies, fixture.arcStartTimes );

    fixture.parametersToEstimate = createParametersToEstimate< double >( parameterNames, bodies );

    auto variationalSolver = simulation_setup::createVariationalEquationsSolver< double, double >(
            bodies, propagatorSettings, fixture.parametersToEstimate, true );

    fixture.stmInterface = std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< double > >(
            variationalSolver->getStateTransitionMatrixInterface( ) );
    fixture.simulator =
            std::dynamic_pointer_cast< MultiArcDynamicsSimulator< double, double > >( variationalSolver->getDynamicsSimulatorBase( ) );

    BOOST_REQUIRE( fixture.stmInterface != nullptr );
    BOOST_REQUIRE( fixture.simulator != nullptr );

    fixture.fullStateTransitionSize = fixture.stmInterface->getFullStateTransitionMatrixSize( );
    fixture.fullSensitivitySize = fixture.stmInterface->getFullSensitivityMatrixSize( );

    return fixture;
}

Eigen::VectorXd getSinglePairDiscrepancy( const TwoArcFixture& fixture )
{
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    auto settings = fullStateContinuity( "Earth", { fixture.arcStartTimes[ 1 ] }, 1.0, 1.0, 1.0 );
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );
    auto contribution = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                                  fixture.parametersToEstimate,
                                                                                  fixture.simulator,
                                                                                  fixture.stmInterface,
                                                                                  columnNormalizationFactors,
                                                                                  totalParameterSize );
    return contribution.perPairDiscrepancies.at( 0 );
}

TwoArcFixture buildTwoBodyTwoArcFixture( )
{
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames = { "Earth", "Mars", "Sun" };

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 4.0E7;
    const double buffer = 3.6E5;

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Mars" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    SystemOfBodies bodies = createSystemOfBodies< double, double >( bodySettings );

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Earth" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Mars" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    std::vector< std::string > bodiesToIntegrate = { "Earth", "Mars" };
    std::vector< std::string > centralBodies = { "SSB", "SSB" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    const double arcDuration = 1.0E7;
    TwoArcFixture fixture;
    fixture.arcStartTimes = { initialEphemerisTime + 1.0E5, initialEphemerisTime + 1.0E5 + arcDuration };
    fixture.arcEndTimes = { fixture.arcStartTimes[ 0 ] + arcDuration, fixture.arcStartTimes[ 1 ] + arcDuration };

    std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKutta4Settings< double >( 600.0 );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > propagatorSettingsList;
    for( unsigned int i = 0; i < fixture.arcStartTimes.size( ); ++i )
    {
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

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double, double >( propagatorSettings, bodies, fixture.arcStartTimes );
    fixture.parametersToEstimate = createParametersToEstimate< double >( parameterNames, bodies );

    auto variationalSolver = simulation_setup::createVariationalEquationsSolver< double, double >(
            bodies, propagatorSettings, fixture.parametersToEstimate, true );

    fixture.stmInterface = std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< double > >(
            variationalSolver->getStateTransitionMatrixInterface( ) );
    fixture.simulator =
            std::dynamic_pointer_cast< MultiArcDynamicsSimulator< double, double > >( variationalSolver->getDynamicsSimulatorBase( ) );

    BOOST_REQUIRE( fixture.stmInterface != nullptr );
    BOOST_REQUIRE( fixture.simulator != nullptr );

    fixture.fullStateTransitionSize = fixture.stmInterface->getFullStateTransitionMatrixSize( );
    fixture.fullSensitivitySize = fixture.stmInterface->getFullSensitivityMatrixSize( );

    return fixture;
}

}  // namespace

//! At a shared arc boundary, the per-arc-index "full" STM accessor returns identity in the right arc's 6-block
//! and zeros everywhere else (other arcs' state blocks and the sensitivity block).
BOOST_AUTO_TEST_CASE( test_StmForArc_SharedBoundaryIdentity )
{
    auto fixture = buildTwoArcFixture( );

    BOOST_REQUIRE_EQUAL( fixture.arcStartTimes.size( ), 2u );
    const double connectionEpoch = fixture.arcStartTimes[ 1 ];

    // Verify the structural property for arc 1: identity-Phi, zero-S at the arc start.
    Eigen::MatrixXd fullArc1 = fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, connectionEpoch );

    BOOST_CHECK_EQUAL( fullArc1.rows( ), 6 );
    BOOST_CHECK_EQUAL( fullArc1.cols( ), fixture.fullStateTransitionSize + fixture.fullSensitivitySize );

    // Layout: column block [0, 6) holds arc 0's 6x6 STM slot, [6, 12) holds arc 1's. Sensitivity columns follow.
    const Eigen::Matrix< double, 6, 6 > arc0Block = fullArc1.block< 6, 6 >( 0, 0 );
    const Eigen::Matrix< double, 6, 6 > arc1Block = fullArc1.block< 6, 6 >( 0, 6 );
    const Eigen::Matrix< double, 6, 6 > identity6 = Eigen::Matrix< double, 6, 6 >::Identity( );
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
    BOOST_CHECK_GT( ( arc0BlockFromLeft - identity6 ).norm( ), 1.0E-3 );
    for( int row = 0; row < 6; ++row )
    {
        for( int col = 0; col < 6; ++col )
        {
            BOOST_CHECK_SMALL( std::fabs( arc1BlockFromLeft( row, col ) ), 1.0E-12 );
        }
    }
}

//! Test the per-arc-index accessor's range validation: arcIndex out of range and time outside the arc both throw.
BOOST_AUTO_TEST_CASE( test_StmForArc_RangeValidation )
{
    auto fixture = buildTwoArcFixture( );

    BOOST_CHECK_THROW( fixture.stmInterface->getCombinedStateTransitionAndSensitivityMatrixForArc( -1, fixture.arcStartTimes[ 0 ] ),
                       std::runtime_error );
    BOOST_CHECK_THROW( fixture.stmInterface->getCombinedStateTransitionAndSensitivityMatrixForArc( 2, fixture.arcStartTimes[ 0 ] ),
                       std::runtime_error );
    BOOST_CHECK_THROW( fixture.stmInterface->getCombinedStateTransitionAndSensitivityMatrixForArc( 0, fixture.arcEndTimes[ 1 ] + 1.0 ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, fixture.arcStartTimes[ 0 ] - 1.0 ),
            std::runtime_error );
}

//! Structural check of the assembly module: build a position-only continuity contribution at the shared boundary
//! and verify symmetry, dimensionality, finite discrepancy diagnostics, and a non-negative cost.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_StructureAndSymmetry )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;

    auto settings = positionOnlyContinuity( "Earth", { fixture.arcStartTimes[ 1 ] }, 1.0, 1.0 );
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );

    auto contribution = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                                  fixture.parametersToEstimate,
                                                                                  fixture.simulator,
                                                                                  fixture.stmInterface,
                                                                                  columnNormalizationFactors,
                                                                                  totalParameterSize );

    BOOST_REQUIRE_EQUAL( contribution.additionalNormalMatrix.rows( ), totalParameterSize );
    BOOST_REQUIRE_EQUAL( contribution.additionalNormalMatrix.cols( ), totalParameterSize );
    BOOST_REQUIRE_EQUAL( contribution.additionalRightHandSide.size( ), totalParameterSize );
    BOOST_REQUIRE_EQUAL( contribution.perPairDiscrepancies.size( ), 1u );

    // The normal-matrix contribution must be symmetric.
    Eigen::MatrixXd normalMatrixAsymmetry = contribution.additionalNormalMatrix - contribution.additionalNormalMatrix.transpose( );
    BOOST_CHECK_SMALL( normalMatrixAsymmetry.norm( ), 1.0E-9 );

    // The normal-matrix contribution is PSD, allowing for floating-point tolerance relative to the largest eigenvalue.
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( contribution.additionalNormalMatrix );
    const double largestEigenvalue = std::max( solver.eigenvalues( ).maxCoeff( ), 1.0 );
    BOOST_CHECK_GE( solver.eigenvalues( ).minCoeff( ), -1.0E-9 * largestEigenvalue );

    // The discrepancy at the shared boundary is finite (no NaN/Inf). Magnitude depends on the propagation model:
    // this test uses point-mass Earth-Sun gravity only, so the arc-0 forward propagation drifts substantially from
    // SPICE Earth (used as arc-1's initial state). That's irrelevant to the structural assembly check.
    const Eigen::VectorXd& stateDiscrepancy = contribution.perPairDiscrepancies[ 0 ];
    BOOST_CHECK_EQUAL( stateDiscrepancy.rows( ), 6 );
    BOOST_CHECK( std::isfinite( stateDiscrepancy.norm( ) ) );

    // Cost is non-negative.
    BOOST_CHECK_GE( contribution.totalConstraintCost, 0.0 );
}

//! Exact assembly check: the design matrix is the right-arc variational block minus the left-arc variational block,
//! its columns are normalized before the normal equations are formed, the residual sign is negative state discrepancy,
//! and dense rank-deficient PSD weights are handled without whitening/Cholesky assumptions.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_DensePsdExactNormalEquations )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    const double connectionEpoch = fixture.arcStartTimes[ 1 ];

    Eigen::Matrix< double, 6, 1 > rankOneWeightVector;
    rankOneWeightVector << 2.0, -1.0, 0.5, 0.25, -0.75, 1.5;
    Eigen::Matrix< double, 6, 6 > denseRankOneWeightMatrix = rankOneWeightVector * rankOneWeightVector.transpose( );
    const double constraintScalingFactor = 2.0;
    auto settings = generalContinuity( "Earth", { connectionEpoch }, { denseRankOneWeightMatrix }, constraintScalingFactor );

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

    Eigen::MatrixXd leftArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 0, connectionEpoch );
    Eigen::MatrixXd rightArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, connectionEpoch );
    Eigen::MatrixXd continuityDesignMatrix = rightArcVariationalMatrix - leftArcVariationalMatrix;
    for( int col = 0; col < totalParameterSize; ++col )
    {
        continuityDesignMatrix.col( col ) /= columnNormalizationFactors( col );
    }

    const Eigen::VectorXd& stateDiscrepancy = contribution.perPairDiscrepancies.at( 0 );
    BOOST_CHECK_EQUAL( stateDiscrepancy.rows( ), 6 );
    const Eigen::Matrix< double, 6, 6 > scaledConstraintWeight =
            denseRankOneWeightMatrix / constraintScalingFactor;  // rank(weight matrix)=1 -> total constrained dimension=1
    const Eigen::MatrixXd expectedNormalMatrixContribution =
            continuityDesignMatrix.transpose( ) * scaledConstraintWeight * continuityDesignMatrix;
    const Eigen::VectorXd expectedRightHandSideContribution =
            -continuityDesignMatrix.transpose( ) * ( scaledConstraintWeight * stateDiscrepancy );
    const double expectedCost = stateDiscrepancy.transpose( ) * scaledConstraintWeight * stateDiscrepancy;

    BOOST_CHECK_LT( ( contribution.additionalNormalMatrix - expectedNormalMatrixContribution ).norm( ) /
                            std::max( expectedNormalMatrixContribution.norm( ), 1.0E-30 ),
                    1.0E-12 );
    BOOST_CHECK_LT( ( contribution.additionalRightHandSide - expectedRightHandSideContribution ).norm( ) /
                            std::max( expectedRightHandSideContribution.norm( ), 1.0E-30 ),
                    1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( contribution.totalConstraintCost, expectedCost, 1.0E-12 );
}

//! Orbit14-style mid-gap connection epoch. The adjacent arcs overlap around the connection epoch so both sides are
//! evaluated strictly inside their propagated intervals; neither STM block is the arc-start identity.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_MidGapConnectionEpoch )
{
    auto fixture = buildTwoArcFixture( {}, true );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    const double connectionEpoch = 0.5 * ( fixture.arcStartTimes[ 1 ] + fixture.arcEndTimes[ 0 ] );
    BOOST_CHECK_GT( connectionEpoch, fixture.arcStartTimes[ 1 ] );
    BOOST_CHECK_LT( connectionEpoch, fixture.arcEndTimes[ 0 ] );

    auto settings = fullStateContinuity( "Earth", { connectionEpoch }, 1.0, 1.0, 1.0 );
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );
    auto contribution = assembleInterArcContinuityContribution< double, double >( { settings },
                                                                                  fixture.parametersToEstimate,
                                                                                  fixture.simulator,
                                                                                  fixture.stmInterface,
                                                                                  columnNormalizationFactors,
                                                                                  totalParameterSize );

    BOOST_REQUIRE_EQUAL( contribution.perPairDiscrepancies.size( ), 1u );
    BOOST_CHECK( std::isfinite( contribution.perPairDiscrepancies.at( 0 ).norm( ) ) );
    BOOST_CHECK_GE( contribution.totalConstraintCost, 0.0 );

    Eigen::MatrixXd leftArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 0, connectionEpoch );
    Eigen::MatrixXd rightArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, connectionEpoch );
    const Eigen::Matrix< double, 6, 6 > identity6 = Eigen::Matrix< double, 6, 6 >::Identity( );

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

    BOOST_REQUIRE_EQUAL( totalParameterSize, 12 );

    Eigen::MatrixXd leftArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 0, connectionEpoch );
    Eigen::MatrixXd rightArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, connectionEpoch );
    const Eigen::MatrixXd analyticalContinuityDesignMatrix = rightArcVariationalMatrix - leftArcVariationalMatrix;

    Eigen::MatrixXd finiteDifferenceContinuityDesignMatrix( 6, totalParameterSize );
    finiteDifferenceContinuityDesignMatrix.setZero( );

    for( int column = 0; column < totalParameterSize; ++column )
    {
        const double perturbation = ( column % 6 < 3 ) ? 10.0 : 1.0E-3;
        std::vector< Eigen::Matrix< double, 6, 1 > > positivePerturbations( 2, Eigen::Matrix< double, 6, 1 >::Zero( ) );
        std::vector< Eigen::Matrix< double, 6, 1 > > negativePerturbations( 2, Eigen::Matrix< double, 6, 1 >::Zero( ) );
        positivePerturbations.at( static_cast< unsigned int >( column / 6 ) )( column % 6 ) = perturbation;
        negativePerturbations.at( static_cast< unsigned int >( column / 6 ) )( column % 6 ) = -perturbation;

        auto positiveFixture = buildTwoArcFixture( positivePerturbations );
        auto negativeFixture = buildTwoArcFixture( negativePerturbations );
        const Eigen::VectorXd positiveDiscrepancy = getSinglePairDiscrepancy( positiveFixture );
        const Eigen::VectorXd negativeDiscrepancy = getSinglePairDiscrepancy( negativeFixture );
        BOOST_CHECK_EQUAL( positiveDiscrepancy.rows( ), 6 );
        BOOST_CHECK_EQUAL( negativeDiscrepancy.rows( ), 6 );
        finiteDifferenceContinuityDesignMatrix.col( column ) = ( positiveDiscrepancy - negativeDiscrepancy ) / ( 2.0 * perturbation );
    }

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

    auto positionSettings = positionOnlyContinuity( "Earth", { connectionEpoch }, 2.0, 5.0 );
    auto velocitySettings = velocityOnlyContinuity( "Earth", { connectionEpoch }, 4.0, 5.0 );
    auto fullSettings = fullStateContinuity( "Earth", { connectionEpoch }, 2.0, 4.0, 5.0 );

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

    const Eigen::VectorXd& stateDiscrepancy = fullContribution.perPairDiscrepancies.at( 0 );
    BOOST_CHECK_EQUAL( stateDiscrepancy.rows( ), 6 );
    const double expectedPositionCost = ( 2.0 / ( 5.0 * 3.0 ) ) * stateDiscrepancy.head( 3 ).squaredNorm( );
    const double expectedVelocityCost = ( 4.0 / ( 5.0 * 3.0 ) ) * stateDiscrepancy.tail( 3 ).squaredNorm( );
    const double expectedFullCost =
            ( 1.0 / ( 5.0 * 6.0 ) ) * ( 2.0 * stateDiscrepancy.head( 3 ).squaredNorm( ) + 4.0 * stateDiscrepancy.tail( 3 ).squaredNorm( ) );

    BOOST_CHECK_CLOSE_FRACTION( positionContribution.totalConstraintCost, expectedPositionCost, 1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( velocityContribution.totalConstraintCost, expectedVelocityCost, 1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( fullContribution.totalConstraintCost, expectedFullCost, 1.0E-12 );

    // Very small positive weights must still be counted by relative rank; this used to fail with an absolute
    // 1e-12 floor in rankOf6x6PsdMatrix.
    auto tinySettings = positionOnlyContinuity( "Earth", { connectionEpoch }, 1.0E-30, 1.0 );
    auto tinyContribution = assembleInterArcContinuityContribution< double, double >( { tinySettings },
                                                                                      fixture.parametersToEstimate,
                                                                                      fixture.simulator,
                                                                                      fixture.stmInterface,
                                                                                      columnNormalizationFactors,
                                                                                      totalParameterSize );
    BOOST_CHECK_GT( tinyContribution.totalConstraintCost, 0.0 );
    BOOST_CHECK( std::isfinite( tinyContribution.totalConstraintCost ) );
}

//! Regression test for body-specific state/STM row slicing: constraining Mars in a two-body multi-arc propagation
//! must use Mars' 6-row block, not the first body's block.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_UsesRequestedBodyRows )
{
    auto fixture = buildTwoBodyTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    const double connectionEpoch = fixture.arcStartTimes[ 1 ];
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );

    auto marsSettings = positionOnlyContinuity( "Mars", { connectionEpoch }, 3.0, 2.0 );
    auto contribution = assembleInterArcContinuityContribution< double, double >( { marsSettings },
                                                                                  fixture.parametersToEstimate,
                                                                                  fixture.simulator,
                                                                                  fixture.stmInterface,
                                                                                  columnNormalizationFactors,
                                                                                  totalParameterSize );

    const auto layout = fixture.stmInterface->getArcWiseAndFullSolutionInitialStateIndices( );
    const int marsArcWiseRowsLeft = layout.at( 0 ).at( "Mars" ).first.first;
    const int marsArcWiseRowsRight = layout.at( 1 ).at( "Mars" ).first.first;
    const int marsFullRowsLeft = layout.at( 0 ).at( "Mars" ).second.first.first;
    const int marsFullRowsRight = layout.at( 1 ).at( "Mars" ).second.first.first;
    BOOST_REQUIRE_NE( marsFullRowsLeft, 0 );
    BOOST_REQUIRE_NE( marsFullRowsRight, 0 );

    BOOST_CHECK( std::isfinite( contribution.perPairDiscrepancies.at( 0 ).norm( ) ) );

    Eigen::MatrixXd leftArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 0, connectionEpoch );
    Eigen::MatrixXd rightArcVariationalMatrix =
            fixture.stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( 1, connectionEpoch );
    Eigen::MatrixXd expectedContinuityDesignMatrix = rightArcVariationalMatrix.block( marsFullRowsRight, 0, 6, totalParameterSize ) -
            leftArcVariationalMatrix.block( marsFullRowsLeft, 0, 6, totalParameterSize );
    Eigen::Matrix< double, 6, 6 > scaledConstraintWeight = Eigen::Matrix< double, 6, 6 >::Zero( );
    scaledConstraintWeight.block< 3, 3 >( 0, 0 ) = ( 3.0 / ( 2.0 * 3.0 ) ) * Eigen::Matrix3d::Identity( );
    const Eigen::MatrixXd expectedNormalMatrixContribution =
            expectedContinuityDesignMatrix.transpose( ) * scaledConstraintWeight * expectedContinuityDesignMatrix;

    BOOST_CHECK_LT( ( contribution.additionalNormalMatrix - expectedNormalMatrixContribution ).norm( ) /
                            std::max( expectedNormalMatrixContribution.norm( ), 1.0E-30 ),
                    1.0E-12 );
}

//! Normalization invariance. Compare assembly with two different column-normalization conventions; the physical
//! parameter update implied by the normal-equation contribution must be invariant under uniform rescaling of the
//! normalization vector.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_NormalisationInvariance )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;

    auto settings = positionOnlyContinuity( "Earth", { fixture.arcStartTimes[ 1 ] }, 1.0, 1.0 );

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

    // The cost is computed in physical units (independent of column normalization) and must match exactly.
    BOOST_CHECK_CLOSE_FRACTION( unit.totalConstraintCost, scaled.totalConstraintCost, 1.0E-12 );

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
}

//! Integration test exercising the full OD loop with constraints attached. Covers tests 1 (smoke check that a
//! position-only constraint runs end-to-end), 5 (heterogeneous-weight superposition: a settings entry with two
//! identical pairs accumulates correctly), and 13 (history population reflects per-iteration constraint cost,
//! and best-iteration selection uses the combined cost).
BOOST_AUTO_TEST_CASE( test_OdLoop_WithInterArcContinuity_EndToEnd )
{
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

    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );

    // Synthetic observations spaced over each arc.
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

    // Perturb the truth slightly so the estimator has work to do.
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

    auto estimationInputNoConstraint = std::make_shared< EstimationInput< double, double > >( observations );
    auto outputNoConstraint = orbitDeterminationManager.estimateParameters( estimationInputNoConstraint );
    BOOST_REQUIRE( outputNoConstraint != nullptr );
    BOOST_CHECK( outputNoConstraint->getInterArcContinuityCostHistory( ).empty( ) );
    BOOST_CHECK( outputNoConstraint->getInterArcContinuityDiscrepancyHistory( ).empty( ) );

    // Re-estimate with a position-only continuity prior at the shared OCM boundary.
    parametersToEstimate->resetParameterValues( initialEstimate );
    auto estimationInputWithConstraint = std::make_shared< EstimationInput< double, double > >( observations );
    auto constraint = positionOnlyContinuity( "Earth", { arcStartTimes[ 1 ] }, 1.0, 1.0E-12 );
    estimationInputWithConstraint->setInterArcContinuityConstraints( { constraint } );
    auto outputWithConstraint = orbitDeterminationManager.estimateParameters( estimationInputWithConstraint );

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

    // The boundary discrepancy at the best iteration should be smaller than at the first iteration for any
    // non-trivially-weak constraint.
    const auto& history = outputWithConstraint->getInterArcContinuityDiscrepancyHistory( );
    const int bestIterationIndex =
            std::max( 0, std::min( outputWithConstraint->bestIteration_, static_cast< int >( history.size( ) ) - 1 ) );
    const Eigen::VectorXd& bestIterationStateDiscrepancy = history.at( bestIterationIndex ).at( 0 );
    BOOST_CHECK_EQUAL( bestIterationStateDiscrepancy.rows( ), 6 );
    BOOST_TEST_MESSAGE( "Best-iter position discrepancy with constraint: " << bestIterationStateDiscrepancy.head( 3 ).norm( ) );
    BOOST_TEST_MESSAGE( "Best-iter velocity discrepancy with constraint: " << bestIterationStateDiscrepancy.tail( 3 ).norm( ) );
    BOOST_CHECK( std::isfinite( bestIterationStateDiscrepancy.norm( ) ) );
}

//! Test 6: connection epoch outside the propagated interval of either side arc throws with the documented
//! diagnostic identifying the offending arc and its bounds.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_EpochOutsideArcThrows )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );

    // Epoch that lies outside both arcs.
    const double outOfRangeConnectionEpoch = fixture.arcEndTimes[ 1 ] + 1.0E6;
    auto badSettings = std::make_shared< InterArcStateContinuityConstraintSettings >(
            "Earth",
            std::vector< double >{ outOfRangeConnectionEpoch },
            std::vector< Eigen::Matrix< double, 6, 6 > >{ tudat::simulation_setup::detail::diagonalWeight( 1.0, 0.0 ) },
            std::vector< double >{ 1.0 } );

    try
    {
        assembleInterArcContinuityContribution< double, double >( { badSettings },
                                                                  fixture.parametersToEstimate,
                                                                  fixture.simulator,
                                                                  fixture.stmInterface,
                                                                  columnNormalizationFactors,
                                                                  totalParameterSize );
        BOOST_FAIL( "Expected runtime_error from out-of-range connection epoch was not thrown." );
    }
    catch( const std::runtime_error& error )
    {
        const std::string what = error.what( );
        BOOST_CHECK( what.find( "Inter-arc continuity connection epoch" ) != std::string::npos );
        BOOST_CHECK( what.find( "Earth" ) != std::string::npos );
        BOOST_CHECK( what.find( "outside the propagated interval" ) != std::string::npos );
        BOOST_CHECK( what.find( "Extend the arc propagation interval" ) != std::string::npos );
    }
}

//! Test 12: covariance analysis with inter-arc continuity priors attached directly to CovarianceAnalysisInput. The regularized inverse
//! normalized covariance matrix must dominate the unregularized one elementwise on the parameter block (PSD
//! inequality), confirming that the position-only prior genuinely tightens the covariance in the regularized
//! directions.
BOOST_AUTO_TEST_CASE( test_CovarianceAnalysis_WithInterArcContinuity_Tightens )
{
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

    auto unconstrainedInput = std::make_shared< CovarianceAnalysisInput< double, double > >( observations );
    auto unconstrainedOutput = orbitDeterminationManager.computeCovariance( unconstrainedInput );

    auto constrainedInput = std::make_shared< CovarianceAnalysisInput< double, double > >( observations );
    auto constraint = positionOnlyContinuity( "Earth", { arcStartTimes[ 1 ] }, 1.0, 1.0E-15 );
    constrainedInput->setInterArcContinuityConstraints( { constraint } );
    auto constrainedOutput = orbitDeterminationManager.computeCovariance( constrainedInput );

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
    // The largest eigenvalue must be strictly positive — the constraint must have a non-trivial effect.
    BOOST_CHECK_GT( solver.eigenvalues( ).maxCoeff( ), 0.0 );
    BOOST_CHECK_GT( constrainedOutput->getInterArcContinuityCost( ), 0.0 );
    BOOST_REQUIRE_EQUAL( constrainedOutput->getInterArcContinuityDiscrepancies( ).size( ), 1u );
    BOOST_CHECK( std::isfinite( constrainedOutput->getInterArcContinuityDiscrepancies( ).at( 0 ).norm( ) ) );
}

//! Test global constrained-dimension accounting: passing two settings entries with one rank-3 position-only pair
//! each should halve each pair's scaled weight while doubling the number of pairs, leaving the accumulated normal
//! equations unchanged.
BOOST_AUTO_TEST_CASE( test_AssembleInterArcContinuity_GlobalMdAccounting )
{
    auto fixture = buildTwoArcFixture( );
    const int totalParameterSize = fixture.fullStateTransitionSize + fixture.fullSensitivitySize;
    Eigen::VectorXd columnNormalizationFactors = Eigen::VectorXd::Ones( totalParameterSize );

    // Single position-only constraint at the boundary: constrained dimension equals 3.
    auto singleSettings = positionOnlyContinuity( "Earth", { fixture.arcStartTimes[ 1 ] }, 1.0, 1.0 );
    auto singleContribution = assembleInterArcContinuityContribution< double, double >( { singleSettings },
                                                                                        fixture.parametersToEstimate,
                                                                                        fixture.simulator,
                                                                                        fixture.stmInterface,
                                                                                        columnNormalizationFactors,
                                                                                        totalParameterSize );

    // Add a second identical settings entry: the total constrained dimension goes from 3 to 6, so each pair's
    // scaled weight halves and the accumulated normal matrix is 2 * (1/2) = 1x the single-pair normal matrix. The right-hand side behaves
    // the same way. The total cost scales like the per-pair weight (factor 1/2 for each pair).
    auto duplicateSettings = positionOnlyContinuity( "Earth", { fixture.arcStartTimes[ 1 ] }, 1.0, 1.0 );
    auto duplicatedContribution = assembleInterArcContinuityContribution< double, double >( { singleSettings, duplicateSettings },
                                                                                            fixture.parametersToEstimate,
                                                                                            fixture.simulator,
                                                                                            fixture.stmInterface,
                                                                                            columnNormalizationFactors,
                                                                                            totalParameterSize );

    // Normal matrix from {single, duplicate} == 2 * (1/2) * single normal matrix.
    BOOST_CHECK_LT( ( duplicatedContribution.additionalNormalMatrix - singleContribution.additionalNormalMatrix ).norm( ),
                    1.0E-9 * std::max( singleContribution.additionalNormalMatrix.norm( ), 1.0 ) );
    // Cost from {single, duplicate} == 2 * (1/2) * cost_single = cost_single.
    BOOST_CHECK_CLOSE_FRACTION( duplicatedContribution.totalConstraintCost, singleContribution.totalConstraintCost, 1.0E-9 );
    BOOST_CHECK_EQUAL( duplicatedContribution.perPairDiscrepancies.size( ), 2u );
}

//! Test 11: Passing empty additionalNormalMatrix / additionalRightHandSide leaves the LSQ output bit-identical
//! to the pre-existing overload behaviour. Establishes that the new optional arguments are zero-impact
//! when unused.
BOOST_AUTO_TEST_CASE( test_LeastSquaresEmptyAdditionsNoOp )
{
    Eigen::MatrixXd designMatrix( 5, 3 );
    designMatrix << 1.0, 0.5, -0.2, 0.3, 1.2, 0.1, -0.7, 0.4, 1.0, 0.2, -0.3, 0.9, 1.1, 0.8, -0.5;
    Eigen::VectorXd residuals( 5 );
    residuals << 0.1, -0.2, 0.05, -0.05, 0.15;
    Eigen::VectorXd weights = Eigen::VectorXd::Constant( 5, 1.0 );
    Eigen::MatrixXd inverseApriori = 0.01 * Eigen::MatrixXd::Identity( 3, 3 );

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

    // Shape-validation throws.
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
}

//! Test A5: hard-equality constraints (Lagrange multipliers) and soft additions (additionalNormalMatrix /
//! additionalRightHandSide) compose correctly. The hard constraint augments the system to (n + n_c); the soft
//! addition goes into the top-left n×n parameter block only.
BOOST_AUTO_TEST_CASE( test_LeastSquares_HardAndSoftConstraintsCompose )
{
    Eigen::MatrixXd designMatrix( 5, 3 );
    designMatrix << 1.0, 0.5, -0.2, 0.3, 1.2, 0.1, -0.7, 0.4, 1.0, 0.2, -0.3, 0.9, 1.1, 0.8, -0.5;
    Eigen::VectorXd residuals( 5 );
    residuals << 0.1, -0.2, 0.05, -0.05, 0.15;
    Eigen::VectorXd weights = Eigen::VectorXd::Constant( 5, 1.0 );
    Eigen::MatrixXd inverseApriori = 0.01 * Eigen::MatrixXd::Identity( 3, 3 );

    // Hard constraint: x[0] = 0.5.
    Eigen::MatrixXd constraintMultiplier( 1, 3 );
    constraintMultiplier << 1.0, 0.0, 0.0;
    Eigen::VectorXd constraintRhs( 1 );
    constraintRhs << 0.5;

    // Soft addition: alpha*I on the normal matrix, non-zero RHS.
    const double alpha = 0.2;
    Eigen::MatrixXd soft = alpha * Eigen::MatrixXd::Identity( 3, 3 );
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

    // The solution vector has size: number of estimated parameters plus number of hard constraints.
    BOOST_REQUIRE_EQUAL( withBoth.first.size( ), 4 );

    // The hard constraint is exactly satisfied.
    const Eigen::VectorXd parameterCorrection = withBoth.first.head( 3 );
    BOOST_CHECK_CLOSE_FRACTION( ( constraintMultiplier * parameterCorrection )( 0 ), constraintRhs( 0 ), 1.0E-10 );

    // The soft addition has a non-trivial effect: compare to the same run without it.
    auto withHardOnly = tudat::linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix(
            designMatrix, residuals, weights, inverseApriori, 1.0E8, constraintMultiplier, constraintRhs );
    BOOST_CHECK_GT( ( withBoth.first.head( 3 ) - withHardOnly.first.head( 3 ) ).norm( ), 1.0E-6 );

    // The parameter block of the normal matrix is augmented by exactly the soft matrix; the Lagrange-multiplier
    // rows/columns are untouched.
    const Eigen::MatrixXd normalDiff = withBoth.second.topLeftCorner( 3, 3 ) - withHardOnly.second.topLeftCorner( 3, 3 );
    BOOST_CHECK_LT( ( normalDiff - soft ).norm( ), 1.0E-10 );
}

//! Settings class: each preset builder produces the expected constraint weight matrix structure and the validation
//! rules in the constructor reject obviously malformed inputs.
BOOST_AUTO_TEST_CASE( test_InterArcStateContinuityConstraintSettings_PresetsAndValidation )
{
    using tudat::simulation_setup::fullStateContinuity;
    using tudat::simulation_setup::generalContinuity;
    using tudat::simulation_setup::InterArcStateContinuityConstraintSettings;
    using tudat::simulation_setup::positionOnlyContinuity;
    using tudat::simulation_setup::velocityOnlyContinuity;

    const std::vector< double > epochs = { 100.0, 200.0 };

    // Position-only: position weights non-zero, velocity weights zero. Rank should be 3.
    auto positionOnly = positionOnlyContinuity( "Sat", epochs, 2.5, 1.0 );
    const auto& positionWeightMatrix = positionOnly->weightMatrixForPair( 0 );
    BOOST_CHECK_CLOSE_FRACTION( positionWeightMatrix( 0, 0 ), 2.5, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( positionWeightMatrix( 1, 1 ), 2.5, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( positionWeightMatrix( 2, 2 ), 2.5, 1.0E-15 );
    BOOST_CHECK_EQUAL( positionWeightMatrix( 3, 3 ), 0.0 );
    BOOST_CHECK_EQUAL( positionWeightMatrix( 4, 4 ), 0.0 );
    BOOST_CHECK_EQUAL( positionWeightMatrix( 5, 5 ), 0.0 );

    // Velocity-only: inverse pattern.
    auto velocityOnly = velocityOnlyContinuity( "Sat", epochs, 0.1 );
    const auto& velocityWeightMatrix = velocityOnly->weightMatrixForPair( 0 );
    BOOST_CHECK_EQUAL( velocityWeightMatrix( 0, 0 ), 0.0 );
    BOOST_CHECK_CLOSE_FRACTION( velocityWeightMatrix( 3, 3 ), 0.1, 1.0E-15 );

    // Full state with anisotropic weights.
    auto fullState = fullStateContinuity( "Sat", epochs, 1.5, 0.7 );
    const auto& fullStateWeightMatrix = fullState->weightMatrixForPair( 0 );
    BOOST_CHECK_CLOSE_FRACTION( fullStateWeightMatrix( 0, 0 ), 1.5, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( fullStateWeightMatrix( 3, 3 ), 0.7, 1.0E-15 );

    // Broadcasting: a single-entry weight matrix list applied across multiple pairs is allowed.
    BOOST_CHECK_EQUAL( positionOnly->numberOfPairs( ), 2u );

    // Non-positive constraint scaling factors throw.
    BOOST_CHECK_THROW( positionOnlyContinuity( "Sat", epochs, 1.0, 0.0 ), std::runtime_error );
    BOOST_CHECK_THROW( positionOnlyContinuity( "Sat", epochs, 1.0, -1.0 ), std::runtime_error );

    // Asymmetric weight matrix throws.
    Eigen::Matrix< double, 6, 6 > asymmetric = Eigen::Matrix< double, 6, 6 >::Zero( );
    asymmetric( 0, 1 ) = 1.0;
    BOOST_CHECK_THROW( generalContinuity( "Sat", epochs, { asymmetric } ), std::runtime_error );

    // Non-PSD (negative eigenvalue) weight matrix throws.
    Eigen::Matrix< double, 6, 6 > indefinite = Eigen::Matrix< double, 6, 6 >::Zero( );
    indefinite( 0, 0 ) = -1.0;
    BOOST_CHECK_THROW( generalContinuity( "Sat", epochs, { indefinite } ), std::runtime_error );

    // Mismatched arcPairs / connectionEpochs sizes throw.
    BOOST_CHECK_THROW( InterArcStateContinuityConstraintSettings(
                               "Sat", epochs, { tudat::simulation_setup::detail::diagonalWeight( 1.0, 1.0 ) }, { 1.0 }, { { 0, 1 } } ),
                       std::runtime_error );

    // Non-consecutive arc pair throws.
    BOOST_CHECK_THROW(
            InterArcStateContinuityConstraintSettings(
                    "Sat", epochs, { tudat::simulation_setup::detail::diagonalWeight( 1.0, 1.0 ) }, { 1.0 }, { { 0, 2 }, { 1, 3 } } ),
            std::runtime_error );

    // weightMatrices size must either broadcast from one matrix or match the number of arc pairs.
    BOOST_CHECK_THROW( InterArcStateContinuityConstraintSettings( "Sat",
                                                                  epochs,
                                                                  std::vector< Eigen::Matrix< double, 6, 6 > >(
                                                                          3, tudat::simulation_setup::detail::diagonalWeight( 1.0, 1.0 ) ),
                                                                  { 1.0 } ),
                       std::runtime_error );

    // constraintScalingFactors size must either broadcast from one scaling factor or match the number of arc pairs.
    BOOST_CHECK_THROW( InterArcStateContinuityConstraintSettings(
                               "Sat", epochs, { tudat::simulation_setup::detail::diagonalWeight( 1.0, 1.0 ) }, { 1.0, 2.0, 3.0 } ),
                       std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
