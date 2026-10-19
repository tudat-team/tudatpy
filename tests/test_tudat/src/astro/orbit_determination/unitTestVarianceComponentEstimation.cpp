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

#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include <Eigen/Cholesky>
#include <Eigen/SparseCore>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/varianceComponentEstimation.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_variance_component_estimation )

using namespace tudat::estimatable_parameters;
using namespace tudat::numerical_integrators;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;

//! Small linear test problem: y = A x + noise, with a different noise level per group.
struct LinearProblem {
    Eigen::MatrixXd designMatrix;
    Eigen::VectorXd observations;
    std::vector< std::vector< int > > groupRowIndices;
    std::vector< double > trueSigmas;
};

LinearProblem createLinearProblem( const std::vector< int >& groupSizes,
                                   const std::vector< double >& trueSigmas,
                                   const int numberOfParameters )
{
    std::mt19937 generator( 42 );
    std::normal_distribution< double > normal( 0.0, 1.0 );

    LinearProblem problem;
    problem.trueSigmas = trueSigmas;
    int numberOfObservations = 0;
    for( int size : groupSizes )
    {
        numberOfObservations += size;
    }

    problem.designMatrix.resize( numberOfObservations, numberOfParameters );
    for( int i = 0; i < numberOfObservations; i++ )
    {
        for( int j = 0; j < numberOfParameters; j++ )
        {
            problem.designMatrix( i, j ) = normal( generator );
        }
    }
    Eigen::VectorXd trueParameters = Eigen::VectorXd::LinSpaced( numberOfParameters, 1.0, 2.0 );
    problem.observations = problem.designMatrix * trueParameters;

    int row = 0;
    for( unsigned int k = 0; k < groupSizes.size( ); k++ )
    {
        problem.groupRowIndices.push_back( std::vector< int >( ) );
        for( int i = 0; i < groupSizes.at( k ); i++ )
        {
            problem.observations( row ) += trueSigmas.at( k ) * normal( generator );
            problem.groupRowIndices.back( ).push_back( row );
            row++;
        }
    }
    return problem;
}

//! Weighted least-squares solution with an optional prior (centred on zero). Returns the residuals and sets the normal matrix.
Eigen::VectorXd solveLinearProblem( const LinearProblem& problem,
                                    const Eigen::VectorXd& weights,
                                    const Eigen::MatrixXd& inverseAprioriCovariance,
                                    Eigen::MatrixXd& normalMatrix )
{
    const Eigen::MatrixXd& designMatrix = problem.designMatrix;
    normalMatrix = designMatrix.transpose( ) * weights.asDiagonal( ) * designMatrix + inverseAprioriCovariance;
    Eigen::VectorXd estimate = normalMatrix.ldlt( ).solve( designMatrix.transpose( ) * weights.cwiseProduct( problem.observations ) );
    return problem.observations - designMatrix * estimate;
}

//! Iterate variance component estimation on a linear problem until the factors are 1. Returns the final weights.
Eigen::VectorXd iterateLinearProblem( const LinearProblem& problem,
                                      Eigen::VectorXd weights,
                                      const Eigen::MatrixXd& inverseAprioriCovariance,
                                      const bool useSimplifiedEstimator,
                                      int& numberOfIterations )
{
    Eigen::MatrixXd normalMatrix;
    for( numberOfIterations = 1; numberOfIterations <= 200; numberOfIterations++ )
    {
        Eigen::VectorXd residuals = solveLinearProblem( problem, weights, inverseAprioriCovariance, normalMatrix );
        VarianceComponentEstimate estimate =
                computeHelmertVarianceComponents( problem.designMatrix, weights, normalMatrix, residuals, problem.groupRowIndices );

        double maximumChange = 0.0;
        for( unsigned int k = 0; k < problem.groupRowIndices.size( ); k++ )
        {
            const double factor = useSimplifiedEstimator ? estimate.groupResults.at( k ).simplifiedVarianceFactor
                                                         : estimate.groupResults.at( k ).varianceFactor;
            maximumChange = std::max( maximumChange, std::fabs( factor - 1.0 ) );
            for( int row : problem.groupRowIndices.at( k ) )
            {
                weights( row ) /= factor;
            }
        }
        if( maximumChange < 1.0E-12 )
        {
            break;
        }
    }
    return weights;
}

//! Weight vector with a constant sigma per group, and the given weight for rows in no group.
Eigen::VectorXd getGroupWeights( const LinearProblem& problem, const std::vector< double >& sigmas, const double ungroupedWeight = 1.0 )
{
    Eigen::VectorXd weights = Eigen::VectorXd::Constant( problem.designMatrix.rows( ), ungroupedWeight );
    for( unsigned int k = 0; k < problem.groupRowIndices.size( ); k++ )
    {
        for( int row : problem.groupRowIndices.at( k ) )
        {
            weights( row ) = 1.0 / ( sigmas.at( k ) * sigmas.at( k ) );
        }
    }
    return weights;
}

//! The redundancy of each group must equal the sum of the diagonal of I - W A N^-1 A^T over its rows, and the
//! redundancies must add up to m - p (no prior) or m - p + tr(N^-1 P) (with prior).
BOOST_AUTO_TEST_CASE( test_RedundancyNumbers )
{
    const int numberOfParameters = 6;
    LinearProblem problem = createLinearProblem( { 300, 100 }, { 0.5, 20.0 }, numberOfParameters );
    const Eigen::VectorXd weights = getGroupWeights( problem, { 2.0, 3.0 } );
    const int numberOfObservations = problem.designMatrix.rows( );

    for( int usePrior = 0; usePrior < 2; usePrior++ )
    {
        Eigen::MatrixXd inverseAprioriCovariance = Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters );
        if( usePrior )
        {
            inverseAprioriCovariance.diagonal( ) = Eigen::VectorXd::LinSpaced( numberOfParameters, 1.0, 50.0 );
        }

        Eigen::MatrixXd normalMatrix;
        Eigen::VectorXd residuals = solveLinearProblem( problem, weights, inverseAprioriCovariance, normalMatrix );
        VarianceComponentEstimate estimate =
                computeHelmertVarianceComponents( problem.designMatrix, weights, normalMatrix, residuals, problem.groupRowIndices );

        // Explicit redundancy numbers from the hat matrix.
        const Eigen::MatrixXd hatMatrix =
                weights.asDiagonal( ) * problem.designMatrix * normalMatrix.inverse( ) * problem.designMatrix.transpose( );
        const Eigen::VectorXd redundancyNumbers = Eigen::VectorXd::Ones( numberOfObservations ) - hatMatrix.diagonal( );

        double totalRedundancy = 0.0;
        for( unsigned int k = 0; k < problem.groupRowIndices.size( ); k++ )
        {
            double expectedRedundancy = 0.0;
            double expectedSquareSum = 0.0;
            for( int row : problem.groupRowIndices.at( k ) )
            {
                expectedRedundancy += redundancyNumbers( row );
                expectedSquareSum += weights( row ) * residuals( row ) * residuals( row );
            }
            BOOST_CHECK_CLOSE_FRACTION( estimate.groupResults.at( k ).redundancy, expectedRedundancy, 1.0E-10 );
            BOOST_CHECK_CLOSE_FRACTION( estimate.groupResults.at( k ).weightedSquareSum, expectedSquareSum, 1.0E-12 );
            BOOST_CHECK_EQUAL( estimate.groupResults.at( k ).numberOfObservations,
                               static_cast< int >( problem.groupRowIndices.at( k ).size( ) ) );
            totalRedundancy += estimate.groupResults.at( k ).redundancy;
        }

        const double priorShare = ( normalMatrix.inverse( ) * inverseAprioriCovariance ).trace( );
        BOOST_CHECK_CLOSE_FRACTION( totalRedundancy, numberOfObservations - numberOfParameters + priorShare, 1.0E-10 );
    }
}

//! Scaling the parameters (columns of the design matrix) must not change any result.
BOOST_AUTO_TEST_CASE( test_NormalizationInvariance )
{
    const int numberOfParameters = 6;
    LinearProblem problem = createLinearProblem( { 300, 100 }, { 0.5, 20.0 }, numberOfParameters );
    const Eigen::VectorXd weights = getGroupWeights( problem, { 2.0, 3.0 } );
    Eigen::MatrixXd inverseAprioriCovariance = Eigen::MatrixXd::Identity( numberOfParameters, numberOfParameters );

    Eigen::MatrixXd normalMatrix;
    Eigen::VectorXd residuals = solveLinearProblem( problem, weights, inverseAprioriCovariance, normalMatrix );
    VarianceComponentEstimate estimate =
            computeHelmertVarianceComponents( problem.designMatrix, weights, normalMatrix, residuals, problem.groupRowIndices );

    const Eigen::VectorXd scaling = Eigen::VectorXd::LinSpaced( numberOfParameters, 1.0E-4, 1.0E4 );
    const Eigen::MatrixXd scaledDesignMatrix = problem.designMatrix * scaling.asDiagonal( );
    const Eigen::MatrixXd scaledNormalMatrix = scaling.asDiagonal( ) * normalMatrix * scaling.asDiagonal( );
    VarianceComponentEstimate scaledEstimate =
            computeHelmertVarianceComponents( scaledDesignMatrix, weights, scaledNormalMatrix, residuals, problem.groupRowIndices );

    for( unsigned int k = 0; k < problem.groupRowIndices.size( ); k++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( scaledEstimate.groupResults.at( k ).redundancy, estimate.groupResults.at( k ).redundancy, 1.0E-9 );
        BOOST_CHECK_CLOSE_FRACTION(
                scaledEstimate.groupResults.at( k ).varianceFactor, estimate.groupResults.at( k ).varianceFactor, 1.0E-9 );
    }
}

//! Starting from weights that are wrong by a factor 10 in sigma, the iteration must recover the noise levels. Helmert's
//! method and the simplified estimator must converge to the same weights, with and without a prior, and with
//! observations that are in no group.
BOOST_AUTO_TEST_CASE( test_ConvergenceToTrueSigmas )
{
    const int numberOfParameters = 6;
    const std::vector< double > trueSigmas = { 0.5, 20.0 };
    LinearProblem problem = createLinearProblem( { 400, 100 }, trueSigmas, numberOfParameters );

    for( int testCase = 0; testCase < 3; testCase++ )
    {
        LinearProblem currentProblem = problem;
        Eigen::MatrixXd inverseAprioriCovariance = Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters );
        if( testCase == 1 )
        {
            inverseAprioriCovariance.diagonal( ) = Eigen::VectorXd::Constant( numberOfParameters, 10.0 );
        }
        else if( testCase == 2 )
        {
            // Last 50 rows of the second group are in no group, and keep their weight.
            currentProblem.groupRowIndices.at( 1 ).resize( 50 );
        }

        const Eigen::VectorXd initialWeights = getGroupWeights( currentProblem, { 5.0, 2.0 }, 1.0 / ( 20.0 * 20.0 ) );
        int helmertIterations = 0, simplifiedIterations = 0;
        const Eigen::VectorXd helmertWeights =
                iterateLinearProblem( currentProblem, initialWeights, inverseAprioriCovariance, false, helmertIterations );
        const Eigen::VectorXd simplifiedWeights =
                iterateLinearProblem( currentProblem, initialWeights, inverseAprioriCovariance, true, simplifiedIterations );

        BOOST_CHECK_LT( helmertIterations, 200 );
        BOOST_CHECK_LT( simplifiedIterations, 200 );
        for( unsigned int k = 0; k < currentProblem.groupRowIndices.size( ); k++ )
        {
            const int row = currentProblem.groupRowIndices.at( k ).front( );
            BOOST_CHECK_CLOSE_FRACTION( helmertWeights( row ), simplifiedWeights( row ), 1.0E-8 );

            const double estimatedSigma = 1.0 / std::sqrt( helmertWeights( row ) );
            BOOST_CHECK_CLOSE_FRACTION( estimatedSigma, trueSigmas.at( k ), 0.15 );
        }
    }
}

//! Group 0 has noise that is correlated in pairs of observations, with block weights that match it. The iteration must
//! recover both noise levels, with the same result for Helmert's method and the simplified estimator. A correlation between
//! two groups must be rejected.
BOOST_AUTO_TEST_CASE( test_BlockWeights )
{
    const int numberOfParameters = 6;
    const std::vector< double > trueSigmas = { 0.5, 20.0 };
    const double correlation = 0.8;
    LinearProblem problem = createLinearProblem( { 400, 100 }, { 0.0, trueSigmas.at( 1 ) }, numberOfParameters );

    // Correlated noise in pairs for group 0: e = sigma * L z, with L L^T = R.
    Eigen::Matrix2d correlationMatrix;
    correlationMatrix << 1.0, correlation, correlation, 1.0;
    const Eigen::Matrix2d choleskyFactor = correlationMatrix.llt( ).matrixL( );
    std::mt19937 generator( 7 );
    std::normal_distribution< double > normal( 0.0, 1.0 );
    for( int pair = 0; pair < 200; pair++ )
    {
        const Eigen::Vector2d noise = trueSigmas.at( 0 ) * choleskyFactor * Eigen::Vector2d( normal( generator ), normal( generator ) );
        problem.observations.segment( 2 * pair, 2 ) += noise;
    }

    // Start 10x off in sigma: group 0 weight blocks R^-1 / 5^2, group 1 weights 1 / 2^2.
    const Eigen::Matrix2d inverseCorrelation = correlationMatrix.inverse( );
    auto createWeights = [ & ]( const double sigma0, const double sigma1 ) {
        std::vector< Eigen::Triplet< double > > entries;
        for( int pair = 0; pair < 200; pair++ )
        {
            for( int i = 0; i < 2; i++ )
            {
                for( int j = 0; j < 2; j++ )
                {
                    entries.emplace_back( 2 * pair + i, 2 * pair + j, inverseCorrelation( i, j ) / ( sigma0 * sigma0 ) );
                }
            }
        }
        for( int row = 400; row < 500; row++ )
        {
            entries.emplace_back( row, row, 1.0 / ( sigma1 * sigma1 ) );
        }
        Eigen::SparseMatrix< double > weights( 500, 500 );
        weights.setFromTriplets( entries.begin( ), entries.end( ) );
        return weights;
    };

    for( int useSimplified = 0; useSimplified < 2; useSimplified++ )
    {
        Eigen::SparseMatrix< double > weights = createWeights( 5.0, 2.0 );
        Eigen::VectorXd totalFactors = Eigen::VectorXd::Ones( 2 );
        for( int iteration = 0; iteration < 200; iteration++ )
        {
            const Eigen::MatrixXd normalMatrix = problem.designMatrix.transpose( ) * ( weights * problem.designMatrix );
            const Eigen::VectorXd estimate =
                    normalMatrix.ldlt( ).solve( problem.designMatrix.transpose( ) * ( weights * problem.observations ) );
            const Eigen::VectorXd residuals = problem.observations - problem.designMatrix * estimate;
            const Eigen::VectorXd diagonal = weights.diagonal( );
            VarianceComponentEstimate result = computeHelmertVarianceComponents(
                    problem.designMatrix, diagonal, weights, normalMatrix, residuals, problem.groupRowIndices );

            double maximumChange = 0.0;
            for( int k = 0; k < 2; k++ )
            {
                const double factor =
                        useSimplified ? result.groupResults.at( k ).simplifiedVarianceFactor : result.groupResults.at( k ).varianceFactor;
                totalFactors( k ) *= factor;
                maximumChange = std::max( maximumChange, std::fabs( factor - 1.0 ) );
            }
            weights = createWeights( 5.0 * std::sqrt( totalFactors( 0 ) ), 2.0 * std::sqrt( totalFactors( 1 ) ) );
            if( maximumChange < 1.0E-12 )
            {
                break;
            }
        }
        BOOST_CHECK_CLOSE_FRACTION( 5.0 * std::sqrt( totalFactors( 0 ) ), trueSigmas.at( 0 ), 0.15 );
        BOOST_CHECK_CLOSE_FRACTION( 2.0 * std::sqrt( totalFactors( 1 ) ), trueSigmas.at( 1 ), 0.15 );
    }

    // A correlation between the two groups is rejected.
    Eigen::SparseMatrix< double > crossGroupWeights = createWeights( 1.0, 1.0 );
    crossGroupWeights.coeffRef( 0, 450 ) = 0.1;
    crossGroupWeights.coeffRef( 450, 0 ) = 0.1;
    const Eigen::MatrixXd normalMatrix = problem.designMatrix.transpose( ) * ( crossGroupWeights * problem.designMatrix );
    const Eigen::VectorXd diagonal = crossGroupWeights.diagonal( );
    BOOST_CHECK_THROW( computeHelmertVarianceComponents( problem.designMatrix,
                                                         diagonal,
                                                         crossGroupWeights,
                                                         normalMatrix,
                                                         Eigen::VectorXd::Zero( 500 ),
                                                         problem.groupRowIndices ),
                       std::runtime_error );
}

//! A group with too little redundancy is kept fixed. Overlapping groups and bad input are rejected.
BOOST_AUTO_TEST_CASE( test_GroupChecks )
{
    const int numberOfParameters = 6;
    LinearProblem problem = createLinearProblem( { 300, 3 }, { 0.5, 20.0 }, numberOfParameters );
    const Eigen::VectorXd weights = getGroupWeights( problem, { 0.5, 20.0 } );
    Eigen::MatrixXd normalMatrix;
    Eigen::VectorXd residuals =
            solveLinearProblem( problem, weights, Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters ), normalMatrix );

    VarianceComponentEstimate estimate =
            computeHelmertVarianceComponents( problem.designMatrix, weights, normalMatrix, residuals, problem.groupRowIndices, 5.0 );
    BOOST_CHECK( !estimate.groupResults.at( 0 ).isFrozen );
    BOOST_CHECK( estimate.groupResults.at( 1 ).isFrozen );
    BOOST_CHECK_EQUAL( estimate.groupResults.at( 1 ).varianceFactor, 1.0 );
    BOOST_CHECK_EQUAL( estimate.helmertMatrix.rows( ), 1 );

    // Rows with zero weight are not counted.
    Eigen::VectorXd weightsWithZeros = weights;
    weightsWithZeros.segment( 0, 10 ).setZero( );
    Eigen::VectorXd residualsWithZeros =
            solveLinearProblem( problem, weightsWithZeros, Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters ), normalMatrix );
    estimate = computeHelmertVarianceComponents(
            problem.designMatrix, weightsWithZeros, normalMatrix, residualsWithZeros, problem.groupRowIndices );
    BOOST_CHECK_EQUAL( estimate.groupResults.at( 0 ).numberOfObservations, 290 );

    std::vector< std::vector< int > > overlappingGroups = problem.groupRowIndices;
    overlappingGroups.at( 1 ).push_back( 0 );
    BOOST_CHECK_THROW( computeHelmertVarianceComponents( problem.designMatrix, weights, normalMatrix, residuals, overlappingGroups ),
                       std::runtime_error );

    std::vector< std::vector< int > > outOfRangeGroups = { { 0, 1, 1000 } };
    BOOST_CHECK_THROW( computeHelmertVarianceComponents( problem.designMatrix, weights, normalMatrix, residuals, outOfRangeGroups ),
                       std::runtime_error );
}

//! Orbit estimation with range and position observations, in one or two arcs. Each arc has its own range noise level.
//! The variance factors are estimated per observable type (and per arc, for the multi-arc case).
void runOrbitVarianceComponentEstimation( const bool useMultiArc )
{
    const double gravitationalParameter = 3.986004418E14;
    const double arcDuration = 4.0 * 3600.0;
    const int numberOfArcs = useMultiArc ? 2 : 1;
    std::vector< double > arcStartTimes, arcEndTimes;
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        arcStartTimes.push_back( arc * arcDuration );
        arcEndTimes.push_back( ( arc + 1 ) * arcDuration );
    }

    // Earth with a fixed position and orientation, and a vehicle whose orbit is estimated.
    BodyListSettings bodySettings( "SSB", "J2000" );
    bodySettings.addSettings( "Earth" );
    bodySettings.at( "Earth" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings.at( "Earth" )->gravityFieldSettings = centralGravitySettings( gravitationalParameter );
    bodySettings.at( "Earth" )->shapeModelSettings = sphericalBodyShapeSettings( 6.371E6 );
    bodySettings.at( "Earth" )->rotationModelSettings =
            constantRotationModelSettings( "J2000", "IAU_Earth", Eigen::Quaterniond::Identity( ) );
    for( int i = 0; i < 2; i++ )
    {
        bodySettings.at( "Earth" )->groundStationSettings.push_back( groundStationSettings( "Station" + std::to_string( i ),
                                                                                            Eigen::Vector3d( 0.0, 0.3 * i, 1.0 + 2.0 * i ),
                                                                                            coordinate_conversions::geodetic_position ) );
    }
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->ephemerisSettings =
            std::make_shared< TabulatedEphemerisSettings >( std::map< double, Eigen::Vector6d >( ), "Earth", "J2000" );
    if( useMultiArc )
    {
        bodySettings.at( "Vehicle" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    }
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    const std::vector< std::string > bodiesToIntegrate = { "Vehicle" };
    const std::vector< std::string > centralBodies = { "Earth" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    Eigen::Vector6d keplerianElements;
    keplerianElements << 7200.0E3, 0.05, 1.2, 0.5, 0.3, 0.0;
    const Eigen::Vector6d initialState =
            orbital_element_conversions::convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKutta4Settings< double >( 10.0 );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > arcPropagatorSettings;
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        // Each arc starts from a slightly different state, so the arcs are independent.
        arcPropagatorSettings.push_back( std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                centralBodies,
                accelerationModelMap,
                bodiesToIntegrate,
                initialState * ( 1.0 + 0.01 * arc ),
                arcStartTimes.at( arc ),
                integratorSettings,
                propagationTimeTerminationSettings( arcEndTimes.at( arc ) ) ) );
    }
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings;
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    if( useMultiArc )
    {
        auto multiArcSettings = std::make_shared< MultiArcPropagatorSettings< double, double > >( arcPropagatorSettings );
        parameterNames = getInitialMultiArcParameterSettings< double, double >( multiArcSettings, bodies, arcStartTimes );
        propagatorSettings = multiArcSettings;
    }
    else
    {
        parameterNames = getInitialStateParameterSettings< double, double >( arcPropagatorSettings.at( 0 ), bodies );
        propagatorSettings = arcPropagatorSettings.at( 0 );
    }
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies );

    // Range from two stations and the position of the vehicle.
    std::vector< LinkEnds > rangeLinkEnds;
    for( int i = 0; i < 2; i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Earth", "Station" + std::to_string( i ) );
        linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
        rangeLinkEnds.push_back( linkEnds );
    }
    LinkEnds positionLinkEnds;
    positionLinkEnds[ observed_body ] = LinkEndId( "Vehicle", "" );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( const LinkEnds& linkEnds : rangeLinkEnds )
    {
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, linkEnds ) );
    }
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, positionLinkEnds ) );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettingsList, propagatorSettings );
    const Eigen::VectorXd trueParameters = parametersToEstimate->getFullParameterValues< double >( );

    // True noise: range sigma differs per arc, position sigma is the same for all arcs.
    const std::vector< double > trueRangeSigmas = { 1.0, 3.0 };
    const double truePositionSigma = 50.0;
    Eigen::Matrix3d positionCorrelation;
    positionCorrelation << 1.0, 0.5, 0.0, 0.5, 1.0, 0.3, 0.0, 0.3, 1.0;
    const Eigen::Matrix3d positionCholesky = positionCorrelation.llt( ).matrixL( );
    auto generator = std::make_shared< std::mt19937 >( 1234 );
    auto noiseFunction = [ generator ]( const double sigma, const int size ) {
        return std::function< Eigen::VectorXd( const double ) >( [ generator, sigma, size ]( const double ) {
            std::normal_distribution< double > normal( 0.0, sigma );
            Eigen::VectorXd noise( size );
            for( int i = 0; i < size; i++ )
            {
                noise( i ) = normal( *generator );
            }
            return noise;
        } );
    };

    // One set of simulation settings per arc, so that each observation set lies within one arc.
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > simulationSettings;
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        std::vector< double > rangeTimes, positionTimes;
        for( double time = arcStartTimes.at( arc ) + 60.0; time < arcEndTimes.at( arc ) - 60.0; time += 30.0 )
        {
            rangeTimes.push_back( time );
        }
        for( double time = arcStartTimes.at( arc ) + 60.0; time < arcEndTimes.at( arc ) - 60.0; time += 600.0 )
        {
            positionTimes.push_back( time );
        }
        for( const LinkEnds& linkEnds : rangeLinkEnds )
        {
            auto settings =
                    std::make_shared< TabulatedObservationSimulationSettings< double > >( one_way_range, linkEnds, rangeTimes, receiver );
            settings->setObservationNoiseFunction( noiseFunction( trueRangeSigmas.at( arc ), 1 ) );
            simulationSettings.push_back( settings );
        }
        auto settings = std::make_shared< TabulatedObservationSimulationSettings< double > >(
                position_observable, positionLinkEnds, positionTimes, observed_body );
        const auto independentNoise = noiseFunction( truePositionSigma, 3 );
        settings->setObservationNoiseFunction( std::function< Eigen::VectorXd( const double ) >(
                [ independentNoise, positionCholesky ]( const double time ) -> Eigen::VectorXd {
                    return positionCholesky * independentNoise( time );
                } ) );
        simulationSettings.push_back( settings );
    }
    std::shared_ptr< ObservationDataset< double, double > > observations = simulateObservationDataset< double, double >(
            simulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Groups: observable type, per arc.
    using Condition = ObservationSelectionCondition< double, double >;
    std::vector< Condition > groups;
    std::vector< double > trueSigmas, initialSigmas;
    for( ObservableType observableType : { one_way_range, position_observable } )
    {
        for( int arc = 0; arc < numberOfArcs; arc++ )
        {
            groups.push_back( Condition::observableType( observableType ) &&
                              Condition::timeBounds( arcStartTimes.at( arc ), arcEndTimes.at( arc ) ) );
            trueSigmas.push_back( observableType == one_way_range ? trueRangeSigmas.at( arc ) : truePositionSigma );

            // Start a factor 10 off: range weights too low, position weights too high.
            initialSigmas.push_back( observableType == one_way_range ? 10.0 : 5.0 );
            if( observableType == one_way_range )
            {
                observations->setConstantSingleObservationScalarWeight( groups.back( ), 1.0 / ( 10.0 * 10.0 ) );
            }
            else
            {
                // Correlated position components, with a weight matrix that matches the noise correlation.
                observations->setConstantSingleObservationMatrixWeight( groups.back( ), positionCorrelation.inverse( ) / ( 5.0 * 5.0 ) );
            }
        }
    }

    // Perturb the initial state of each arc.
    Eigen::VectorXd perturbedParameters = trueParameters;
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        perturbedParameters.segment( 6 * arc, 3 ) += Eigen::Vector3d::Constant( 100.0 );
        perturbedParameters.segment( 6 * arc + 3, 3 ) += Eigen::Vector3d::Constant( 0.1 );
    }
    parametersToEstimate->resetParameterValues( perturbedParameters );

    auto estimationInput = std::make_shared< EstimationInput< double, double > >(
            observations, Eigen::MatrixXd::Zero( 0, 0 ), estimationConvergenceChecker( 6 ) );
    estimationInput->defineEstimationSettings( true, true, true, false, false, false );

    auto settings = varianceComponentSettings< double, double >( groups, 0.01, 15 );
    auto output = estimateParametersWithVarianceComponents( orbitDeterminationManager, estimationInput, settings );

    BOOST_CHECK( output->isConverged );
    BOOST_CHECK_LT( output->iterationResults.size( ), 15 );
    BOOST_CHECK_EQUAL( output->iterationResults.size( ), output->totalVarianceFactorHistory.size( ) );

    const Eigen::VectorXd totalFactors = output->getTotalVarianceFactors( );
    for( unsigned int k = 0; k < groups.size( ); k++ )
    {
        const double estimatedSigma = initialSigmas.at( k ) * std::sqrt( totalFactors( k ) );
        const bool isRangeGroup = ( k < static_cast< unsigned int >( numberOfArcs ) );
        BOOST_CHECK_CLOSE_FRACTION( estimatedSigma, trueSigmas.at( k ), isRangeGroup ? 0.08 : 0.3 );

        // The weights in the observation dataset are the final, estimated weights. Position weights keep their correlation.
        const Eigen::VectorXd groupWeights = observations->getWeightDiagonal( groups.at( k ) );
        const double expectedWeight = ( isRangeGroup ? 1.0 : positionCorrelation.inverse( )( 0, 0 ) ) / ( estimatedSigma * estimatedSigma );
        BOOST_CHECK_CLOSE_FRACTION( groupWeights( 0 ), expectedWeight, 1.0E-12 );
    }

    // The estimated orbit is close to the truth (within a few range sigmas).
    const Eigen::VectorXd parameterError = output->estimationOutput->parameterEstimate_ - trueParameters;
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        BOOST_CHECK_SMALL( parameterError.segment( 6 * arc, 3 ).norm( ), 10.0 );
    }
}

BOOST_AUTO_TEST_CASE( test_SingleArcOrbitVarianceComponents )
{
    runOrbitVarianceComponentEstimation( false );
}

BOOST_AUTO_TEST_CASE( test_MultiArcOrbitVarianceComponents )
{
    runOrbitVarianceComponentEstimation( true );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
