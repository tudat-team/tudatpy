/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INTERARCCONTINUITYCONSTRAINT_H
#define TUDAT_INTERARCCONTINUITYCONSTRAINT_H

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/propagators/stateTransitionMatrixInterface.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/simulation/estimation_setup/interArcStateContinuityConstraintSettings.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolverBase.h"
#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"

namespace tudat
{
namespace simulation_setup
{

//! Aggregated normal-equation contribution of all soft inter-arc continuity priors.
struct InterArcConstraintContribution {
    Eigen::MatrixXd additionalNormalMatrix;
    Eigen::VectorXd additionalRightHandSide;
    double totalConstraintCost = 0.0;
    std::vector< Eigen::VectorXd > perPairDiscrepancies;
};

namespace detail
{

//! Numerical rank of a PSD constraint weight matrix using its eigenvalue spectrum.
inline int computePsdMatrixRank( const Eigen::MatrixXd& constraintWeightMatrix )
{
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( 0.5 * ( constraintWeightMatrix + constraintWeightMatrix.transpose( ) ) );
    if( solver.info( ) != Eigen::Success )
    {
        throw std::runtime_error( "Error in computePsdMatrixRank: eigen-decomposition failed." );
    }
    const double largest = std::max( solver.eigenvalues( ).maxCoeff( ), 0.0 );
    if( largest <= 0.0 )
    {
        return 0;
    }
    const double threshold = 1.0E-12 * largest;
    int rank = 0;
    for( int i = 0; i < solver.eigenvalues( ).size( ); ++i )
    {
        if( solver.eigenvalues( )( i ) > threshold )
        {
            ++rank;
        }
    }
    return rank;
}

//! Single-arc state evaluator at an arbitrary time inside [arcInitialTime, arcFinalTime].
//! Short-circuits at the arc endpoints (avoids interpolation noise at shared OCM boundaries) and otherwise uses
//! an 8th-order Lagrange interpolation through a local neighbourhood of the propagated solution, mirroring
//! getArcInitialStateFromPreviousArcResult in multiArcDynamicsSimulator.h.
//!
//! Adaptive-order fallback: when fewer than 8 samples are available on one side of the evaluation time (i.e.
//! evaluationTime is closer than ~4 integrator steps to either arc edge), the interpolation order silently
//! drops to the largest even value supported by the local neighbourhood — down to 2 (linear) in the limit.
//! This is normally invisible because the boundary short-circuits catch arc-start/arc-end evaluations, but
//! callers using mid-gap connection epochs near an arc boundary should be aware that accuracy degrades there.
template< typename StateScalarType, typename TimeType >
Eigen::VectorXd evaluateArcStateAtTime( const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& arcSolution,
                                        const double evaluationTime,
                                        const double arcInitialTime,
                                        const double arcFinalTime,
                                        const int stateOffset,
                                        const int stateBlockSize )
{
    if( arcSolution.empty( ) )
    {
        throw std::runtime_error( "Error in evaluateArcStateAtTime: arc solution is empty." );
    }

    auto extract = [ stateOffset, stateBlockSize ]( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& fullState ) {
        Eigen::VectorXd result( stateBlockSize );
        for( int i = 0; i < stateBlockSize; ++i )
        {
            result( i ) = static_cast< double >( fullState( stateOffset + i ) );
        }
        return result;
    };

    // Use a tolerance proportional to the arc span; the fixed-step integrator can stop a few wall-clock units
    // short of the nominal arc end time. Also short-circuit when the evaluation time is at or beyond either
    // actual extremum of the sampled solution (Lagrange interpolation isn't well-posed at the bounds anyway).
    const double arcSpan = std::max( arcFinalTime - arcInitialTime, 1.0 );
    const double boundaryTolerance = std::max( 1.0E-7 * arcSpan, 1.0E-6 );
    const double firstSampleTime = static_cast< double >( arcSolution.begin( )->first );
    const double lastSampleTime = static_cast< double >( arcSolution.rbegin( )->first );

    if( evaluationTime <= firstSampleTime || std::fabs( evaluationTime - arcInitialTime ) < boundaryTolerance ||
        std::fabs( evaluationTime - firstSampleTime ) < boundaryTolerance )
    {
        return extract( arcSolution.begin( )->second );
    }
    if( evaluationTime >= lastSampleTime || std::fabs( evaluationTime - arcFinalTime ) < boundaryTolerance ||
        std::fabs( evaluationTime - lastSampleTime ) < boundaryTolerance )
    {
        return extract( arcSolution.rbegin( )->second );
    }

    // Find the local window for 8th-order Lagrange interpolation. The earlier short-circuits guarantee that
    // upper_bound lies strictly inside the sampled solution, so it is neither begin() nor end().
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > interpolationMap;
    auto upper = arcSolution.upper_bound( static_cast< TimeType >( evaluationTime ) );
    auto lower = std::prev( upper );

    // Take ~4 samples on each side of the evaluation time (Lagrange order 8).
    const int half = 4;
    auto leftStart = lower;
    for( int i = 0; i < half - 1; ++i )
    {
        if( leftStart == arcSolution.begin( ) ) break;
        --leftStart;
    }
    auto rightEnd = upper;
    for( int i = 0; i < half - 1; ++i )
    {
        if( std::next( rightEnd ) == arcSolution.end( ) ) break;
        ++rightEnd;
    }

    auto it = leftStart;
    while( true )
    {
        interpolationMap[ it->first ] = it->second;
        if( it == rightEnd )
        {
            break;
        }
        ++it;
    }

    // Drop the interpolation order if not enough surrounding samples are available. tudat's LagrangeInterpolator
    // requires (a) an even number of stages and (b) dataMap.size() >= numberOfStages — so we round the available
    // sample count down to the nearest even value, with a minimum of 2 (linear).
    const int availableSamples = static_cast< int >( interpolationMap.size( ) );
    if( availableSamples < 2 )
    {
        const double leftGap = std::fabs( evaluationTime - firstSampleTime );
        const double rightGap = std::fabs( evaluationTime - lastSampleTime );
        return extract( leftGap < rightGap ? arcSolution.begin( )->second : arcSolution.rbegin( )->second );
    }
    int order = std::min( 8, availableSamples );
    if( order % 2 != 0 )
    {
        --order;
    }
    if( order < 2 )
    {
        order = 2;
    }
    auto interpolator = std::make_shared<
            interpolators::LagrangeInterpolator< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >, long double > >(
            interpolationMap, order );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > interpolated =
            interpolator->interpolate( static_cast< TimeType >( evaluationTime ) );
    return extract( interpolated );
}

//! Arc-wise state and full-STM row slices for a body in an explicit arc of the multi-arc layout.
struct BodyArcRows {
    int arcWiseStateStart = -1;
    int fullStateStart = -1;
    int size = 0;
};

template< typename StateScalarType >
BodyArcRows getBodyRowsInMultiArcLayout(
        const std::shared_ptr< propagators::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >& stmInterface,
        const int arcIndex,
        const std::string& body )
{
    const auto indexMapPerArc = stmInterface->getArcWiseAndFullSolutionInitialStateIndices( );
    const auto arcIterator = indexMapPerArc.find( arcIndex );
    if( arcIterator == indexMapPerArc.end( ) )
    {
        throw std::runtime_error( "Error in assembleInterArcContinuityContribution: no STM layout information found for arc " +
                                  std::to_string( arcIndex ) + "." );
    }
    const auto bodyIterator = arcIterator->second.find( body );
    if( bodyIterator == arcIterator->second.end( ) )
    {
        throw std::runtime_error( "Error in assembleInterArcContinuityContribution: body \"" + body + "\" is not present in arc " +
                                  std::to_string( arcIndex ) + " of the multi-arc STM layout." );
    }
    BodyArcRows rows;
    rows.arcWiseStateStart = bodyIterator->second.first.first;
    rows.fullStateStart = bodyIterator->second.second.first.first;
    rows.size = bodyIterator->second.second.second;
    return rows;
}

}  // namespace detail

//! Assemble the normal-equation contribution of all inter-arc continuity priors at the current linearization point.
//!
//! For each constrained boundary, this evaluates the left- and right-arc propagated states at the connection
//! epoch, forms their state discrepancy, extracts the corresponding rows of the full variational matrices, and
//! adds the resulting soft prior to the normal equations. Per Lari et al. (2021) Eq. 28, this is:
//!   stateDiscrepancy = rightArcState(connectionEpoch) - leftArcState(connectionEpoch)
//!   continuityDesignMatrix = rightArcVariationalRows - leftArcVariationalRows
//!   constraintWeight = constraintWeightMatrix / (constraintScalingFactor * totalConstraintDimension)
//!   normal += continuityDesignMatrix^T * constraintWeight * continuityDesignMatrix
//!   rhs    -= continuityDesignMatrix^T * constraintWeight * stateDiscrepancy
//! where totalConstraintDimension is the rank sum across every settings entry passed in. Column-normalization of
//! continuityDesignMatrix uses the same factors that the OD loop applies to the observation design matrix.
template< typename ObservationScalarType, typename TimeType >
InterArcConstraintContribution assembleInterArcContinuityContribution(
        const std::vector< std::shared_ptr< InterArcStateContinuityConstraintSettings > >& constraintSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >& parametersToEstimate,
        const std::shared_ptr< propagators::MultiArcDynamicsSimulator< ObservationScalarType, TimeType > >& multiArcSimulator,
        const std::shared_ptr< propagators::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< ObservationScalarType > >&
                stmInterface,
        const Eigen::VectorXd& columnNormalizationFactors,
        const int totalParameterSize )
{
    InterArcConstraintContribution contribution;
    contribution.additionalNormalMatrix = Eigen::MatrixXd::Zero( totalParameterSize, totalParameterSize );
    contribution.additionalRightHandSide = Eigen::VectorXd::Zero( totalParameterSize );
    contribution.totalConstraintCost = 0.0;
    contribution.perPairDiscrepancies.clear( );

    if( constraintSettings.empty( ) )
    {
        return contribution;
    }
    if( stmInterface == nullptr || multiArcSimulator == nullptr )
    {
        throw std::runtime_error(
                "Error in assembleInterArcContinuityContribution: multi-arc STM interface or "
                "dynamics simulator is null." );
    }
    if( columnNormalizationFactors.size( ) != totalParameterSize )
    {
        throw std::runtime_error( "Error in assembleInterArcContinuityContribution: columnNormalizationFactors size (" +
                                  std::to_string( columnNormalizationFactors.size( ) ) + ") does not match totalParameterSize (" +
                                  std::to_string( totalParameterSize ) + ")." );
    }

    const int fullStmCols = stmInterface->getFullStateTransitionMatrixSize( ) + stmInterface->getFullSensitivityMatrixSize( );
    if( fullStmCols != totalParameterSize )
    {
        throw std::runtime_error( "Error in assembleInterArcContinuityContribution: full STM column count (" +
                                  std::to_string( fullStmCols ) + ") does not match LSQ parameter size (" +
                                  std::to_string( totalParameterSize ) +
                                  "). Inter-arc continuity priors currently require a pure multi-arc parameter set "
                                  "(spec section 9: hybrid-arc out of scope for v1)." );
    }

    // First pass: compute the total constraint dimension from the rank of every active weight matrix.
    int totalConstraintDimension = 0;
    for( const auto& settings : constraintSettings )
    {
        for( const auto& body : settings->bodies( ) )
        {
            for( std::size_t pairIndex = 0; pairIndex < settings->numberOfPairsForBody( body ); ++pairIndex )
            {
                totalConstraintDimension += detail::computePsdMatrixRank( settings->weightMatrixForBodyAndPair( body, pairIndex ) );
            }
        }
    }
    if( totalConstraintDimension == 0 )
    {
        throw std::runtime_error(
                "Error in assembleInterArcContinuityContribution: total constraint dimension is zero across all "
                "settings (every weight matrix is zero). Specify a non-zero constraint weight matrix." );
    }

    const std::vector< double > arcStartTimes = multiArcSimulator->getArcStartTimes( );
    const std::vector< double > arcEndTimes = multiArcSimulator->getArcEndTimes( );
    auto multiArcResults = multiArcSimulator->getMultiArcPropagationResults( );
    const auto& perArcResults = multiArcResults->getSingleArcResults( );

    auto multiArcStateParameters = parametersToEstimate->getEstimatedMultiArcInitialStateParameters( );

    // Second pass: build each pair's continuity design matrix and accumulate its normal-equation contribution.
    for( const auto& settings : constraintSettings )
    {
        for( const auto& body : settings->bodies( ) )
        {
            // Locate the body's arc-wise initial-state parameter.
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > >
                    bodyParameter;
            for( const auto& candidate : multiArcStateParameters )
            {
                if( candidate->getParameterName( ).second.first == body )
                {
                    bodyParameter = candidate;
                    break;
                }
            }
            if( bodyParameter == nullptr )
            {
                throw std::runtime_error( "Error in assembleInterArcContinuityContribution: body \"" + body +
                                          "\" referenced by an inter-arc continuity prior does not have an arc-wise "
                                          "initial state parameter." );
            }

            const auto& epochs = settings->connectionEpochsForBody( body );
            const auto& arcPairs = settings->arcPairsForBody( body );

            if( epochs.empty( ) )
            {
                throw std::runtime_error( "Error in assembleInterArcContinuityContribution: body \"" + body +
                                          "\" has no inter-arc continuity epochs." );
            }

            const int representativeArcIndex = arcPairs.empty( ) ? 0 : arcPairs.front( ).first;
            const detail::BodyArcRows firstBodyRows = detail::getBodyRowsInMultiArcLayout( stmInterface, representativeArcIndex, body );
            if( firstBodyRows.size <= 0 )
            {
                throw std::runtime_error( "Error in assembleInterArcContinuityContribution: body \"" + body +
                                          "\" has an invalid state-block size in the multi-arc STM layout." );
            }
            if( bodyParameter->getParameterSize( ) % firstBodyRows.size != 0 )
            {
                throw std::runtime_error( "Error in assembleInterArcContinuityContribution: body \"" + body +
                                          "\" parameter size is not an integer multiple of its state-block size." );
            }
            const int numberOfArcs = static_cast< int >( bodyParameter->getParameterSize( ) ) / firstBodyRows.size;

            for( std::size_t pairIndex = 0; pairIndex < epochs.size( ); ++pairIndex )
            {
                // If the user did not provide explicit arc pairs, connect consecutive arcs in epoch order.
                const std::pair< int, int > arcPair = arcPairs.empty( )
                        ? std::make_pair( static_cast< int >( pairIndex ), static_cast< int >( pairIndex + 1 ) )
                        : arcPairs.at( pairIndex );
                if( arcPair.first < 0 || arcPair.second >= numberOfArcs )
                {
                    throw std::runtime_error( "Error in assembleInterArcContinuityContribution for body \"" + body + "\": arc pair (" +
                                              std::to_string( arcPair.first ) + ", " + std::to_string( arcPair.second ) +
                                              ") is out of range [0, " + std::to_string( numberOfArcs ) + ")." );
                }

                const double connectionEpoch = epochs[ pairIndex ];
                const detail::BodyArcRows leftBodyRows = detail::getBodyRowsInMultiArcLayout( stmInterface, arcPair.first, body );
                const detail::BodyArcRows rightBodyRows = detail::getBodyRowsInMultiArcLayout( stmInterface, arcPair.second, body );
                if( leftBodyRows.size != firstBodyRows.size || rightBodyRows.size != firstBodyRows.size )
                {
                    throw std::runtime_error( "Error in assembleInterArcContinuityContribution: body \"" + body +
                                              "\" has inconsistent state-block sizes across the multi-arc STM layout." );
                }

                // Both arcs must contain the connection epoch so that state and variational data can be evaluated
                // from the explicitly requested left and right sides of the boundary.
                for( int side = 0; side < 2; ++side )
                {
                    const int arcIndex = ( side == 0 ) ? arcPair.first : arcPair.second;
                    const double arcStart = arcStartTimes.at( arcIndex );
                    const double arcEnd = arcEndTimes.at( arcIndex );
                    if( connectionEpoch < arcStart || connectionEpoch > arcEnd )
                    {
                        throw std::runtime_error( "Inter-arc continuity connection epoch " + std::to_string( connectionEpoch ) +
                                                  " for body " + body + " and arc pair (" + std::to_string( arcPair.first ) + ", " +
                                                  std::to_string( arcPair.second ) + ") is outside the propagated interval of arc " +
                                                  std::to_string( arcIndex ) + " [" + std::to_string( arcStart ) + ", " +
                                                  std::to_string( arcEnd ) +
                                                  "]. Extend the arc propagation interval or change the "
                                                  "connection epoch." );
                    }
                }

                // Evaluate the constrained body's state from both sides of the connection epoch.
                const Eigen::VectorXd leftArcState =
                        detail::evaluateArcStateAtTime( perArcResults.at( arcPair.first )->getEquationsOfMotionNumericalSolution( ),
                                                        connectionEpoch,
                                                        arcStartTimes.at( arcPair.first ),
                                                        arcEndTimes.at( arcPair.first ),
                                                        leftBodyRows.arcWiseStateStart,
                                                        leftBodyRows.size );
                const Eigen::VectorXd rightArcState =
                        detail::evaluateArcStateAtTime( perArcResults.at( arcPair.second )->getEquationsOfMotionNumericalSolution( ),
                                                        connectionEpoch,
                                                        arcStartTimes.at( arcPair.second ),
                                                        arcEndTimes.at( arcPair.second ),
                                                        rightBodyRows.arcWiseStateStart,
                                                        rightBodyRows.size );
                const Eigen::VectorXd stateDiscrepancy = rightArcState - leftArcState;

                // Retrieve each side's full variational matrix explicitly; time-based arc lookup alone cannot
                // distinguish the left and right arcs at a shared boundary.
                Eigen::MatrixXd leftArcVariationalMatrix =
                        stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( arcPair.first, connectionEpoch );
                Eigen::MatrixXd rightArcVariationalMatrix =
                        stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( arcPair.second, connectionEpoch );
                if( leftArcVariationalMatrix.rows( ) != stmInterface->getFullStateSize( ) ||
                    leftArcVariationalMatrix.cols( ) != totalParameterSize ||
                    rightArcVariationalMatrix.rows( ) != stmInterface->getFullStateSize( ) ||
                    rightArcVariationalMatrix.cols( ) != totalParameterSize )
                {
                    throw std::runtime_error(
                            "Error in assembleInterArcContinuityContribution: per-arc STM has unexpected "
                            "shape (rows=" +
                            std::to_string( leftArcVariationalMatrix.rows( ) ) +
                            ", cols=" + std::to_string( leftArcVariationalMatrix.cols( ) ) + "), expected (" +
                            std::to_string( stmInterface->getFullStateSize( ) ) + ", " + std::to_string( totalParameterSize ) + ")." );
                }
                Eigen::MatrixXd continuityDesignMatrix =
                        rightArcVariationalMatrix.block( rightBodyRows.fullStateStart, 0, rightBodyRows.size, totalParameterSize ) -
                        leftArcVariationalMatrix.block( leftBodyRows.fullStateStart, 0, leftBodyRows.size, totalParameterSize );

                // Column-normalise the continuity design matrix using the same factors applied to the observation
                // design matrix before solving the least-squares normal equations.
                for( int col = 0; col < continuityDesignMatrix.cols( ); ++col )
                {
                    const double factor = columnNormalizationFactors( col );
                    if( factor == 0.0 )
                    {
                        throw std::runtime_error(
                                "Error in assembleInterArcContinuityContribution: column "
                                "normalization factor is zero at parameter index " +
                                std::to_string( col ) + "." );
                    }
                    continuityDesignMatrix.col( col ) /= factor;
                }

                const Eigen::MatrixXd& constraintWeightMatrix = settings->weightMatrixForBodyAndPair( body, pairIndex );
                if( constraintWeightMatrix.rows( ) != stateDiscrepancy.rows( ) )
                {
                    throw std::runtime_error( "Error in assembleInterArcContinuityContribution: weight matrix for body \"" + body +
                                              "\" has size " + std::to_string( constraintWeightMatrix.rows( ) ) +
                                              ", but the propagated state block has size " + std::to_string( stateDiscrepancy.rows( ) ) +
                                              "." );
                }
                const Eigen::MatrixXd constraintWeight =
                        ( 1.0 / ( settings->constraintScalingFactor( ) * static_cast< double >( totalConstraintDimension ) ) ) *
                        constraintWeightMatrix;

                // Symmetrise the normal-matrix contribution explicitly. Mathematically the product is symmetric, but
                // Eigen evaluates the triple product left-to-right and small numerical asymmetries can appear when
                // columns have very different magnitudes.
                const Eigen::MatrixXd pairNormalMatrixContribution =
                        continuityDesignMatrix.transpose( ) * constraintWeight * continuityDesignMatrix;
                contribution.additionalNormalMatrix.noalias( ) +=
                        0.5 * ( pairNormalMatrixContribution + pairNormalMatrixContribution.transpose( ) );
                contribution.additionalRightHandSide.noalias( ) -=
                        continuityDesignMatrix.transpose( ) * ( constraintWeight * stateDiscrepancy );
                contribution.totalConstraintCost += ( stateDiscrepancy.transpose( ) * constraintWeight * stateDiscrepancy )( 0, 0 );
                contribution.perPairDiscrepancies.push_back( stateDiscrepancy );
            }
        }
    }

    return contribution;
}

//! Assemble soft inter-arc continuity-prior terms from generic OD manager interfaces.
/*!
 *  This helper centralizes the pure multi-arc casts and diagnostics used by both parameter estimation and covariance
 *  analysis. It returns an empty no-op contribution when no continuity priors are configured.
 */
template< typename ObservationScalarType, typename TimeType >
InterArcConstraintContribution assembleInterArcContinuityContributionFromManagerInterfaces(
        const std::vector< std::shared_ptr< InterArcStateContinuityConstraintSettings > >& constraintSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >& parametersToEstimate,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface >& stateTransitionInterface,
        const std::shared_ptr< propagators::VariationalEquationsSolver< ObservationScalarType, TimeType > >& variationalEquationsSolver,
        const Eigen::VectorXd& columnNormalizationFactors,
        const int totalParameterSize,
        const std::string& context )
{
    InterArcConstraintContribution contribution;
    contribution.additionalNormalMatrix = Eigen::MatrixXd( 0, 0 );
    contribution.additionalRightHandSide = Eigen::VectorXd( 0 );
    contribution.totalConstraintCost = 0.0;
    contribution.perPairDiscrepancies.clear( );

    if( constraintSettings.empty( ) )
    {
        return contribution;
    }

    auto multiArcStmInterface =
            std::dynamic_pointer_cast< propagators::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< ObservationScalarType > >(
                    stateTransitionInterface );
    if( multiArcStmInterface == nullptr )
    {
        throw std::runtime_error( "Error when applying inter-arc continuity priors in " + context +
                                  ": state-transition matrix interface is not a "
                                  "MultiArcCombinedStateTransitionAndSensitivityMatrixInterface. Inter-arc continuity priors "
                                  "are only supported for pure multi-arc estimators." );
    }
    if( variationalEquationsSolver == nullptr )
    {
        throw std::runtime_error( "Error when applying inter-arc continuity priors in " + context +
                                  ": variational equations solver is null." );
    }

    auto multiArcSimulator = std::dynamic_pointer_cast< propagators::MultiArcDynamicsSimulator< ObservationScalarType, TimeType > >(
            variationalEquationsSolver->getDynamicsSimulatorBase( ) );
    if( multiArcSimulator == nullptr )
    {
        throw std::runtime_error( "Error when applying inter-arc continuity priors in " + context +
                                  ": dynamics simulator is not a MultiArcDynamicsSimulator." );
    }

    return assembleInterArcContinuityContribution< ObservationScalarType, TimeType >( constraintSettings,
                                                                                      parametersToEstimate,
                                                                                      multiArcSimulator,
                                                                                      multiArcStmInterface,
                                                                                      columnNormalizationFactors,
                                                                                      totalParameterSize );
}

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_INTERARCCONTINUITYCONSTRAINT_H
