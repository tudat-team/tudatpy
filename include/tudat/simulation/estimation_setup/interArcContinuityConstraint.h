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
#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"

namespace tudat
{
namespace simulation_setup
{

//! Per-iteration aggregated contribution of all inter-arc continuity constraints to the normal equations.
struct InterArcConstraintContribution {
    Eigen::MatrixXd additionalNormalMatrix;
    Eigen::VectorXd additionalRightHandSide;
    double totalConstraintCost = 0.0;
    std::vector< Eigen::Matrix< double, 6, 1 > > perPairDiscrepancies;
};

namespace detail
{

//! Numerical rank of a 6x6 PSD matrix using its eigenvalue spectrum.
inline int rankOf6x6PsdMatrix( const Eigen::Matrix< double, 6, 6 >& C )
{
    Eigen::SelfAdjointEigenSolver< Eigen::Matrix< double, 6, 6 > > solver( 0.5 * ( C + C.transpose( ) ) );
    if( solver.info( ) != Eigen::Success )
    {
        throw std::runtime_error( "Error in rankOf6x6PsdMatrix: eigen-decomposition failed." );
    }
    const double largest = std::max( solver.eigenvalues( ).maxCoeff( ), 0.0 );
    if( largest <= 0.0 )
    {
        return 0;
    }
    const double threshold = 1.0E-12 * largest;
    int rank = 0;
    for( int i = 0; i < 6; ++i )
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
Eigen::Matrix< double, 6, 1 > evaluateArcStateAtTime(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& arcSolution,
        const double evaluationTime,
        const double arcInitialTime,
        const double arcFinalTime,
        const int stateOffset )
{
    if( arcSolution.empty( ) )
    {
        throw std::runtime_error( "Error in evaluateArcStateAtTime: arc solution is empty." );
    }

    auto extract = [ stateOffset ]( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& fullState ) {
        Eigen::Matrix< double, 6, 1 > result;
        for( int i = 0; i < 6; ++i )
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

//! Assemble the normal-equation contribution of all inter-arc continuity constraints for the current iteration.
//!
//! Per Lari et al. (2021) Eq. 28: for each constrained boundary (k_left, k_right) at epoch t_c,
//!   d = x_right(t_c) - x_left(t_c)                              (6-vector, physical units)
//!   D = M_right(t_c) - M_left(t_c)                              (6 x N_params, full multi-arc layout)
//!   W_d = (1 / (mu_pair * m_d_total)) * C_pair
//!   H_constraint += D_norm^T W_d D_norm,  g_constraint += - D_norm^T W_d d
//! where m_d_total is a single global rank sum across every settings entry passed in. Column-normalisation of D
//! uses the same factors that the OD loop applies to the observation design matrix (spec section 2.4).
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
                                  "). Inter-arc continuity constraints currently require a pure multi-arc parameter set "
                                  "(spec section 9: hybrid-arc out of scope for v1)." );
    }

    // First pass: compute the global m_d_total = sum over every (settings, pair) of rank(C_pair).
    int mDTotal = 0;
    for( const auto& settings : constraintSettings )
    {
        for( std::size_t pairIndex = 0; pairIndex < settings->numberOfPairs( ); ++pairIndex )
        {
            mDTotal += detail::rankOf6x6PsdMatrix( settings->weightMatrixForPair( pairIndex ) );
        }
    }
    if( mDTotal == 0 )
    {
        throw std::runtime_error(
                "Error in assembleInterArcContinuityContribution: total m_d is zero across all "
                "settings (every weight matrix is zero). Specify a non-zero C." );
    }

    const std::vector< double > arcStartTimes = multiArcSimulator->getArcStartTimes( );
    const std::vector< double > arcEndTimes = multiArcSimulator->getArcEndTimes( );
    auto multiArcResults = multiArcSimulator->getMultiArcPropagationResults( );
    const auto& perArcResults = multiArcResults->getSingleArcResults( );

    auto multiArcStateParameters = parametersToEstimate->getEstimatedMultiArcInitialStateParameters( );

    // Second pass: build D and accumulate contributions per (settings, pair).
    for( const auto& settings : constraintSettings )
    {
        // Locate the body's arc-wise translational state parameter.
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > >
                bodyParameter;
        for( const auto& candidate : multiArcStateParameters )
        {
            if( candidate->getParameterName( ).second.first == settings->body( ) )
            {
                bodyParameter = candidate;
                break;
            }
        }
        if( bodyParameter == nullptr )
        {
            throw std::runtime_error( "Error in assembleInterArcContinuityContribution: body \"" + settings->body( ) +
                                      "\" referenced by an inter-arc continuity constraint does not have an arc-wise "
                                      "translational initial state parameter." );
        }

        const int nArcs = static_cast< int >( bodyParameter->getParameterSize( ) ) / 6;

        const auto arcPairs = settings->arcPairs( );
        const auto& epochs = settings->connectionEpochs( );

        for( std::size_t pairIndex = 0; pairIndex < epochs.size( ); ++pairIndex )
        {
            const std::pair< int, int > pair = arcPairs.empty( )
                    ? std::make_pair( static_cast< int >( pairIndex ), static_cast< int >( pairIndex + 1 ) )
                    : arcPairs.at( pairIndex );
            if( pair.first < 0 || pair.second >= nArcs )
            {
                throw std::runtime_error( "Error in assembleInterArcContinuityContribution for body \"" + settings->body( ) +
                                          "\": arc pair (" + std::to_string( pair.first ) + ", " + std::to_string( pair.second ) +
                                          ") is out of range [0, " + std::to_string( nArcs ) + ")." );
            }

            const double tC = epochs[ pairIndex ];
            const detail::BodyArcRows leftBodyRows = detail::getBodyRowsInMultiArcLayout( stmInterface, pair.first, settings->body( ) );
            const detail::BodyArcRows rightBodyRows = detail::getBodyRowsInMultiArcLayout( stmInterface, pair.second, settings->body( ) );
            if( leftBodyRows.size != 6 || rightBodyRows.size != 6 )
            {
                throw std::runtime_error( "Error in assembleInterArcContinuityContribution: body \"" + settings->body( ) +
                                          "\" does not have a 6-state translational block in the multi-arc STM layout." );
            }
            for( int side = 0; side < 2; ++side )
            {
                const int arcIndex = ( side == 0 ) ? pair.first : pair.second;
                const double arcStart = arcStartTimes.at( arcIndex );
                const double arcEnd = arcEndTimes.at( arcIndex );
                if( tC < arcStart || tC > arcEnd )
                {
                    throw std::runtime_error( "Inter-arc continuity connection epoch " + std::to_string( tC ) + " for body " +
                                              settings->body( ) + " and arc pair (" + std::to_string( pair.first ) + ", " +
                                              std::to_string( pair.second ) + ") is outside the propagated interval of arc " +
                                              std::to_string( arcIndex ) + " [" + std::to_string( arcStart ) + ", " +
                                              std::to_string( arcEnd ) +
                                              "]. Extend the arc propagation interval or change the "
                                              "connection epoch." );
                }
            }

            const Eigen::Matrix< double, 6, 1 > xLeft =
                    detail::evaluateArcStateAtTime( perArcResults.at( pair.first )->getEquationsOfMotionNumericalSolution( ),
                                                    tC,
                                                    arcStartTimes.at( pair.first ),
                                                    arcEndTimes.at( pair.first ),
                                                    leftBodyRows.arcWiseStateStart );
            const Eigen::Matrix< double, 6, 1 > xRight =
                    detail::evaluateArcStateAtTime( perArcResults.at( pair.second )->getEquationsOfMotionNumericalSolution( ),
                                                    tC,
                                                    arcStartTimes.at( pair.second ),
                                                    arcEndTimes.at( pair.second ),
                                                    rightBodyRows.arcWiseStateStart );
            const Eigen::Matrix< double, 6, 1 > d = xRight - xLeft;

            Eigen::MatrixXd mLeft = stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( pair.first, tC );
            Eigen::MatrixXd mRight = stmInterface->getFullCombinedStateTransitionAndSensitivityMatrixForArc( pair.second, tC );
            if( mLeft.rows( ) != stmInterface->getFullStateSize( ) || mLeft.cols( ) != totalParameterSize ||
                mRight.rows( ) != stmInterface->getFullStateSize( ) || mRight.cols( ) != totalParameterSize )
            {
                throw std::runtime_error(
                        "Error in assembleInterArcContinuityContribution: per-arc STM has unexpected "
                        "shape (rows=" +
                        std::to_string( mLeft.rows( ) ) + ", cols=" + std::to_string( mLeft.cols( ) ) + "), expected (" +
                        std::to_string( stmInterface->getFullStateSize( ) ) + ", " + std::to_string( totalParameterSize ) + ")." );
            }
            Eigen::MatrixXd D = mRight.block( rightBodyRows.fullStateStart, 0, 6, totalParameterSize ) -
                    mLeft.block( leftBodyRows.fullStateStart, 0, 6, totalParameterSize );  // 6 x N_params

            // Column-normalise D using the same factors applied to the observation design matrix.
            for( int col = 0; col < D.cols( ); ++col )
            {
                const double factor = columnNormalizationFactors( col );
                if( factor == 0.0 )
                {
                    throw std::runtime_error(
                            "Error in assembleInterArcContinuityContribution: column "
                            "normalisation factor is zero at parameter index " +
                            std::to_string( col ) + "." );
                }
                D.col( col ) /= factor;
            }

            const Eigen::Matrix< double, 6, 6 >& C = settings->weightMatrixForPair( pairIndex );
            const double mu = settings->muForPair( pairIndex );
            const Eigen::Matrix< double, 6, 6 > weight = ( 1.0 / ( mu * static_cast< double >( mDTotal ) ) ) * C;

            // Symmetrise the rank-one normal-matrix contribution explicitly. Mathematically D^T W D is symmetric
            // (W is symmetric), but Eigen evaluates the triple product left-to-right and the resulting numerical
            // asymmetry can be non-trivial when D has columns of disparate magnitudes (e.g. position and velocity
            // STM entries mixed together).
            const Eigen::MatrixXd Hpair = D.transpose( ) * weight * D;
            contribution.additionalNormalMatrix.noalias( ) += 0.5 * ( Hpair + Hpair.transpose( ) );
            contribution.additionalRightHandSide.noalias( ) -= D.transpose( ) * ( weight * d );
            contribution.totalConstraintCost += d.transpose( ) * weight * d;
            contribution.perPairDiscrepancies.push_back( d );
        }
    }

    return contribution;
}

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_INTERARCCONTINUITYCONSTRAINT_H
