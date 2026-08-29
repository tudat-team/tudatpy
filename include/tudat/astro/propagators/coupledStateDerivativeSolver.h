/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#ifndef TUDAT_COUPLED_STATE_DERIVATIVE_SOLVER_H
#define TUDAT_COUPLED_STATE_DERIVATIVE_SOLVER_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>

#include <Eigen/Core>
#include <Eigen/LU>

namespace tudat
{

namespace propagators
{

//! Policy used when the coupled state-derivative equations cannot be solved to tolerance.
enum CoupledStateDerivativeFailureHandling {
    throw_exception_on_coupled_derivative_failure = 0,
    accept_last_coupled_derivative_iteration = 1
};

//! Settings for solving algebraic dependencies between propagated state derivatives.
class CoupledStateDerivativeSolverSettings
{
public:
    CoupledStateDerivativeSolverSettings(
            const bool useDirectAffineSolution = true,
            const double relativeTolerance = 1.0e-11,
            const double absoluteTolerance = 1.0e-13,
            const unsigned int maximumIterations = 25,
            const CoupledStateDerivativeFailureHandling failureHandling = throw_exception_on_coupled_derivative_failure ):
        useDirectAffineSolution_( useDirectAffineSolution ), relativeTolerance_( relativeTolerance ),
        absoluteTolerance_( absoluteTolerance ), maximumIterations_( maximumIterations ), failureHandling_( failureHandling )
    {
        validate( );
    }

    void validate( ) const
    {
        if( relativeTolerance_ < 0.0 || absoluteTolerance_ < 0.0 )
        {
            throw std::invalid_argument( "Coupled state-derivative solver tolerances must be non-negative." );
        }
        if( maximumIterations_ == 0 )
        {
            throw std::invalid_argument( "Coupled state-derivative solver maximum iteration count must be positive." );
        }
    }

    bool useDirectAffineSolution_;
    double relativeTolerance_;
    double absoluteTolerance_;
    unsigned int maximumIterations_;
    CoupledStateDerivativeFailureHandling failureHandling_;
};

//! Result of a coupled state-derivative solve.
template< typename ScalarType >
struct CoupledStateDerivativeSolution {
    Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > stateDerivative_;
    Eigen::Matrix< ScalarType, Eigen::Dynamic, Eigen::Dynamic > implicitDerivativeMultiplier_;
    bool usedDirectSolution_ = false;
    bool converged_ = false;
    bool implicitDerivativeMultiplierIsValid_ = true;
    unsigned int numberOfIterations_ = 0;
};

namespace detail
{

template< typename ScalarType >
bool isCoupledDerivativeSolutionConverged( const Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >& candidate,
                                           const Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >& mappedCandidate,
                                           const Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >& componentScales,
                                           const CoupledStateDerivativeSolverSettings& settings )
{
    for( int i = 0; i < candidate.rows( ); ++i )
    {
        const ScalarType relativeScale = std::max( std::abs( candidate( i ) ), std::abs( mappedCandidate( i ) ) );
        const ScalarType tolerance = static_cast< ScalarType >( settings.absoluteTolerance_ ) * componentScales( i ) +
                static_cast< ScalarType >( settings.relativeTolerance_ ) * relativeScale;
        if( std::abs( candidate( i ) - mappedCandidate( i ) ) > tolerance )
        {
            return false;
        }
    }
    return true;
}

}  // namespace detail

//! Solve a coupled derivative fixed point, using a direct affine solve when requested.
/*!
 * The supplied function evaluates the mapping F(d) from a trial coupled derivative d to
 * the derivative obtained after updating all dependent environment and derivative models.
 * Component scales nondimensionalise the mixed gravity-coefficient-rate and angular-
 * acceleration entries. If the direct affine reconstruction is not valid to the requested
 * tolerance, a scaled fixed-point iteration is used as a fallback.
 */
template< typename ScalarType >
CoupledStateDerivativeSolution< ScalarType > solveCoupledStateDerivative(
        const std::function< Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >( const Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >& ) >&
                derivativeMapping,
        const Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >& initialDerivative,
        const Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >& componentScales,
        const CoupledStateDerivativeSolverSettings& settings )
{
    using VectorType = Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >;
    using MatrixType = Eigen::Matrix< ScalarType, Eigen::Dynamic, Eigen::Dynamic >;

    settings.validate( );
    if( initialDerivative.rows( ) != componentScales.rows( ) )
    {
        throw std::invalid_argument( "Coupled state-derivative values and scales have incompatible sizes." );
    }
    if( ( componentScales.array( ) <= static_cast< ScalarType >( 0 ) ).any( ) || !componentScales.allFinite( ) )
    {
        throw std::invalid_argument( "Coupled state-derivative scales must be finite and strictly positive." );
    }

    CoupledStateDerivativeSolution< ScalarType > solution;
    const int numberOfComponents = initialDerivative.rows( );
    solution.implicitDerivativeMultiplier_ = MatrixType::Identity( numberOfComponents, numberOfComponents );
    if( numberOfComponents == 0 )
    {
        solution.stateDerivative_ = initialDerivative;
        solution.converged_ = true;
        return solution;
    }

    const auto calculateScaledCouplingMatrix = [ & ]( const VectorType& centre ) {
        MatrixType scaledCouplingMatrix( numberOfComponents, numberOfComponents );
        for( int column = 0; column < numberOfComponents; ++column )
        {
            VectorType positivePerturbation = centre;
            VectorType negativePerturbation = centre;
            positivePerturbation( column ) += componentScales( column );
            negativePerturbation( column ) -= componentScales( column );
            scaledCouplingMatrix.col( column ) = ( derivativeMapping( positivePerturbation ) - derivativeMapping( negativePerturbation ) )
                                                         .cwiseQuotient( componentScales ) /
                    static_cast< ScalarType >( 2 );
        }
        return scaledCouplingMatrix;
    };

    if( settings.useDirectAffineSolution_ )
    {
        const VectorType zeroDerivative = VectorType::Zero( numberOfComponents );
        const VectorType constantTerm = derivativeMapping( zeroDerivative );
        const MatrixType scaledCouplingMatrix = calculateScaledCouplingMatrix( zeroDerivative );

        const MatrixType scaledSystemMatrix = MatrixType::Identity( numberOfComponents, numberOfComponents ) - scaledCouplingMatrix;
        Eigen::FullPivLU< MatrixType > decomposition( scaledSystemMatrix );
        if( decomposition.isInvertible( ) )
        {
            const VectorType scaledConstantTerm = constantTerm.cwiseQuotient( componentScales );
            const VectorType scaledDerivative = decomposition.solve( scaledConstantTerm );
            const VectorType directDerivative = scaledDerivative.cwiseProduct( componentScales );
            const VectorType mappedDirectDerivative = derivativeMapping( directDerivative );

            if( directDerivative.allFinite( ) && mappedDirectDerivative.allFinite( ) &&
                detail::isCoupledDerivativeSolutionConverged( directDerivative, mappedDirectDerivative, componentScales, settings ) )
            {
                solution.stateDerivative_ = directDerivative;
                solution.implicitDerivativeMultiplier_ =
                        componentScales.asDiagonal( ) * decomposition.inverse( ) * componentScales.cwiseInverse( ).asDiagonal( );
                solution.usedDirectSolution_ = true;
                solution.converged_ = true;
                return solution;
            }
        }
    }

    VectorType currentDerivative = initialDerivative;
    for( unsigned int iteration = 0; iteration < settings.maximumIterations_; ++iteration )
    {
        const VectorType mappedDerivative = derivativeMapping( currentDerivative );
        solution.numberOfIterations_ = iteration + 1;
        if( detail::isCoupledDerivativeSolutionConverged( currentDerivative, mappedDerivative, componentScales, settings ) )
        {
            solution.stateDerivative_ = mappedDerivative;
            solution.converged_ = true;
            const MatrixType scaledCouplingMatrix = calculateScaledCouplingMatrix( mappedDerivative );
            const MatrixType scaledSystemMatrix = MatrixType::Identity( numberOfComponents, numberOfComponents ) - scaledCouplingMatrix;
            Eigen::FullPivLU< MatrixType > decomposition( scaledSystemMatrix );
            if( decomposition.isInvertible( ) )
            {
                solution.implicitDerivativeMultiplier_ =
                        componentScales.asDiagonal( ) * decomposition.inverse( ) * componentScales.cwiseInverse( ).asDiagonal( );
            }
            else
            {
                solution.implicitDerivativeMultiplierIsValid_ = false;
            }
            return solution;
        }
        currentDerivative = mappedDerivative;
    }

    solution.stateDerivative_ = currentDerivative;
    solution.implicitDerivativeMultiplierIsValid_ = false;
    if( settings.failureHandling_ == throw_exception_on_coupled_derivative_failure )
    {
        throw std::runtime_error( "Coupled state-derivative solution did not converge in " + std::to_string( settings.maximumIterations_ ) +
                                  " iterations." );
    }
    return solution;
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_COUPLED_STATE_DERIVATIVE_SOLVER_H
