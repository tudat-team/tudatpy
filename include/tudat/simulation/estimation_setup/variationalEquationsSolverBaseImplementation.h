/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VARIATIONALEQUATIONSSOLVER_BASE_IMPLEMENTATION_H
#define TUDAT_VARIATIONALEQUATIONSSOLVER_BASE_IMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/variationalEquationsSolverBase.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< StateScalarType, TimeType >::value, int >::type Dummy >
typename VariationalEquationsSolver< StateScalarType, TimeType, Dummy >::MatrixType
VariationalEquationsSolver< StateScalarType, TimeType, Dummy >::createInitialConditions( const VectorType initialStateEstimate )
{
    if( stateTransitionMatrixSize_ != initialStateEstimate.rows( ) )
    {
        throw std::runtime_error( "Error when getting initial condition for variational equations, sizes are incompatible: " +
                                  std::to_string( stateTransitionMatrixSize_ ) + ", " + std::to_string( initialStateEstimate.rows( ) ) );
    }

    // Initialize initial conditions to zeros.
    MatrixType varSystemInitialState = MatrixType( stateTransitionMatrixSize_, propagatedParameterVectorSize_ + 1 ).setZero( );

    // Set initial state transition matrix to identity
    varSystemInitialState.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).setIdentity( );

    // Set initial body states to current estimate of initial body states.
    varSystemInitialState.block( 0, propagatedParameterVectorSize_, stateTransitionMatrixSize_, 1 ) = initialStateEstimate;

    return varSystemInitialState;
}

template< typename StateScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< StateScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd VariationalEquationsSolver< StateScalarType, TimeType, Dummy >::createInitialVariationalEquationsSolution( )
{
    // Initialize initial conditions to zeros.
    Eigen::MatrixXd varSystemInitialState = Eigen::MatrixXd::Zero( stateTransitionMatrixSize_, propagatedParameterVectorSize_ );

    // Set initial state transition matrix to identity
    varSystemInitialState.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).setIdentity( );

    return varSystemInitialState;
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_VARIATIONALEQUATIONSSOLVER_BASE_IMPLEMENTATION_H
