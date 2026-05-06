/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/createInverseAprioriCovariance.h"

#include <stdexcept>
#include <string>

namespace tudat
{
namespace simulation_setup
{

//! Validate matrix dimensions or create a correctly sized zero matrix from a 0x0 input.
Eigen::MatrixXd getValidatedCovarianceInput(
        const Eigen::MatrixXd& covarianceMatrix,
        const int numberOfEstimatedParameters )
{
    if( covarianceMatrix.rows( ) == 0 && covarianceMatrix.cols( ) == 0 )
    {
        return Eigen::MatrixXd::Zero( numberOfEstimatedParameters, numberOfEstimatedParameters );
    }

    if( covarianceMatrix.rows( ) != numberOfEstimatedParameters ||
        covarianceMatrix.cols( ) != numberOfEstimatedParameters )
    {
        throw std::runtime_error( "Error when creating/updating covariance matrix entries: provided matrix has size " +
                                  std::to_string( covarianceMatrix.rows( ) ) + "x" +
                                  std::to_string( covarianceMatrix.cols( ) ) + ", expected " +
                                  std::to_string( numberOfEstimatedParameters ) + "x" +
                                  std::to_string( numberOfEstimatedParameters ) + "." );
    }

    return covarianceMatrix;
}

//! Convert scalar/vector covariance diagonal input to per-component values and validate non-negativity.
Eigen::VectorXd getCovarianceDiagonalValuesForParameter(
        const Eigen::VectorXd& covarianceDiagonalValues,
        const int parameterSize,
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier )
{
    Eigen::VectorXd currentCovarianceDiagonalValues;

    if( covarianceDiagonalValues.size( ) == 0 )
    {
        throw std::runtime_error(
                "Error when applying covariance diagonal entries for parameter type '" +
                estimatable_parameters::getParameterTypeString( parameterIdentifier.first ) + "' with identifiers ('" +
                parameterIdentifier.second.first + "', '" + parameterIdentifier.second.second + "'): input is empty." );
    }

    if( covarianceDiagonalValues.size( ) == 1 )
    {
        currentCovarianceDiagonalValues = Eigen::VectorXd::Constant( parameterSize, covarianceDiagonalValues( 0 ) );
    }
    else if( covarianceDiagonalValues.size( ) == parameterSize )
    {
        currentCovarianceDiagonalValues = covarianceDiagonalValues;
    }
    else
    {
        throw std::runtime_error(
                "Error when applying covariance diagonal entries for parameter type '" +
                estimatable_parameters::getParameterTypeString( parameterIdentifier.first ) + "' with identifiers ('" +
                parameterIdentifier.second.first + "', '" + parameterIdentifier.second.second +
                "'): value vector size is incompatible with parameter size. "
                "Received vector of size " +
                std::to_string( covarianceDiagonalValues.size( ) ) + " for parameter size " +
                std::to_string( parameterSize ) + ". Use a scalar value or a vector matching the parameter size." );
    }

    for( int i = 0; i < currentCovarianceDiagonalValues.size( ); i++ )
    {
        if( !( currentCovarianceDiagonalValues( i ) >= 0.0 ) )
        {
            throw std::runtime_error(
                    "Error when applying covariance diagonal entries for parameter type '" +
                    estimatable_parameters::getParameterTypeString( parameterIdentifier.first ) + "' with identifiers ('" +
                    parameterIdentifier.second.first + "', '" + parameterIdentifier.second.second +
                    "'): entries must be non-negative." );
        }
    }

    return currentCovarianceDiagonalValues;
}

}  // namespace simulation_setup
}  // namespace tudat
