/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATE_INVERSE_APRIORI_COVARIANCE_H
#define TUDAT_CREATE_INVERSE_APRIORI_COVARIANCE_H

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"

namespace tudat
{
namespace simulation_setup
{
//! Validate/initialize a covariance-like matrix input.
/*!
 * If the input matrix is 0x0, a zero matrix with the required dimensions is created.
 * Otherwise, the function verifies that the input matrix is square and matches the
 * full estimated-parameter dimension.
 * \param covarianceMatrix Input covariance-like matrix.
 * \param numberOfEstimatedParameters Required matrix size based on the parameter set.
 * \return Validated covariance-like matrix.
 */
Eigen::MatrixXd getValidatedCovarianceInput(
        const Eigen::MatrixXd& covarianceMatrix,
        const int numberOfEstimatedParameters );

//! Expand and validate covariance diagonal values for a single estimatable-parameter block.
/*!
 * Accepts either a scalar value (size-1 vector) or one value per block component.
 * Scalar input is broadcast to the block size. The resulting vector is checked to be
 * non-empty and non-negative.
 * \param covarianceDiagonalValues Input covariance diagonal values for a parameter.
 * \param parameterSize Size of the matched parameter block.
 * \param parameterIdentifier Identifier of the matched parameter (for error reporting).
 * \return Covariance diagonal values of size \p parameterSize.
 */
Eigen::VectorXd getCovarianceDiagonalValuesForParameter(
        const Eigen::VectorXd& covarianceDiagonalValues,
        const int parameterSize,
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier );

//! Add/update diagonal entries in a covariance-like matrix from per-parameter vectors.
/*!
 * For each entry in \p covarianceDiagonalEntriesPerParameter, this function resolves the matching
 * parameter block(s) in \p parameterSet and sets the corresponding diagonal terms in the
 * covariance-like matrix directly to the provided values.
 *
 * If \p covarianceMatrix is 0x0, a zero matrix with the full parameter-set size
 * is created first.
 *
 * This function is representation-agnostic: it can be used to build a covariance matrix,
 * an inverse covariance matrix, or any diagonal matrix with parameter-aligned entries.
 *
 * \param covarianceMatrix Existing covariance-like matrix (or 0x0 matrix).
 * \param parameterSet Estimatable-parameter set used for block index lookup.
 * \param covarianceDiagonalEntriesPerParameter Per-parameter diagonal-entry vectors.
 * \return Updated covariance-like matrix.
 */
template< typename InitialStateParameterType = double >
Eigen::MatrixXd addCovarianceDiagonalEntries(
        const Eigen::MatrixXd& covarianceMatrix,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                covarianceDiagonalEntriesPerParameter )
{
    if( parameterSet == nullptr )
    {
        throw std::runtime_error( "Error when creating/updating covariance matrix entries: parameter_set is null." );
    }

    Eigen::MatrixXd updatedCovarianceMatrix =
            getValidatedCovarianceInput( covarianceMatrix, parameterSet->getEstimatedParameterSetSize( ) );

    for( const std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd >& aprioriEntry :
         covarianceDiagonalEntriesPerParameter )
    {
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier = aprioriEntry.first;
        const Eigen::VectorXd& covarianceDiagonalValues = aprioriEntry.second;

        const std::vector< std::pair< int, int > > parameterIndices =
                parameterSet->getParameterIndicesForParameterIdentifier( parameterIdentifier );

        for( const std::pair< int, int >& indexAndSize : parameterIndices )
        {
            const int startIndex = indexAndSize.first;
            const int parameterSize = indexAndSize.second;
            const Eigen::VectorXd currentCovarianceDiagonalValues =
                    getCovarianceDiagonalValuesForParameter( covarianceDiagonalValues, parameterSize, parameterIdentifier );

            for( int i = 0; i < parameterSize; i++ )
            {
                updatedCovarianceMatrix( startIndex + i, startIndex + i ) = currentCovarianceDiagonalValues( i );
            }
        }
    }

    return updatedCovarianceMatrix;
}

//! Create a covariance-like matrix from per-parameter diagonal-entry vectors.
/*!
 * Convenience wrapper around addCovarianceDiagonalEntries that starts from an empty
 * matrix and creates a correctly sized matrix automatically.
 *
 * This function is representation-agnostic: it can create a covariance matrix,
 * an inverse covariance matrix, or any parameter-aligned diagonal matrix.
 *
 * \param parameterSet Estimatable-parameter set used for block index lookup.
 * \param covarianceDiagonalEntriesPerParameter Per-parameter diagonal-entry vectors.
 * \return Covariance-like matrix with constrained diagonal entries set.
 */
template< typename InitialStateParameterType = double >
Eigen::MatrixXd createCovarianceFromDiagonalEntries(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                covarianceDiagonalEntriesPerParameter )
{
    return addCovarianceDiagonalEntries< InitialStateParameterType >(
            Eigen::MatrixXd::Zero( 0, 0 ), parameterSet, covarianceDiagonalEntriesPerParameter );
}

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_CREATE_INVERSE_APRIORI_COVARIANCE_H
