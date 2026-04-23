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
//! Validate/initialize the inverse a priori covariance matrix input.
/*!
 * If the input matrix is 0x0, a zero matrix with the required dimensions is created.
 * Otherwise, the function verifies that the input matrix is square and matches the
 * full estimated-parameter dimension.
 * \param inverseAprioriCovariance Input inverse a priori covariance matrix.
 * \param numberOfEstimatedParameters Required matrix size based on the parameter set.
 * \return Validated inverse a priori covariance matrix.
 */
Eigen::MatrixXd getValidatedInverseAprioriCovarianceInput(
        const Eigen::MatrixXd& inverseAprioriCovariance,
        const int numberOfEstimatedParameters );

//! Expand and validate uncertainties for a single estimatable-parameter block.
/*!
 * Accepts either a scalar uncertainty (size-1 vector) or one uncertainty per block component.
 * Scalar input is broadcast to the block size. The resulting uncertainties are checked to be
 * non-empty and strictly positive.
 * \param uncertaintyValues Input uncertainty values for a parameter.
 * \param parameterSize Size of the matched parameter block.
 * \param parameterIdentifier Identifier of the matched parameter (for error reporting).
 * \return Uncertainty vector of size \p parameterSize.
 */
Eigen::VectorXd getUncertaintyValuesForParameter(
        const Eigen::VectorXd& uncertaintyValues,
        const int parameterSize,
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier );

//! Add/update inverse a priori covariance entries from per-parameter uncertainties.
/*!
 * For each entry in \p aprioriUncertaintyPerParameter, this function resolves the matching
 * parameter block(s) in \p parameterSet and sets the corresponding diagonal terms in the
 * inverse a priori covariance matrix to :math:`1/\sigma^2`.
 *
 * If \p inverseAprioriCovariance is 0x0, a zero matrix with the full parameter-set size
 * is created first.
 *
 * \param inverseAprioriCovariance Existing inverse a priori covariance matrix (or 0x0 matrix).
 * \param parameterSet Estimatable-parameter set used for block index lookup.
 * \param aprioriUncertaintyPerParameter Per-parameter uncertainty entries.
 * \return Updated inverse a priori covariance matrix.
 */
template< typename InitialStateParameterType = double >
Eigen::MatrixXd addInverseAprioriCovarianceEntries(
        const Eigen::MatrixXd& inverseAprioriCovariance,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                aprioriUncertaintyPerParameter )
{
    if( parameterSet == nullptr )
    {
        throw std::runtime_error( "Error when creating/updating inverse apriori covariance: parameter_set is null." );
    }

    Eigen::MatrixXd updatedInverseAprioriCovariance =
            getValidatedInverseAprioriCovarianceInput( inverseAprioriCovariance, parameterSet->getEstimatedParameterSetSize( ) );

    for( const std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd >& aprioriEntry :
         aprioriUncertaintyPerParameter )
    {
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier = aprioriEntry.first;
        const Eigen::VectorXd& uncertaintyValues = aprioriEntry.second;

        const std::vector< std::pair< int, int > > parameterIndices =
                parameterSet->getParameterIndicesForParameterIdentifier( parameterIdentifier );

        for( const std::pair< int, int >& indexAndSize : parameterIndices )
        {
            const int startIndex = indexAndSize.first;
            const int parameterSize = indexAndSize.second;
            const Eigen::VectorXd currentUncertaintyValues =
                    getUncertaintyValuesForParameter( uncertaintyValues, parameterSize, parameterIdentifier );

            for( int i = 0; i < parameterSize; i++ )
            {
                updatedInverseAprioriCovariance( startIndex + i, startIndex + i ) =
                        1.0 / ( currentUncertaintyValues( i ) * currentUncertaintyValues( i ) );
            }
        }
    }

    return updatedInverseAprioriCovariance;
}

//! Create a full inverse a priori covariance matrix from per-parameter uncertainties.
/*!
 * Convenience wrapper around addInverseAprioriCovarianceEntries that starts from an empty
 * matrix and creates a correctly sized matrix automatically.
 * \param parameterSet Estimatable-parameter set used for block index lookup.
 * \param aprioriUncertaintyPerParameter Per-parameter uncertainty entries.
 * \return Inverse a priori covariance matrix with constrained diagonal entries set.
 */
template< typename InitialStateParameterType = double >
Eigen::MatrixXd createInverseAprioriCovariance(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                aprioriUncertaintyPerParameter )
{
    return addInverseAprioriCovarianceEntries< InitialStateParameterType >(
            Eigen::MatrixXd::Zero( 0, 0 ), parameterSet, aprioriUncertaintyPerParameter );
}

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_CREATE_INVERSE_APRIORI_COVARIANCE_H
