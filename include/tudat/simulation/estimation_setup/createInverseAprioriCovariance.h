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
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/io/trackingSupplementaryData.h"

namespace tudat
{
namespace simulation_setup
{

//! Extract the per-picture camera pointing inverse variances from tracking supplementary data of any source format.
std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
getCameraPointingInverseAprioriCovarianceEntries(
        const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData );

//! Function to get the inverse variances of reference point positions (e.g. landmarks), for use in an a-priori covariance.
/*!
 *  Gives one entry per reference point with a position uncertainty, for its ground_station_position parameter.
 *  \param supplementaryData Tracking supplementary data holding the position uncertainties.
 *  \return Inverse variance of each body-fixed position component, per reference point.
 */
std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
getReferencePointPositionInverseAprioriCovarianceEntries(
        const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData );

//! Function to get the variances of reference point positions (e.g. landmarks).
/*!
 *  Same as getReferencePointPositionInverseAprioriCovarianceEntries, but gives the variances themselves.
 *  \param supplementaryData Tracking supplementary data holding the position uncertainties.
 *  \return Variance of each body-fixed position component, per reference point.
 */
std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
getReferencePointPositionAprioriCovarianceEntries(
        const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData );

//! Validate/initialize a covariance-like matrix input.
/*!
 * If the input matrix is 0x0, a zero matrix with the required dimensions is created.
 * Otherwise, the function verifies that the input matrix is square and matches the
 * full estimated-parameter dimension.
 * \param covarianceMatrix Input covariance-like matrix.
 * \param numberOfEstimatedParameters Required matrix size based on the parameter set.
 * \return Validated covariance-like matrix.
 */
Eigen::MatrixXd getValidatedCovarianceInput( const Eigen::MatrixXd& covarianceMatrix, const int numberOfEstimatedParameters );

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
 * \param requireAllEntriesToMatch If true, each provided parameter identifier must match at least one parameter block.
 * \return Updated covariance-like matrix.
 */
template< typename InitialStateParameterType = double >
Eigen::MatrixXd addCovarianceDiagonalEntries(
        const Eigen::MatrixXd& covarianceMatrix,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                covarianceDiagonalEntriesPerParameter,
        const bool requireAllEntriesToMatch = true )
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

        const std::vector< std::pair< int, int > > parameterIndices = parameterSet->getIndicesForParameterType( parameterIdentifier );

        if( parameterIndices.size( ) == 0 && requireAllEntriesToMatch )
        {
            throw std::runtime_error(
                    "Error when applying covariance diagonal entries: no parameters found for identifier with parameter type '" +
                    estimatable_parameters::getParameterTypeString( parameterIdentifier.first ) + "' and strings ('" +
                    parameterIdentifier.second.first + "', '" + parameterIdentifier.second.second +
                    "'). Disable strict checking by setting requireAllEntriesToMatch to false." );
        }

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
 * \param requireAllEntriesToMatch If true, each provided parameter identifier must match at least one parameter block.
 * \return Covariance-like matrix with constrained diagonal entries set.
 */
template< typename InitialStateParameterType = double >
Eigen::MatrixXd createCovarianceFromDiagonalEntries(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                covarianceDiagonalEntriesPerParameter,
        const bool requireAllEntriesToMatch = true )
{
    return addCovarianceDiagonalEntries< InitialStateParameterType >(
            Eigen::MatrixXd::Zero( 0, 0 ), parameterSet, covarianceDiagonalEntriesPerParameter, requireAllEntriesToMatch );
}

//! Function to create the consider covariance for reference points (e.g. landmarks) that are held fixed.
/*!
 *  Each ground_station_position consider parameter gets its position variance on the diagonal. Other consider
 *  parameters keep their values from the input matrix.
 *  \param supplementaryData Tracking supplementary data holding the position uncertainties.
 *  \param parameterSet Parameter set that holds the consider parameters.
 *  \param considerCovariance Consider covariance to add to, or a 0x0 matrix to start from zeros.
 *  \return Covariance of the consider parameters, in the order of the consider parameter set.
 */
template< typename InitialStateParameterType = double >
Eigen::MatrixXd createReferencePointConsiderCovariance(
        const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameterSet,
        const Eigen::MatrixXd& considerCovariance = Eigen::MatrixXd::Zero( 0, 0 ) )
{
    const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > considerParameters =
            parameterSet->getConsiderParameters( );
    if( considerParameters == nullptr )
    {
        throw std::runtime_error( "Error when creating reference point consider covariance: parameter set has no consider parameters." );
    }

    const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > > varianceEntries =
            getReferencePointPositionAprioriCovarianceEntries( supplementaryData );
    std::set< estimatable_parameters::EstimatebleParameterIdentifier > referencePointsWithUncertainty;
    for( const auto& entry : varianceEntries )
    {
        referencePointsWithUncertainty.insert( entry.first );
    }
    for( const auto& parameter : considerParameters->getEstimatedVectorParameters( ) )
    {
        const estimatable_parameters::EstimatebleParameterIdentifier parameterName = parameter->getParameterName( );
        if( parameterName.first == estimatable_parameters::ground_station_position &&
            referencePointsWithUncertainty.count( parameterName ) == 0 )
        {
            throw std::runtime_error( "Error when creating reference point consider covariance: reference point " +
                                      parameterName.second.second + " on " + parameterName.second.first +
                                      " is a consider parameter but has no position uncertainty." );
        }
    }

    return addCovarianceDiagonalEntries< InitialStateParameterType >( considerCovariance, considerParameters, varianceEntries, false );
}

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_CREATE_INVERSE_APRIORI_COVARIANCE_H
