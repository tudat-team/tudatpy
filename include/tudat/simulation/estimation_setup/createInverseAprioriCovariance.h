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
std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > aprioriUncertaintyEntry(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier,
        const double uncertainty );

std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > aprioriUncertaintyEntry(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier,
        const Eigen::VectorXd& uncertaintyValues );

Eigen::MatrixXd getValidatedInverseAprioriCovarianceInput(
        const Eigen::MatrixXd& inverseAprioriCovariance,
        const int numberOfEstimatedParameters );

Eigen::VectorXd getUncertaintyValuesForParameter(
        const Eigen::VectorXd& uncertaintyValues,
        const int parameterSize,
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier );

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
                parameterSet->getIndicesForParameterIdentifier( parameterIdentifier );

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
