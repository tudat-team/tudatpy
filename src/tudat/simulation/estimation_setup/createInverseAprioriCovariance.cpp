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

namespace
{

Eigen::MatrixXd getValidatedInverseAprioriCovarianceInput(
        const Eigen::MatrixXd& inverseAprioriCovariance,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >& parameterSet )
{
    if( parameterSet == nullptr )
    {
        throw std::runtime_error( "Error when creating/updating inverse apriori covariance: parameter_set is null." );
    }

    const int numberOfEstimatedParameters = parameterSet->getEstimatedParameterSetSize( );
    if( inverseAprioriCovariance.rows( ) == 0 && inverseAprioriCovariance.cols( ) == 0 )
    {
        return Eigen::MatrixXd::Zero( numberOfEstimatedParameters, numberOfEstimatedParameters );
    }

    if( inverseAprioriCovariance.rows( ) != numberOfEstimatedParameters ||
        inverseAprioriCovariance.cols( ) != numberOfEstimatedParameters )
    {
        throw std::runtime_error( "Error when creating/updating inverse apriori covariance: provided matrix has size " +
                                  std::to_string( inverseAprioriCovariance.rows( ) ) + "x" +
                                  std::to_string( inverseAprioriCovariance.cols( ) ) + ", expected " +
                                  std::to_string( numberOfEstimatedParameters ) + "x" +
                                  std::to_string( numberOfEstimatedParameters ) + "." );
    }

    return inverseAprioriCovariance;
}

Eigen::VectorXd getUncertaintyValuesForParameter(
        const Eigen::VectorXd& uncertaintyValues,
        const int parameterSize,
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier )
{
    Eigen::VectorXd currentUncertaintyValues;

    if( uncertaintyValues.size( ) == 1 )
    {
        currentUncertaintyValues = Eigen::VectorXd::Constant( parameterSize, uncertaintyValues( 0 ) );
    }
    else if( uncertaintyValues.size( ) == parameterSize )
    {
        currentUncertaintyValues = uncertaintyValues;
    }
    else
    {
        throw std::runtime_error(
                "Error when applying apriori uncertainty entries for parameter type '" +
                estimatable_parameters::getParameterTypeString( parameterIdentifier.first ) + "' with identifiers ('" +
                parameterIdentifier.second.first + "', '" + parameterIdentifier.second.second +
                "'): uncertainty size is incompatible with parameter size. "
                "Received uncertainty vector of size " +
                std::to_string( uncertaintyValues.size( ) ) + " for parameter size " + std::to_string( parameterSize ) +
                ". Use a scalar uncertainty or a vector matching the parameter size." );
    }

    if( currentUncertaintyValues.size( ) == 0 )
    {
        throw std::runtime_error(
                "Error when applying apriori uncertainty entries for parameter type '" +
                estimatable_parameters::getParameterTypeString( parameterIdentifier.first ) + "' with identifiers ('" +
                parameterIdentifier.second.first + "', '" + parameterIdentifier.second.second + "'): uncertainty is empty." );
    }

    for( int i = 0; i < currentUncertaintyValues.size( ); i++ )
    {
        if( !( currentUncertaintyValues( i ) > 0.0 ) )
        {
            throw std::runtime_error(
                    "Error when applying apriori uncertainty entries for parameter type '" +
                    estimatable_parameters::getParameterTypeString( parameterIdentifier.first ) + "' with identifiers ('" +
                    parameterIdentifier.second.first + "', '" + parameterIdentifier.second.second +
                    "'): uncertainty entries must be strictly positive." );
        }
    }

    return currentUncertaintyValues;
}

}  // namespace

std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > aprioriUncertaintyEntry(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier,
        const double uncertainty )
{
    return std::make_pair( parameterIdentifier, Eigen::VectorXd::Constant( 1, uncertainty ) );
}

std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > aprioriUncertaintyEntry(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier,
        const Eigen::VectorXd& uncertaintyValues )
{
    return std::make_pair( parameterIdentifier, uncertaintyValues );
}

Eigen::MatrixXd addInverseAprioriCovarianceEntries(
        const Eigen::MatrixXd& inverseAprioriCovariance,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                aprioriUncertaintyPerParameter )
{
    Eigen::MatrixXd updatedInverseAprioriCovariance =
            getValidatedInverseAprioriCovarianceInput( inverseAprioriCovariance, parameterSet );

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

Eigen::MatrixXd createInverseAprioriCovariance(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                aprioriUncertaintyPerParameter )
{
    return addInverseAprioriCovarianceEntries( Eigen::MatrixXd::Zero( 0, 0 ), parameterSet, aprioriUncertaintyPerParameter );
}

}  // namespace simulation_setup
}  // namespace tudat
