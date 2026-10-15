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

std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
getCameraPointingInverseAprioriCovarianceEntries(
        const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData )
{
    std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > > entries;
    for( const std::shared_ptr< data::TrackingSupplementaryData >& currentData : supplementaryData )
    {
        for( const std::shared_ptr< data::InstrumentSupplementaryData >& instrument : currentData->getInstrumentSupplementaryData( ) )
        {
            const auto camera = std::dynamic_pointer_cast< data::CameraInstrumentSupplementaryData >( instrument );
            if( camera == nullptr || camera->getPointingUncertainty( ).array( ).isNaN( ).all( ) )
            {
                continue;
            }
            const Eigen::Vector3d& uncertainty = camera->getPointingUncertainty( );
            if( !uncertainty.array( ).isFinite( ).all( ) || ( uncertainty.array( ) <= 0.0 ).any( ) )
            {
                throw std::runtime_error( "Error when creating pointing a-priori for camera " + camera->getCameraId( ) +
                                          ": uncertainties must be finite and positive." );
            }
            entries.push_back( std::make_pair( estimatable_parameters::EstimatebleParameterIdentifier(
                                                       estimatable_parameters::camera_pointing_correction,
                                                       std::make_pair( currentData->getBodyName( ), camera->getCameraId( ) ) ),
                                               Eigen::VectorXd( uncertainty.array( ).square( ).inverse( ) ) ) );
        }
    }
    return entries;
}

std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
getReferencePointPositionInverseAprioriCovarianceEntries(
        const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData )
{
    std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > > entries;
    for( const std::shared_ptr< data::TrackingSupplementaryData >& currentData : supplementaryData )
    {
        const std::shared_ptr< data::ReferencePointSupplementaryData >& referencePoint = currentData->getReferencePointSupplementaryData( );
        if( referencePoint == nullptr || referencePoint->getPositionUncertainty( ).array( ).isNaN( ).all( ) )
        {
            continue;
        }
        const Eigen::Vector3d& uncertainty = referencePoint->getPositionUncertainty( );
        if( !uncertainty.array( ).isFinite( ).all( ) || ( uncertainty.array( ) <= 0.0 ).any( ) )
        {
            throw std::runtime_error( "Error when creating position a-priori for reference point " + currentData->getReferencePointName( ) +
                                      " on " + currentData->getBodyName( ) + ": uncertainties must be finite and positive." );
        }
        entries.push_back( std::make_pair( estimatable_parameters::EstimatebleParameterIdentifier(
                                                   estimatable_parameters::ground_station_position,
                                                   std::make_pair( currentData->getBodyName( ), currentData->getReferencePointName( ) ) ),
                                           Eigen::VectorXd( uncertainty.array( ).square( ).inverse( ) ) ) );
    }
    return entries;
}

//! Validate matrix dimensions or create a correctly sized zero matrix from a 0x0 input.
Eigen::MatrixXd getValidatedCovarianceInput( const Eigen::MatrixXd& covarianceMatrix, const int numberOfEstimatedParameters )
{
    if( covarianceMatrix.rows( ) == 0 && covarianceMatrix.cols( ) == 0 )
    {
        return Eigen::MatrixXd::Zero( numberOfEstimatedParameters, numberOfEstimatedParameters );
    }

    if( covarianceMatrix.rows( ) != numberOfEstimatedParameters || covarianceMatrix.cols( ) != numberOfEstimatedParameters )
    {
        throw std::runtime_error( "Error when creating/updating covariance matrix entries: provided matrix has size " +
                                  std::to_string( covarianceMatrix.rows( ) ) + "x" + std::to_string( covarianceMatrix.cols( ) ) +
                                  ", expected " + std::to_string( numberOfEstimatedParameters ) + "x" +
                                  std::to_string( numberOfEstimatedParameters ) + "." );
    }

    return covarianceMatrix;
}

//! Convert scalar/vector covariance diagonal input to per-component values and validate non-negativity.
Eigen::VectorXd getCovarianceDiagonalValuesForParameter( const Eigen::VectorXd& covarianceDiagonalValues,
                                                         const int parameterSize,
                                                         const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier )
{
    Eigen::VectorXd currentCovarianceDiagonalValues;

    if( covarianceDiagonalValues.size( ) == 0 )
    {
        throw std::runtime_error( "Error when applying covariance diagonal entries for parameter type '" +
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
        throw std::runtime_error( "Error when applying covariance diagonal entries for parameter type '" +
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
            throw std::runtime_error( "Error when applying covariance diagonal entries for parameter type '" +
                                      estimatable_parameters::getParameterTypeString( parameterIdentifier.first ) +
                                      "' with identifiers ('" + parameterIdentifier.second.first + "', '" +
                                      parameterIdentifier.second.second + "'): entries must be non-negative." );
        }
    }

    return currentCovarianceDiagonalValues;
}

}  // namespace simulation_setup
}  // namespace tudat
