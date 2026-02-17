/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/orbit_determination/acceleration_partials/directFromMetricProperTimeDerivativePartials.h"

#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

Eigen::Matrix< double, 1, Eigen::Dynamic > calculateDirectProperTimeDerivativePartial(
        const double properTimeDerivative,
        const int parameterSize,
        const Eigen::Vector3d& coordinateVelocity,
        Eigen::Matrix< double, 4, 4 >& covariantMetricPerturbation,
        const Eigen::Matrix< double, 3, Eigen::Dynamic >& coordinateVelocityPartials,
        const std::vector< Eigen::Matrix< double, 4, 4 > >& covariantMetricPerturbationPartials )
{
    const Eigen::Vector4d coordinateFourVelocity =
            ( Eigen::Vector4d( ) << physical_constants::SPEED_OF_LIGHT, coordinateVelocity ).finished( );

    Eigen::Matrix< double, Eigen::Dynamic, 1 > partialVector =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( parameterSize, 1 );

    if( coordinateVelocityPartials.cols( ) > 0 )
    {
        if( coordinateVelocityPartials.cols( ) != parameterSize )
        {
            throw std::runtime_error(
                        "Error when calculating direct proper-time derivative partials: "
                        "coordinate-velocity partials dimension mismatch." );
        }

        Eigen::Matrix< double, 4, Eigen::Dynamic > coordinateFourVelocityPartials =
                Eigen::Matrix< double, 4, Eigen::Dynamic >::Zero( 4, parameterSize );
        coordinateFourVelocityPartials.block( 1, 0, 3, parameterSize ) = coordinateVelocityPartials;

        partialVector -= 2.0 * coordinateVelocity.transpose( ) * coordinateVelocityPartials;
        partialVector -= 2.0 * coordinateFourVelocity.transpose( ) *
                covariantMetricPerturbation * coordinateFourVelocityPartials;
    }

    if( !covariantMetricPerturbationPartials.empty( ) )
    {
        if( static_cast< int >( covariantMetricPerturbationPartials.size( ) ) != parameterSize )
        {
            throw std::runtime_error(
                        "Error when calculating direct proper-time derivative partials: "
                        "metric-perturbation partials dimension mismatch." );
        }

        for( int i = 0; i < parameterSize; ++i )
        {
            partialVector( i ) -=
                    coordinateFourVelocity.transpose( ) *
                    covariantMetricPerturbationPartials.at( i ) * coordinateFourVelocity;
        }
    }

    return -partialVector.transpose( ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT /
            ( 2.0 * ( properTimeDerivative + 1.0 ) );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
DirectFromMetricProperTimeDerivativePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int > metricPartialFunctionPair =
            metricPartial_->getParameterPartialFunction( parameter );

    if( metricPartialFunctionPair.second > 0 )
    {
        partialFunction = [ this, metricPartialFunction = metricPartialFunctionPair.first ]( Eigen::MatrixXd& partial )
        {
            partial = calculatePartialWrtDoubleParameter(
                        metricPartialFunction,
                        1,
                        Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, 0 ) );
        };
        numberOfColumns = 1;
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
DirectFromMetricProperTimeDerivativePartial::getParameterPartialFunction(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int > metricPartialFunctionPair =
            metricPartial_->getParameterPartialFunction( parameter );

    if( metricPartialFunctionPair.second > 0 )
    {
        const int parameterSize = parameter->getParameterSize( );
        partialFunction = [ this,
                            metricPartialFunction = metricPartialFunctionPair.first,
                            parameterSize ]( Eigen::MatrixXd& partial )
        {
            partial = calculatePartialWrtVectorParameter(
                        metricPartialFunction,
                        parameterSize,
                        Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, parameterSize ) );
        };
        numberOfColumns = parameter->getParameterSize( );
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat
