/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#ifndef TUDAT_DIRECTFROMMETRICPROPERTIMEDERIVATIVEPARTIALS_H
#define TUDAT_DIRECTFROMMETRICPROPERTIMEDERIVATIVEPARTIALS_H

#include <functional>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/acceleration_partials/relativisticTimeDerivativePartial.h"
#include "tudat/astro/propagators/relativisticTimeStateDerivative.h"
#include "tudat/astro/orbit_determination/metric_partials/metricPartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

Eigen::Matrix< double, 1, Eigen::Dynamic > calculateDirectProperTimeDerivativePartial(
        double properTimeDerivative,
        int parameterSize,
        const Eigen::Vector3d& coordinateVelocity,
        Eigen::Matrix< double, 4, 4 >& covariantMetricPerturbation,
        const Eigen::Matrix< double, 3, Eigen::Dynamic >& coordinateVelocityPartials =
            Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, 0 ),
        const std::vector< Eigen::Matrix< double, 4, 4 > >& covariantMetricPerturbationPartials =
            std::vector< Eigen::Matrix< double, 4, 4 > >( ) );

class DirectFromMetricProperTimeDerivativePartial : public RelativisticTimeDerivativePartial
{
public:
    DirectFromMetricProperTimeDerivativePartial(
            const std::shared_ptr< MetricPartial > metricPartial,
            const std::shared_ptr< DirectProperTimeRateStateDerivative< double, double > >
                properTimeStateDerivativeModel,
            const std::pair< std::string, std::string >& referencePoint ) :
        RelativisticTimeDerivativePartial( referencePoint, propagators::direct_from_metric ),
        metricPartial_( metricPartial ),
        properTimeStateDerivativeModel_( properTimeStateDerivativeModel ),
        metric_( properTimeStateDerivativeModel_->getSpaceTimeMetric( ) )
    { }

    ~DirectFromMetricProperTimeDerivativePartial( ) override = default;

    bool isStateDerivativePartialWrtTranslationalStateNonNull( const std::string& bodyName ) override
    {
        return ( bodyName == integrationReferencePoint_.first ) ||
               metricPartial_->isMetricPartialWrtTranslationalStateNonNull( bodyName );
    }

    std::function< Eigen::Matrix< double, 1, 6 >( ) > getDerivativeFunctionWrtTranslationalStateOfBody(
            const std::string& bodyName ) override
    {
        std::function< Eigen::Matrix< double, 1, 6 >( ) > partialFunction =
                [](){ return Eigen::Matrix< double, 1, 6 >::Zero( ); };

        std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int > metricPartials =
                metricPartial_->getDerivativeFunctionWrtStateOfIntegratedBody(
                        std::make_pair( bodyName, "" ), propagators::translational_state );

        if( metricPartials.second > 0 )
        {
            Eigen::Matrix< double, 3, Eigen::Dynamic > coordinateVelocityPartials;
            if( bodyName == integrationReferencePoint_.first )
            {
                coordinateVelocityPartials = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, 6 );
                coordinateVelocityPartials.block( 0, 3, 3, 3 ).setIdentity( );
            }
            else
            {
                coordinateVelocityPartials = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, 0 );
            }

            partialFunction = std::bind(
                        &DirectFromMetricProperTimeDerivativePartial::calculatePartialWrtVectorParameter,
                        this,
                        metricPartials.first,
                        metricPartials.second,
                        coordinateVelocityPartials );
        }

        return partialFunction;
    }

    Eigen::Matrix< double, 1, Eigen::Dynamic > calculatePartialWrtVectorParameter(
            const std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >& metricPartialsFunction,
            const int parameterSize,
            const Eigen::Matrix< double, 3, Eigen::Dynamic >& coordinateVelocityPartials )
    {
        Eigen::Matrix< double, 4, 4 > currentMetricPerturbation = metric_->getCurrentCovariantMetricPeturbation( );
        return calculateDirectProperTimeDerivativePartial(
                    properTimeStateDerivativeModel_->getCurrentProperTimeDerivative( ),
                    parameterSize,
                    metric_->getCurrentEvaluationState( ).segment( 3, 3 ),
                    currentMetricPerturbation,
                    coordinateVelocityPartials,
                    metricPartialsFunction( ) );
    }

    Eigen::Matrix< double, 1, Eigen::Dynamic > calculatePartialWrtDoubleParameter(
            const std::function< Eigen::Matrix< double, 4, 4 >( ) >& metricPartialsFunction,
            const int parameterSize,
            const Eigen::Matrix< double, 3, Eigen::Dynamic >& coordinateVelocityPartials )
    {
        Eigen::Matrix< double, 4, 4 > currentMetricPerturbation = metric_->getCurrentCovariantMetricPeturbation( );
        return calculateDirectProperTimeDerivativePartial(
                    properTimeStateDerivativeModel_->getCurrentProperTimeDerivative( ),
                    parameterSize,
                    metric_->getCurrentEvaluationState( ).segment( 3, 3 ),
                    currentMetricPerturbation,
                    coordinateVelocityPartials,
                    std::vector< Eigen::Matrix< double, 4, 4 > >{ metricPartialsFunction( ) } );
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    void update( const double currentTime ) override
    {
        if( !( currentTime_ == currentTime ) )
        {
            metricPartial_->update( );
            currentTime_ = currentTime;
        }
    }

protected:

    std::shared_ptr< MetricPartial > metricPartial_;

    std::shared_ptr< DirectProperTimeRateStateDerivative< double, double > > properTimeStateDerivativeModel_;

    std::shared_ptr< relativity::Metric > metric_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_DIRECTFROMMETRICPROPERTIMEDERIVATIVEPARTIALS_H
