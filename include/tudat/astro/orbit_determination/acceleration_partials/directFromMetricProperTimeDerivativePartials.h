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

//! Compute partial of direct-from-metric proper-time derivative.
/*!
 *  \param properTimeDerivative Current proper-time derivative value.
 *  \param parameterSize Size of parameter with respect to which derivative is computed.
 *  \param coordinateVelocity Current coordinate velocity.
 *  \param covariantMetricPerturbation Current covariant metric perturbation.
 *  \param coordinateVelocityPartials Optional partials of coordinate velocity.
 *  \param covariantMetricPerturbationPartials Optional partials of metric perturbation.
 *  \return Row-vector partial of proper-time derivative.
 */
Eigen::Matrix< double, 1, Eigen::Dynamic > calculateDirectProperTimeDerivativePartial(
        double properTimeDerivative,
        int parameterSize,
        const Eigen::Vector3d& coordinateVelocity,
        Eigen::Matrix< double, 4, 4 >& covariantMetricPerturbation,
        const Eigen::Matrix< double, 3, Eigen::Dynamic >& coordinateVelocityPartials =
            Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, 0 ),
        const std::vector< Eigen::Matrix< double, 4, 4 > >& covariantMetricPerturbationPartials =
            std::vector< Eigen::Matrix< double, 4, 4 > >( ) );

//! Partial model for direct-from-metric proper-time rate.
class DirectFromMetricProperTimeDerivativePartial : public RelativisticTimeDerivativePartial
{
public:
    //! Constructor.
    /*!
     *  \param metricPartial Metric-partial model used to provide metric perturbation partials.
     *  \param properTimeStateDerivativeModel Direct proper-time-rate state-derivative model.
     *  \param referencePoint Reference point of propagated proper time.
     */
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

    //! Check whether translational-state partial for a body is non-zero.
    /*!
     *  \param bodyName Name of body for which partial is requested.
     *  \return True if partial is non-zero.
     */
    bool isStateDerivativePartialWrtTranslationalStateNonNull( const std::string& bodyName ) override
    {
        return ( bodyName == integrationReferencePoint_.first ) ||
               metricPartial_->isMetricPartialWrtTranslationalStateNonNull( bodyName );
    }

    //! Retrieve function for partial w.r.t. translational state of a body.
    /*!
     *  \param bodyName Name of body for which partial is requested.
     *  \return Function returning a 1x6 partial matrix.
     */
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

    //! Evaluate partial w.r.t. vector parameter using metric-partial callback.
    /*!
     *  \param metricPartialsFunction Function providing metric perturbation partials.
     *  \param parameterSize Size of requested parameter.
     *  \param coordinateVelocityPartials Partials of coordinate velocity.
     *  \return Row-vector partial of proper-time derivative.
     */
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

    //! Evaluate partial w.r.t. scalar parameter using metric-partial callback.
    /*!
     *  \param metricPartialsFunction Function providing metric perturbation partial.
     *  \param parameterSize Size of requested parameter.
     *  \param coordinateVelocityPartials Partials of coordinate velocity.
     *  \return Row-vector partial of proper-time derivative.
     */
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

    //! Retrieve function for partial w.r.t. scalar parameter.
    /*!
     *  \param parameter Scalar estimatable parameter.
     *  \return Pair of (partial function, parameter size).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    //! Retrieve function for partial w.r.t. vector parameter.
    /*!
     *  \param parameter Vector estimatable parameter.
     *  \return Pair of (partial function, parameter size).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    //! Update cached quantities used in partial computation.
    /*!
     *  \param currentTime Current evaluation time.
     */
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
