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


#ifndef TUDAT_METRICPARTIAL_H
#define TUDAT_METRICPARTIAL_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/relativity/metric.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

//! Rotate metric partials from body-centered frame to local ground-station frame.
/*!
 *  \param metricPartialWrtBodyPosition Metric partials w.r.t. body position in inertial/body-centered frame.
 *  \param rotationToLocalFrame Rotation matrix from body-centered frame to local frame.
 *  \return Metric partials w.r.t. ground-station position in local frame.
 */
std::vector< Eigen::Matrix< double, 4, 4 > > getMetricPartialWrtGroundStationPosition(
        const std::vector< Eigen::Matrix< double, 4, 4 > >& metricPartialWrtBodyPosition,
        const Eigen::Matrix3d& rotationToLocalFrame );

//! Base interface for metric partial derivatives.
class MetricPartial
{
public:
    MetricPartial( ) = default;
    virtual ~MetricPartial( ) = default;

    //! Check whether metric partial w.r.t. a body's translational state is non-zero.
    /*!
     *  \param bodyName Name of body for which translational-state partial is requested.
     *  \return True if the partial is non-zero.
     */
    virtual bool isMetricPartialWrtTranslationalStateNonNull( const std::string& bodyName ) = 0;

    //! Evaluate metric partials w.r.t. state of an integrated body.
    /*!
     *  \param stateReferencePoint Pair containing body and optional reference point identifier.
     *  \param integratedStateType Integrated state type for which partial is requested.
     *  \return Vector of metric partial matrices.
     */
    std::vector< Eigen::Matrix< double, 4, 4 > > wrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        return getDerivativeFunctionWrtStateOfIntegratedBody( stateReferencePoint, integratedStateType ).first( );
    }

    //! Evaluate metric partial w.r.t. scaled time variable.
    /*!
     *  \return Metric partial matrix w.r.t. scaled time.
     */
    virtual Eigen::Matrix< double, 4, 4 > wrtScaledTime( ) = 0;

    //! Retrieve function that computes metric partials w.r.t. integrated-body state.
    /*!
     *  \param stateReferencePoint Pair containing body and optional reference point identifier.
     *  \param integratedStateType Integrated state type for which partial is requested.
     *  \return Pair of (partial function, number of output columns).
     */
    virtual std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
        getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) = 0;

    //! Retrieve function for partial w.r.t. scalar estimatable parameter.
    /*!
     *  \param parameter Scalar estimatable parameter.
     *  \return Pair of (partial function, parameter size). Default implementation returns an empty function.
     */
    virtual std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int > getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< Eigen::Matrix< double, 4, 4 >( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Retrieve function for partial w.r.t. vector estimatable parameter.
    /*!
     *  \param parameter Vector estimatable parameter.
     *  \return Pair of (partial function, parameter size). Default implementation returns an empty function.
     */
    virtual std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int > getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Update cached quantities used in metric partial evaluation.
    virtual void update( ) = 0;
};

//! Compute Christoffel symbols from metric partials and metric object at the current state.
/*!
 *  \param currentState Current translational state at metric evaluation point.
 *  \param currentTime Current evaluation time.
 *  \param metricEvaluationPoint Body/reference-point identifier at which metric is evaluated.
 *  \param metric Metric object used to compute inverse metric and related terms.
 *  \param metricPartial Metric-partial model providing derivative terms.
 *  \return Vector of Christoffel-symbol matrices.
 */
std::vector< Eigen::Matrix< double, 4, 4 > > calculateCurrentChristoffelSymbolsFromMetricPartials(
        const Eigen::Matrix< double, 6, 1 >& currentState,
        const double currentTime,
        const std::pair< std::string, std::string >& metricEvaluationPoint,
        const std::shared_ptr< relativity::Metric > metric,
        const std::shared_ptr< MetricPartial > metricPartial );

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_METRICPARTIAL_H
