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


#ifndef TUDAT_SCHWARZSCHILDMETRICPARTIAL_H
#define TUDAT_SCHWARZSCHILDMETRICPARTIAL_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/metric_partials/metricPartial.h"
#include "tudat/astro/relativity/schwarzschildMetric.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

//! Metric-partial model for the harmonic Schwarzschild metric.
class SchwarzschildMetricPartial : public MetricPartial
{
public:
    using MetricPartial::getParameterPartialFunction;

    //! Constructor.
    /*!
     *  \param schwarzschildMetric Harmonic Schwarzschild metric object.
     *  \param evaluationPointName Pair identifying the metric evaluation point (body, reference point).
     */
    SchwarzschildMetricPartial(
            const std::shared_ptr< relativity::HarmonicSchwarzschildMetric > schwarzschildMetric,
            const std::pair< std::string, std::string >& evaluationPointName );

    ~SchwarzschildMetricPartial( ) override = default;

    //! Check whether partial w.r.t. translational state of a body is non-zero.
    /*!
     *  \param bodyName Name of body for which partial is requested.
     *  \return True if partial is non-zero.
     */
    bool isMetricPartialWrtTranslationalStateNonNull( const std::string& bodyName ) override;

    //! Retrieve function for partials w.r.t. integrated-body state.
    /*!
     *  \param stateReferencePoint Body/reference-point identifier.
     *  \param integratedStateType Integrated state type.
     *  \return Pair of (partial function, number of output columns).
     */
    std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
        getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) override;

    //! Evaluate metric partial w.r.t. scaled time.
    /*!
     *  \return Metric partial matrix w.r.t. scaled time.
     */
    Eigen::Matrix< double, 4, 4 > wrtScaledTime( ) override;

    //! Retrieve function for partial w.r.t. scalar parameter.
    /*!
     *  \param parameter Scalar estimatable parameter.
     *  \return Pair of (partial function, parameter size).
     */
    std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    //! Update cached partial terms.
    void update( ) override;

    //! Get metric partials w.r.t. evaluation-point position.
    /*!
     *  \return Vector of metric partial matrices w.r.t. position components.
     */
    std::vector< Eigen::Matrix< double, 4, 4 > > getPartialWrtReferencePointPosition( );

    //! Get metric partial w.r.t. central-body gravitational parameter.
    /*!
     *  \return Metric partial matrix.
     */
    Eigen::Matrix< double, 4, 4 > getPartialWrtCentralBodyGravitationalParameter( );

    //! Get metric partial w.r.t. PPN parameter gamma.
    /*!
     *  \return Metric partial matrix.
     */
    Eigen::Matrix< double, 4, 4 > getPartialWrtPpnParameterGamma( );

    //! Get metric partial w.r.t. PPN parameter beta.
    /*!
     *  \return Metric partial matrix.
     */
    Eigen::Matrix< double, 4, 4 > getPartialWrtPpnParameterBeta( );

protected:
    std::shared_ptr< relativity::HarmonicSchwarzschildMetric > schwarzschildMetric_;

    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

    std::string centralBodyName_;

    std::pair< std::string, std::string > evaluationPointName_;

    std::vector< Eigen::Matrix< double, 4, 4 > > currentPartialWrtEvaluationPointPosition_;

    double currentDistance_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_SCHWARZSCHILDMETRICPARTIAL_H
