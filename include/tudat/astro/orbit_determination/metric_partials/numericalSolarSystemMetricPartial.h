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


#ifndef TUDAT_NUMERICALSOLARSYSTEMMETRICPARTIAL_H
#define TUDAT_NUMERICALSOLARSYSTEMMETRICPARTIAL_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/relativity/solarSystemMetric.h"
#include "tudat/astro/orbit_determination/metric_partials/metricPartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

//! Numerical partial model for the solar-system metric.
class NumericalSolarSystemMetricPartial : public MetricPartial
{
public:
    //! Constructor.
    /*!
     *  \param solarSystemMetric Solar-system metric object.
     *  \param evaluationPointName Pair identifying the metric evaluation point (body, reference point).
     */
    NumericalSolarSystemMetricPartial(
            const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
            const std::pair< std::string, std::string >& evaluationPointName );

    ~NumericalSolarSystemMetricPartial( ) override = default;

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

    //! Retrieve function for partial w.r.t. scalar parameter.
    /*!
     *  \param parameter Scalar estimatable parameter.
     *  \return Pair of (partial function, parameter size).
     */
    std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    //! Retrieve function for partial w.r.t. vector parameter.
    /*!
     *  \param parameter Vector estimatable parameter.
     *  \return Pair of (partial function, parameter size).
     */
    std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    //! Update partial model (no-op for current numerical implementation).
    void update( ) override { }

private:
    std::vector< Eigen::Matrix< double, 4, 4 > > getPartialWrtReferencePointPosition( );
    std::vector< Eigen::Matrix< double, 4, 4 > > getCurrentPartialWrtBodyState( const int bodyIndex );

    std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric_;

    std::pair< std::string, std::string > evaluationPointName_;

    std::vector< std::string > bodyList_;

    std::vector< Eigen::Matrix< double, 4, 4 > > currentPartialWrtEvaluationPointPosition_;

    std::vector< std::vector< Eigen::Matrix< double, 4, 4 > > > currentPartialsWrtMetricBodyStates_;

    std::map< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::Matrix< double, 4, 4 > >
            currentDoubleParameterPartials_;

    std::map< estimatable_parameters::EstimatebleParameterIdentifier, std::vector< Eigen::Matrix< double, 4, 4 > > >
            currentVectorParameterPartials_;

    std::vector< std::pair< estimatable_parameters::EstimatableParameter< double >, double > > requestedDoubleParameters_;

    std::vector< std::pair< estimatable_parameters::EstimatableParameter< Eigen::VectorXd >, Eigen::VectorXd > >
            requestedVectorParameters_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_NUMERICALSOLARSYSTEMMETRICPARTIAL_H
