/*    Copyright (c) 2010-2019,
 *    Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
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

class SchwarzschildMetricPartial : public MetricPartial
{
public:
    using MetricPartial::getParameterPartialFunction;

    SchwarzschildMetricPartial(
            const std::shared_ptr< relativity::HarmonicSchwarzschildMetric > schwarzschildMetric,
            const std::pair< std::string, std::string >& evaluationPointName );

    ~SchwarzschildMetricPartial( ) override = default;

    bool isMetricPartialWrtTranslationalStateNonNull( const std::string& bodyName ) override;

    std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
        getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) override;

    Eigen::Matrix< double, 4, 4 > wrtScaledTime( ) override;

    std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    void update( ) override;

    std::vector< Eigen::Matrix< double, 4, 4 > > getPartialWrtReferencePointPosition( );

    Eigen::Matrix< double, 4, 4 > getPartialWrtCentralBodyGravitationalParameter( );

    Eigen::Matrix< double, 4, 4 > getPartialWrtPpnParameterGamma( );

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
