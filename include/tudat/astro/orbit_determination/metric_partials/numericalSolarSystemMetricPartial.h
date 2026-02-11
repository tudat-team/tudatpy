/*    Copyright (c) 2010-2019,
 *    Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
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

class NumericalSolarSystemMetricPartial : public MetricPartial
{
public:
    NumericalSolarSystemMetricPartial(
            const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
            const std::pair< std::string, std::string >& evaluationPointName );

    ~NumericalSolarSystemMetricPartial( ) override = default;

    bool isMetricPartialWrtTranslationalStateNonNull( const std::string& bodyName ) override;

    std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
        getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) override;

    std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

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
