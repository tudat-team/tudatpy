/*    Copyright (c) 2010-2019,
 *    Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
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

std::vector< Eigen::Matrix< double, 4, 4 > > getMetricPartialWrtGroundStationPosition(
        const std::vector< Eigen::Matrix< double, 4, 4 > >& metricPartialWrtBodyPosition,
        const Eigen::Matrix3d& rotationToLocalFrame );

class MetricPartial
{
public:
    MetricPartial( ) = default;
    virtual ~MetricPartial( ) = default;

    virtual bool isMetricPartialWrtTranslationalStateNonNull( const std::string& bodyName ) = 0;

    std::vector< Eigen::Matrix< double, 4, 4 > > wrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        return getDerivativeFunctionWrtStateOfIntegratedBody( stateReferencePoint, integratedStateType ).first( );
    }

    virtual Eigen::Matrix< double, 4, 4 > wrtScaledTime( ) = 0;

    virtual std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
        getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) = 0;

    virtual std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int > getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< Eigen::Matrix< double, 4, 4 >( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    virtual std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int > getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    virtual void update( ) = 0;
};

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
