/*    Copyright (c) 2010-2019,
 *    Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#ifndef TUDAT_SCHWARZSCHILDMETRICCHRISTOFFELSYMBOLPARTIALS_H
#define TUDAT_SCHWARZSCHILDMETRICCHRISTOFFELSYMBOLPARTIALS_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/metric_partials/christoffelSymbolPartial.h"
#include "tudat/astro/relativity/schwarzschildMetric.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

class SchwarzschildMetricChristoffelSymbolPartial : public ChristoffelSymbolPartial
{
public:
    SchwarzschildMetricChristoffelSymbolPartial(
            const std::shared_ptr< relativity::HarmonicSchwarzschildMetric > metricObject,
            const std::string& centralBody );

    ~SchwarzschildMetricChristoffelSymbolPartial( ) override = default;

    std::pair< std::function< std::vector< std::vector< Eigen::Matrix< double, 4, 4 > > >( ) >, int >
        getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) override;

    std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int >
        getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int >
        getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    void update( const double currentTime ) override;

private:
    std::vector< std::vector< Eigen::Matrix4d > > currentChristoffelSymbolsPositionPartials_;

    std::shared_ptr< relativity::HarmonicSchwarzschildMetric > metricObject_;

    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

    std::string centralBody_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_SCHWARZSCHILDMETRICCHRISTOFFELSYMBOLPARTIALS_H
