/* Copyright (c) 2010-2019,
 *    Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SOLARSYSTEMMETRICPARTIALS_H
#define TUDAT_SOLARSYSTEMMETRICPARTIALS_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/metric_partials/metricPartial.h"
#include "tudat/astro/relativity/solarSystemMetric.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"

namespace tudat
{

class SphericalHarmonicPartialWrapper
{
public:
    SphericalHarmonicPartialWrapper(
            const std::shared_ptr< acceleration_partials::SphericalHarmonicsGravityPartial >
            sphericalHarmonicGravityPartial,
            const std::shared_ptr< SphericalHarmonicWrapper > sphericalHarmonicWrapper ):
        sphericalHarmonicGravityPartial_( sphericalHarmonicGravityPartial ),
        sphericalHarmonicWrapper_( sphericalHarmonicWrapper )
    { }

private:
    std::shared_ptr< acceleration_partials::SphericalHarmonicsGravityPartial >
            sphericalHarmonicGravityPartial_;

    std::shared_ptr< SphericalHarmonicWrapper > sphericalHarmonicWrapper_;
};

namespace orbit_determination
{

namespace partial_derivatives
{

std::vector< Eigen::Matrix< double, 4, 4 > > calculateSolarSystemMetricPartialFromPotentialPartials(
        const Eigen::VectorXd& scalarPotentialPartials,
        const Eigen::Matrix< double, 3, Eigen::Dynamic >& vectorPotentialPartials,
        const double totalScalarPotential,
        const std::shared_ptr< relativity::PPNParameterSet > currentPpnParameterSet );

Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialFromPotentialPartials(
        const double scalarPotentialPartials,
        const Eigen::Vector3d& vectorPotentialPartials,
        const double totalScalarPotential,
        const std::shared_ptr< relativity::PPNParameterSet > currentPpnParameterSet );

std::pair< std::vector< Eigen::Matrix< double, 4, 4 > >, std::vector< Eigen::Matrix< double, 4, 4 > > >
calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
        const std::shared_ptr< SphericalHarmonicWrapper > shWrapper,
        const int bodyIndex,
        const std::map< int, std::pair< int, int > >& coefficientBlockIndices );

std::pair< std::vector< Eigen::Matrix< double, 4, 4 > >, std::vector< Eigen::Matrix< double, 4, 4 > > >
calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
        const int bodyIndex,
        const int maximumDegree,
        const int maximumOrder,
        const int minimumDegree = 2,
        const int minimumOrder = 0 );

std::pair< std::vector< Eigen::Matrix< double, 4, 4 > >, std::vector< Eigen::Matrix< double, 4, 4 > > >
calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
        const int bodyIndex,
        const std::map< int, std::pair< int, int > >& coefficientBlockIndices );

Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialWrtPpnParameterGamma(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric );

Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialWrtPpnParameterBeta(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric );

Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialWrtBodyGravitationalParameter(
        const int bodyIndex,
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric );

class SolarSystemMetricPartial : public MetricPartial
{
public:
    SolarSystemMetricPartial(
            const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
            const std::pair< std::string, std::string >& evaluationPointName,
            const std::function< Eigen::Quaterniond( const double ) > rotationFunctionFromBodyFixedFrame =
                [] ( const double ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const std::map< int, std::shared_ptr< SphericalHarmonicPartialWrapper > >& sphericalHarmonicPartialWrappers =
                std::map< int, std::shared_ptr< SphericalHarmonicPartialWrapper > >( ) );

    ~SolarSystemMetricPartial( ) override = default;

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

    void update( ) override;

    Eigen::Matrix< double, 4, 4 > wrtScaledTime( ) override;

    std::vector< Eigen::Matrix< double, 4, 4 > > getCurrentPartialWrtBodyState( const int bodyIndex );

    std::vector< Eigen::Matrix< double, 4, 4 > > getPartialWrtReferencePointPosition( );

    Eigen::Matrix< double, 4, 4 > getDoubleParameterPartial(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier );

    std::vector< Eigen::Matrix< double, 4, 4 > > getVectorParameterPartial(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier );

    std::vector< Eigen::Matrix< double, 4, 4 > > calculateSolarSystemMetricPartialWrtGroundStationPosition( );

private:
    void updateParameterPartials( );

    std::vector< Eigen::Matrix< double, 4, 4 > > getUpdatedSphericalHarmonicPartials(
            const int bodyIndex,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterId,
            const std::map< int, std::pair< int, int > > coefficientBlockIndices );

    std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric_;

    std::pair< std::string, std::string > evaluationPointName_;

    std::map< int, std::shared_ptr< SphericalHarmonicPartialWrapper > > sphericalHarmonicPartialWrappers_;

    std::vector< std::string > bodyList_;

    std::function< Eigen::Quaterniond( const double ) > rotationFunctionFromBodyFixedFrame_;

    std::vector< Eigen::Matrix< double, 4, 4 > > currentPartialWrtEvaluationPointPosition_;

    std::vector< std::vector< Eigen::Matrix< double, 4, 4 > > > currentPartialsWrtMetricBodyStates_;

    std::map< estimatable_parameters::EstimatebleParameterIdentifier,
              std::function< Eigen::Matrix< double, 4, 4 >( ) > > doubleParameterPartialFunctions_;

    std::map< estimatable_parameters::EstimatebleParameterIdentifier,
              std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) > > vectorParameterPartialFunctions_;

    std::map< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::Matrix< double, 4, 4 > >
            currentDoubleParameterPartials_;

    std::map< estimatable_parameters::EstimatebleParameterIdentifier, std::vector< Eigen::Matrix< double, 4, 4 > > >
            currentVectorParameterPartials_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_SOLARSYSTEMMETRICPARTIALS_H
