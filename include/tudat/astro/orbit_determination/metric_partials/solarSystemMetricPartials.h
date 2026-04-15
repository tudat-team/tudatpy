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
    //! Constructor.
    /*!
     *  \param sphericalHarmonicGravityPartial Spherical-harmonic gravity partial object.
     *  \param sphericalHarmonicWrapper Wrapper containing metric-side spherical-harmonic utilities.
     */
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

//! Convert scalar/vector potential partials to metric partials (vector-valued parameter case).
/*!
 *  \param scalarPotentialPartials Scalar-potential partials.
 *  \param vectorPotentialPartials Vector-potential partials.
 *  \param totalScalarPotential Current total scalar potential.
 *  \param currentPpnParameterSet Current PPN parameter set.
 *  \return Metric partial matrices.
 */
std::vector< Eigen::Matrix< double, 4, 4 > > calculateSolarSystemMetricPartialFromPotentialPartials(
        const Eigen::VectorXd& scalarPotentialPartials,
        const Eigen::Matrix< double, 3, Eigen::Dynamic >& vectorPotentialPartials,
        const double totalScalarPotential,
        const std::shared_ptr< relativity::PPNParameterSet > currentPpnParameterSet );

//! Convert scalar/vector potential partials to metric partial (scalar parameter case).
/*!
 *  \param scalarPotentialPartials Scalar-potential partial.
 *  \param vectorPotentialPartials Vector-potential partial.
 *  \param totalScalarPotential Current total scalar potential.
 *  \param currentPpnParameterSet Current PPN parameter set.
 *  \return Metric partial matrix.
 */
Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialFromPotentialPartials(
        const double scalarPotentialPartials,
        const Eigen::Vector3d& vectorPotentialPartials,
        const double totalScalarPotential,
        const std::shared_ptr< relativity::PPNParameterSet > currentPpnParameterSet );

//! Compute metric partials w.r.t. spherical-harmonic coefficients using wrapper and block map.
/*!
 *  \param shWrapper Spherical-harmonic wrapper for requested body.
 *  \param bodyIndex Index of body in solar-system metric body list.
 *  \param coefficientBlockIndices Requested coefficient blocks.
 *  \return Pair of (cosine partials, sine partials).
 */
std::pair< std::vector< Eigen::Matrix< double, 4, 4 > >, std::vector< Eigen::Matrix< double, 4, 4 > > >
calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
        const std::shared_ptr< SphericalHarmonicWrapper > shWrapper,
        const int bodyIndex,
        const std::map< int, std::pair< int, int > >& coefficientBlockIndices );

//! Compute metric partials w.r.t. spherical-harmonic coefficients for degree/order range.
/*!
 *  \param solarSystemMetric Solar-system metric object.
 *  \param bodyIndex Index of body in solar-system metric body list.
 *  \param maximumDegree Maximum degree.
 *  \param maximumOrder Maximum order.
 *  \param minimumDegree Minimum degree.
 *  \param minimumOrder Minimum order.
 *  \return Pair of (cosine partials, sine partials).
 */
std::pair< std::vector< Eigen::Matrix< double, 4, 4 > >, std::vector< Eigen::Matrix< double, 4, 4 > > >
calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
        const int bodyIndex,
        const int maximumDegree,
        const int maximumOrder,
        const int minimumDegree = 2,
        const int minimumOrder = 0 );

//! Compute metric partials w.r.t. spherical-harmonic coefficients for block map.
/*!
 *  \param solarSystemMetric Solar-system metric object.
 *  \param bodyIndex Index of body in solar-system metric body list.
 *  \param coefficientBlockIndices Requested coefficient blocks.
 *  \return Pair of (cosine partials, sine partials).
 */
std::pair< std::vector< Eigen::Matrix< double, 4, 4 > >, std::vector< Eigen::Matrix< double, 4, 4 > > >
calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
        const int bodyIndex,
        const std::map< int, std::pair< int, int > >& coefficientBlockIndices );

//! Compute metric partial w.r.t. PPN gamma.
/*!
 *  \param solarSystemMetric Solar-system metric object.
 *  \return Metric partial matrix.
 */
Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialWrtPpnParameterGamma(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric );

//! Compute metric partial w.r.t. PPN beta.
/*!
 *  \param solarSystemMetric Solar-system metric object.
 *  \return Metric partial matrix.
 */
Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialWrtPpnParameterBeta(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric );

//! Compute metric partial w.r.t. one body's gravitational parameter.
/*!
 *  \param bodyIndex Index of body in solar-system metric body list.
 *  \param solarSystemMetric Solar-system metric object.
 *  \return Metric partial matrix.
 */
Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialWrtBodyGravitationalParameter(
        const int bodyIndex,
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric );

//! Metric-partial model for solar-system post-Newtonian metric.
class SolarSystemMetricPartial : public MetricPartial
{
public:
    //! Constructor.
    /*!
     *  \param solarSystemMetric Solar-system metric object.
     *  \param evaluationPointName Pair identifying metric evaluation point (body, reference point).
     *  \param rotationFunctionFromBodyFixedFrame Rotation from body-fixed to inertial frame at current time.
     *  \param sphericalHarmonicPartialWrappers Optional map of spherical-harmonic partial wrappers by body index.
     */
    SolarSystemMetricPartial(
            const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
            const std::pair< std::string, std::string >& evaluationPointName,
            const std::function< Eigen::Quaterniond( const double ) > rotationFunctionFromBodyFixedFrame =
                [] ( const double ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const std::map< int, std::shared_ptr< SphericalHarmonicPartialWrapper > >& sphericalHarmonicPartialWrappers =
                std::map< int, std::shared_ptr< SphericalHarmonicPartialWrapper > >( ) );

    ~SolarSystemMetricPartial( ) override = default;

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

    //! Update all cached metric partials.
    void update( ) override;

    //! Evaluate metric partial w.r.t. scaled time.
    /*!
     *  \return Metric partial matrix w.r.t. scaled time.
     */
    Eigen::Matrix< double, 4, 4 > wrtScaledTime( ) override;

    //! Get current partials w.r.t. one body's translational state.
    /*!
     *  \param bodyIndex Index of body in solar-system metric body list.
     *  \return Vector of metric partial matrices.
     */
    std::vector< Eigen::Matrix< double, 4, 4 > > getCurrentPartialWrtBodyState( const int bodyIndex );

    //! Get current partials w.r.t. evaluation-point position.
    /*!
     *  \return Vector of metric partial matrices.
     */
    std::vector< Eigen::Matrix< double, 4, 4 > > getPartialWrtReferencePointPosition( );

    //! Get current scalar-parameter metric partial from cache.
    /*!
     *  \param parameterIdentifier Parameter identifier.
     *  \return Metric partial matrix.
     */
    Eigen::Matrix< double, 4, 4 > getDoubleParameterPartial(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier );

    //! Get current vector-parameter metric partial from cache.
    /*!
     *  \param parameterIdentifier Parameter identifier.
     *  \return Vector of metric partial matrices.
     */
    std::vector< Eigen::Matrix< double, 4, 4 > > getVectorParameterPartial(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier );

    //! Compute metric partials w.r.t. ground-station position at evaluation point.
    /*!
     *  \return Vector of metric partial matrices.
     */
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
