#include "tudat/astro/orbit_determination/metric_partials/schwarzschildMetricPartial.h"

#include <stdexcept>

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

SchwarzschildMetricPartial::SchwarzschildMetricPartial(
        const std::shared_ptr< relativity::HarmonicSchwarzschildMetric > schwarzschildMetric,
        const std::pair< std::string, std::string >& evaluationPointName ):
    schwarzschildMetric_( schwarzschildMetric ),
    evaluationPointName_( evaluationPointName )
{
    centralBodyName_ = schwarzschildMetric_->getCentralBodyName( );
    ppnParameterSet_ = schwarzschildMetric_->getPpnParameterSet( );
    currentPartialWrtEvaluationPointPosition_.resize( 6 );
    for( int i = 3; i < 6; ++i )
    {
        currentPartialWrtEvaluationPointPosition_[ i ].setZero( );
    }
}

bool SchwarzschildMetricPartial::isMetricPartialWrtTranslationalStateNonNull( const std::string& bodyName )
{
    return bodyName == evaluationPointName_.first;
}

std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
SchwarzschildMetricPartial::getDerivativeFunctionWrtStateOfIntegratedBody(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int > partialFunction =
            std::make_pair( std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >( ), 0 );

    if( integratedStateType == propagators::translational_state )
    {
        if( !stateReferencePoint.second.empty( ) )
        {
            throw std::runtime_error(
                        "Schwarzschild metric partial: translational state reference may not include body-fixed point." );
        }

        if( stateReferencePoint.first == centralBodyName_ )
        {
            throw std::runtime_error(
                        "Schwarzschild metric partial: cannot take derivative w.r.t. central-body state; it is assumed static." );
        }
        else if( stateReferencePoint.first == evaluationPointName_.first )
        {
            partialFunction = std::make_pair(
                        std::bind( &SchwarzschildMetricPartial::getPartialWrtReferencePointPosition, this ), 6 );
        }
    }

    return partialFunction;
}

Eigen::Matrix< double, 4, 4 > SchwarzschildMetricPartial::wrtScaledTime( )
{
    return Eigen::Matrix4d::Zero( );
}

std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int >
SchwarzschildMetricPartial::getParameterPartialFunction(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int > partialFunctionPair =
            std::make_pair( std::function< Eigen::Matrix< double, 4, 4 >( ) >( ), 0 );

    switch( parameter->getParameterName( ).first )
    {
    case estimatable_parameters::gravitational_parameter:
        if( parameter->getParameterName( ).second.first == centralBodyName_ )
        {
            partialFunctionPair = std::make_pair(
                        std::bind( &SchwarzschildMetricPartial::getPartialWrtCentralBodyGravitationalParameter, this ), 1 );
        }
        break;
    case estimatable_parameters::ppn_parameter_gamma:
        partialFunctionPair = std::make_pair(
                    std::bind( &SchwarzschildMetricPartial::getPartialWrtPpnParameterGamma, this ), 1 );
        break;
    case estimatable_parameters::ppn_parameter_beta:
        partialFunctionPair = std::make_pair(
                    std::bind( &SchwarzschildMetricPartial::getPartialWrtPpnParameterBeta, this ), 1 );
        break;
    default:
        break;
    }

    return partialFunctionPair;
}

void SchwarzschildMetricPartial::update( )
{
    const Eigen::Vector3d currentEvaluationPosition =
            schwarzschildMetric_->getCurrentEvaluationState( ).segment( 0, 3 );
    currentDistance_ = currentEvaluationPosition.norm( );
    const double currentDistanceSquared = currentDistance_ * currentDistance_;

    Eigen::Matrix< double, 4, 4 > basePartial =
            -schwarzschildMetric_->getCurrentFirstOrderCovariantMetricContributions( ) / currentDistanceSquared;

    const double baseSecondOrderTimeTimeTerm =
            4.0 * ppnParameterSet_->getParameterBeta( ) * physical_constants::INVERSE_QUARTIC_SPEED_OF_LIGHT *
            std::pow( schwarzschildMetric_->getCurrentCentralGravitationalParameter( ), 2 ) /
            std::pow( currentDistanceSquared, 2 );
    basePartial( 0, 0 ) += baseSecondOrderTimeTimeTerm;

    if( schwarzschildMetric_->getIncludeSecondPostNewtonianOrder( ) )
    {
        throw std::runtime_error(
                    "Schwarzschild metric partial: second-order terms not yet implemented." );
    }

    for( int i = 0; i < 3; ++i )
    {
        currentPartialWrtEvaluationPointPosition_[ i ] = basePartial * currentEvaluationPosition( i );
    }
}

std::vector< Eigen::Matrix< double, 4, 4 > > SchwarzschildMetricPartial::getPartialWrtReferencePointPosition( )
{
    return currentPartialWrtEvaluationPointPosition_;
}

Eigen::Matrix< double, 4, 4 > SchwarzschildMetricPartial::getPartialWrtCentralBodyGravitationalParameter( )
{
    return ( schwarzschildMetric_->getCurrentFirstOrderCovariantMetricContributions( ) +
             2.0 * schwarzschildMetric_->getCurrentSecondOrderCovariantMetricContributions( ) ) /
            schwarzschildMetric_->getCurrentCentralGravitationalParameter( );
}

Eigen::Matrix< double, 4, 4 > SchwarzschildMetricPartial::getPartialWrtPpnParameterGamma( )
{
    const double partialTerm = 2.0 * schwarzschildMetric_->getCurrentCentralGravitationalParameter( ) / currentDistance_ *
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
    Eigen::Matrix< double, 4, 4 > metricPartial = Eigen::Matrix< double, 4, 4 >::Zero( );
    for( int i = 1; i < 4; ++i )
    {
        metricPartial( i, i ) = partialTerm;
    }
    return metricPartial;
}

Eigen::Matrix< double, 4, 4 > SchwarzschildMetricPartial::getPartialWrtPpnParameterBeta( )
{
    const double partialTerm = schwarzschildMetric_->getCurrentCentralGravitationalParameter( ) / currentDistance_ *
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
    Eigen::Matrix< double, 4, 4 > metricPartial = Eigen::Matrix< double, 4, 4 >::Zero( );
    metricPartial( 0, 0 ) = -2.0 * partialTerm * partialTerm;
    return metricPartial;
}

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

