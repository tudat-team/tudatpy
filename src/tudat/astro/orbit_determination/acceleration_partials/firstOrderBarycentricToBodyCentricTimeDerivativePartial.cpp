#include "tudat/astro/orbit_determination/acceleration_partials/firstOrderBarycentricToBodyCentricTimeDerivativePartial.h"

#include <algorithm>
#include <stdexcept>

#include "tudat/basics/utilities.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

Eigen::Matrix< double, 1, 3 > getFirstOrderBarycentricToBodyCentricTimeDerivativeWrtBodyPosition(
        const double currentBodyDistance,
        const double externalBodyGravitationalParameter,
        const Eigen::Vector3d& externalBodyRelativePosition )
{
    return physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
            externalBodyRelativePosition.transpose( ) *
            externalBodyGravitationalParameter /
            ( currentBodyDistance * currentBodyDistance * currentBodyDistance );
}

Eigen::Matrix< double, 1, 3 > getFirstOrderBarycentricToBodyCentricTimeDerivativeWrtBodyPosition(
        const double externalBodyGravitationalParameter,
        const Eigen::Vector3d& externalBodyRelativePosition )
{
    return getFirstOrderBarycentricToBodyCentricTimeDerivativeWrtBodyPosition(
                externalBodyRelativePosition.norm( ),
                externalBodyGravitationalParameter,
                externalBodyRelativePosition );
}

FirstOrderBarycentricToBodyCentricTimeDerivativePartial::
FirstOrderBarycentricToBodyCentricTimeDerivativePartial(
        const std::shared_ptr< FirstOrderBarycentricToBodyCentricTimeStateDerivative< double, double > >
        stateDerivativeModel ):
    RelativisticTimeDerivativePartial(
        std::make_pair( stateDerivativeModel->getCentralBody( ), "" ),
        propagators::first_order_barycentric_to_bodycentric ),
    externalBodies_( stateDerivativeModel->getExternalBodies( ) ),
    stateDerivativeModel_( stateDerivativeModel )
{
    currentPartialsWrtExternalBodyPositions_.resize( externalBodies_.size( ) );
    for( unsigned int i = 0; i < externalBodies_.size( ); ++i )
    {
        externalBodyIndexMap_[ externalBodies_.at( i ) ] = static_cast< int >( i );
    }
}

std::function< Eigen::Matrix< double, 1, 6 >( ) >
FirstOrderBarycentricToBodyCentricTimeDerivativePartial::
getDerivativeFunctionWrtTranslationalStateOfBody( const std::string& bodyName )
{
    if( externalBodyIndexMap_.count( bodyName ) > 0 )
    {
        const int bodyIndex = externalBodyIndexMap_.at( bodyName );
        return std::bind(
                    &FirstOrderBarycentricToBodyCentricTimeDerivativePartial::wrtTranslationalStateOfBody,
                    this,
                    bodyIndex );
    }
    else if( bodyName == integrationReferencePoint_.first )
    {
        return std::bind(
                    &FirstOrderBarycentricToBodyCentricTimeDerivativePartial::getPartialWrtCentralBodyState,
                    this );
    }
    else
    {
        return [](){ return Eigen::Matrix< double, 1, 6 >::Zero( ); };
    }
}

bool FirstOrderBarycentricToBodyCentricTimeDerivativePartial::
isStateDerivativePartialWrtTranslationalStateNonNull( const std::string& bodyName )
{
    if( bodyName == integrationReferencePoint_.first )
    {
        return true;
    }
    return std::find( externalBodies_.begin( ), externalBodies_.end( ), bodyName ) != externalBodies_.end( );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
FirstOrderBarycentricToBodyCentricTimeDerivativePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    switch( parameter->getParameterName( ).first )
    {
    case estimatable_parameters::gravitational_parameter:
    {
        auto bodyFindIterator = std::find(
                    externalBodies_.begin( ), externalBodies_.end( ), parameter->getParameterName( ).second.first );
        if( bodyFindIterator != externalBodies_.end( ) )
        {
            const int vectorIndex = static_cast< int >(
                        std::distance( externalBodies_.begin( ), bodyFindIterator ) );
            partialFunction = [ this, vectorIndex ]( Eigen::MatrixXd& partial )
            {
                partial = getPartialWrtBodyGravitationalParameter( vectorIndex );
            };
            numberOfColumns = 1;
        }
        break;
    }
    default:
        break;
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

void FirstOrderBarycentricToBodyCentricTimeDerivativePartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        const Eigen::Vector3d currentCentralBodyPosition =
                stateDerivativeModel_->getCurrentCentralBodyState( ).segment( 0, 3 );

        currentPartialWrtCentralBodyPosition_ = Eigen::Matrix< double, 1, 3 >::Zero( );

        for( unsigned int i = 0; i < externalBodies_.size( ); ++i )
        {
            currentPartialsWrtExternalBodyPositions_[ i ] =
                    getFirstOrderBarycentricToBodyCentricTimeDerivativeWrtBodyPosition(
                        stateDerivativeModel_->getCurrentExternalBodyDistance( i ),
                        stateDerivativeModel_->getCurrentExternalBodyGravitationalParameter( i ),
                        stateDerivativeModel_->getCurrentExternalBodyStates( i ).segment( 0, 3 ) -
                        currentCentralBodyPosition );

            currentPartialWrtCentralBodyPosition_ -= currentPartialsWrtExternalBodyPositions_.at( i );
        }

        currentTime_ = currentTime;
    }
}

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat
