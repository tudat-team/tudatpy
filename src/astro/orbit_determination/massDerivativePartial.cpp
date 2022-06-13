
#include "tudat/astro/orbit_determination/massDerivativePartial.h"


namespace tudat
{

namespace orbit_determination
{

FromThrustMassRatePartial::FromThrustMassRatePartial( const std::string& body,
                                                      const std::shared_ptr< propulsion::FromThrustMassRateModel > massRateModel ):
    MassRatePartial( body, basic_astrodynamics::from_thrust_mass_rate_model )
{
    thrustAccelerations_ = massRateModel->getThrustAccelerations( );

    for( unsigned int i = 0; i < thrustAccelerations_.size( ); i++ )
    {
        std::vector< std::shared_ptr< system_models::EngineModel > > thrustSources =
                thrustAccelerations_.at( i )->getThrustSources( );

        std::vector< std::shared_ptr< system_models::EngineModel > > currentMassDependentThrustSources;
        for( unsigned int j = 0; j < thrustSources.size( ); j++ )
        {
            if( !thrustSources.at( j )->getThrustMagnitudeWrapper( )->modelIsForceBased( ) )
            {
                currentMassDependentThrustSources.push_back( thrustSources.at( j ) );
            }
            engineModelList_[ thrustSources.at( j )->getEngineName( ) ] = thrustSources.at( j );
        }

        if( currentMassDependentThrustSources.size( ) > 0 )
        {
            accelerationBasedThrustSources_[ i ] = currentMassDependentThrustSources;

        }
    }
}

bool FromThrustMassRatePartial::isMassRatePartialWrtMassNonZero( )
{
    if( accelerationBasedThrustSources_.size( ) != 0 )
    {
        return true;
    }
    else
    {
        return false;
    }
}

void FromThrustMassRatePartial::wrtMassOfBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix )
{
    for( auto it : accelerationBasedThrustSources_ )
    {
        for( unsigned int i = 0; i < it.second.size( ); i++ )
        {
            partialMatrix( 0, 0 ) += it.second.at( i )->getCurrentMassRate( ) /
                    thrustAccelerations_.at( i )->getCurrentBodyMass( );
        }
    }
}

void FromThrustMassRatePartial::wrtEngineSpecificImpulse(
        Eigen::MatrixXd& partialMatrix,
        const std::string& engineName )
{
    std::shared_ptr< system_models::EngineModel > engineModel = engineModelList_.at( engineName );
    partialMatrix( 0, 0 ) += engineModel->getCurrentMassRate( currentBodyMass_ ) /
            engineModel->getThrustMagnitudeWrapper( )->getCurrentSpecificImpulse( );
}

void FromThrustMassRatePartial::wrtEngineThrustMagnitude(
        Eigen::MatrixXd& partialMatrix,
        const std::string& engineName )
{
    std::shared_ptr< system_models::EngineModel > engineModel = engineModelList_.at( engineName );
    partialMatrix( 0, 0 ) += -1.0 /
            ( engineModel->getThrustMagnitudeWrapper( )->getCurrentSpecificImpulse( ) *
              physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int > FromThrustMassRatePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction = std::make_pair( nullptr, 0 );

    if( parameter->getParameterName( ).first == estimatable_parameters::constant_thrust_magnitude_parameter &&
            parameter->getParameterName( ).second.first == body_ &&
            engineModelList_.count( parameter->getParameterName( ).second.second ) != 0 )
    {
        partialFunction = std::make_pair(
                    std::bind( &FromThrustMassRatePartial::wrtEngineThrustMagnitude, this,
                               std::placeholders::_1,
                               parameter->getParameterName( ).second.second ), 1 );

    }

    else if( parameter->getParameterName( ).first == estimatable_parameters::constant_specific_impulse &&
             parameter->getParameterName( ).second.first == body_ &&
             engineModelList_.count( parameter->getParameterName( ).second.second ) != 0 )
    {
        partialFunction = std::make_pair(
                    std::bind( &FromThrustMassRatePartial::wrtEngineSpecificImpulse, this,
                               std::placeholders::_1,
                               parameter->getParameterName( ).second.second ), 1 );

    }

    return partialFunction;
}


} // namespace orbit_determination

} // namespace tudat
