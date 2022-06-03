#include "tudat/astro/propulsion/thrustAccelerationModel.h"

namespace tudat
{

namespace propulsion
{

void ThrustAcceleration::updateMembers( const double currentTime )
{
    // Check if update is needed
    if( !( currentTime_ == currentTime ) )
    {
        currentMassRate_ = 0.0;
        currentAcceleration_.setZero( );
        thrustDirectionCalculator_->update( currentTime );

        for( unsigned int i = 0; i < thrustSources_.size( ); i++ )
        {
            thrustSources_.at( i )->updateEngineModel( currentTime );
            if( !saveThrustContributions_ )
            {
                currentMassRate_ -= thrustSources_.at( i )->getCurrentMassRate( );
                currentAcceleration_ += ( thrustSources_.at( i )->getCurrentThrust( ) / bodyMassFunction_( ) )*
                        thrustDirectionCalculator_->getInertialThrustDirection( thrustSources_.at( i ) ) ;
            }
            else
            {
                currentMassRateContributions_[ i ] = thrustSources_.at( i )->getCurrentMassRate( );
                currentMassRate_ -= currentMassRateContributions_[ i ];
                currentThrustAccelerationContributions_[ i ] = ( thrustSources_.at( i )->getCurrentThrust( ) / bodyMassFunction_( ) )*
                        thrustDirectionCalculator_->getInertialThrustDirection( thrustSources_.at( i ) );
                currentAcceleration_ += currentThrustAccelerationContributions_[ i ];
            }
        }

        // Reset current time.
        currentTime_ = currentTime;

    }
}

Eigen::Vector3d ThrustAcceleration::getCurrentThrustAccelerationContribution(
        const unsigned int index )
{
    if( !saveThrustContributions_ )
    {
        throw std::runtime_error( "Error when getting single thrust acceleration contribution, separate contributions not saved" );
    }
    else if( currentTime_ != currentTime_ )
    {
        throw std::runtime_error( "Error when getting single thrust acceleration contribution, thrust model is not updated" );
    }
    else if( index >= currentThrustAccelerationContributions_.size( ) )
    {
        throw std::runtime_error( "Error when getting single thrust acceleration contribution, requested index is out of bounds" );
    }
    else
    {
        return currentThrustAccelerationContributions_.at( index );
    }
}

} // namespace propulsion

} // namespace tudat
