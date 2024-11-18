#include "tudat/astro/propagators/dynamicsStateDerivativeModel.h"

#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
#include "tudat/astro/orbit_determination/acceleration_partials/thirdBodyGravityPartial.h"
#endif

namespace tudat
{
namespace propagators
{

#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
//! Function to retrieve specific acceleration partial object from list of state derivative partials
std::shared_ptr< acceleration_partials::AccelerationPartial > getAccelerationPartialForBody(
        const orbit_determination::StateDerivativePartialsMap& accelerationPartials,
        const basic_astrodynamics::AvailableAcceleration accelerationType,
        const std::string& bodyExertingAcceleration,
        const std::string& bodyUndergoignAcceleration )
{
    std::shared_ptr< acceleration_partials::AccelerationPartial > matchingAccelerationPartial;

    // Iterate over all partials
    for( unsigned int i = 0; i < accelerationPartials.size( ); i++ )
    {
        for( unsigned int j = 0; j < accelerationPartials.at( i ).size( ); j++ )
        {
            // Check id state derivative type is acceleration
            std::shared_ptr< acceleration_partials::AccelerationPartial > accelerationPartial =
                    std::dynamic_pointer_cast< acceleration_partials::AccelerationPartial >(
                        accelerationPartials.at( i ).at( j ) );
            if( accelerationPartial != nullptr )
            {
                // Compare current partial against required data
                if( ( accelerationPartial->getAccelerationType( ) == accelerationType ) &&
                        ( accelerationPartial->getAcceleratedBody( ) == bodyExertingAcceleration ) &&
                        ( accelerationPartial->getAcceleratingBody( ) == bodyUndergoignAcceleration ) )
                {
                    if( matchingAccelerationPartial != nullptr )
                    {
                        throw std::runtime_error( "Error when getting acceleration partial exerted by " + bodyExertingAcceleration +
                                                  " on " + bodyUndergoignAcceleration +
                                                  " of type " + basic_astrodynamics::getAccelerationModelName( accelerationType ) +
                                                  ", found multiple matching accelerations." );
                    }
                    else
                    {
                        matchingAccelerationPartial = accelerationPartial;
                    }
                }
            }
        }
    }

    if( matchingAccelerationPartial == nullptr )
    {
        if( basic_astrodynamics::isAccelerationDirectGravitational( accelerationType ) )
        {
            matchingAccelerationPartial = getAccelerationPartialForBody(
                        accelerationPartials,
                        basic_astrodynamics::getAssociatedThirdBodyAcceleration( accelerationType ),
                        bodyExertingAcceleration, bodyUndergoignAcceleration );
        }
    }

    if( matchingAccelerationPartial == nullptr )
    {
        throw std::runtime_error( "Error when getting acceleration partial exerted by " + bodyExertingAcceleration +
                                  " on " + bodyUndergoignAcceleration +
                                  " of type " + basic_astrodynamics::getAccelerationModelName( accelerationType ) +
                                  ", found no matching accelerations." );
    }

    return matchingAccelerationPartial;
}
#endif

//template class DynamicsStateDerivativeModel< double, double >;


} // namespace propagators

} // namespace tudat
