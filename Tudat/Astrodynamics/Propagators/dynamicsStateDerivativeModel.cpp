#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/thirdBodyGravityPartial.h"

namespace tudat
{
namespace propagators
{

boost::shared_ptr< acceleration_partials::AccelerationPartial > getAccelerationPartialForBody(
        const orbit_determination::StateDerivativePartialsMap& accelerationPartials,
        const basic_astrodynamics::AvailableAcceleration accelerationType,
        const std::string& bodyExertingAcceleration,
        const std::string& bodyUndergoignAcceleration,
        const std::string& thirdBody )
{
    boost::shared_ptr< acceleration_partials::AccelerationPartial > matchingAccelerationPartial;

    for( unsigned int i = 0; i < accelerationPartials.size( ); i++ )
    {
        for( unsigned int j = 0; j < accelerationPartials.at( i ).size( ); j++ )
        {
            boost::shared_ptr< acceleration_partials::AccelerationPartial > accelerationPartial =
                    boost::dynamic_pointer_cast< acceleration_partials::AccelerationPartial >(
                                   accelerationPartials.at( i ).at( j ) );
            if( accelerationPartial == NULL )
            {
                throw std::runtime_error( "Error when getting acceleration partial, input contained non-acceleration partial" );
            }
            else
            {
                bool partialIdentified = false;
                if( ( accelerationPartial->getAccelerationType( ) == accelerationType ) &&
                    ( accelerationPartial->getAcceleratedBody( ) == bodyExertingAcceleration ) &&
                        ( accelerationPartial->getAcceleratingBody( ) == bodyUndergoignAcceleration ) )
                {
                    partialIdentified = true;
                }

                if( basic_astrodynamics::isAccelerationFromThirdBody( accelerationType ) &&
                        getThirdBodyFromAccelerationPartial( accelerationPartial ) == thirdBody )
                {
                    partialIdentified = true;
                }

                if( partialIdentified )
                {
                    if( matchingAccelerationPartial != NULL )
                    {
                        throw std::runtime_error( "Error when getting acceleration partial, found multiple matching accelerations." );
                    }
                    else
                    {
                        matchingAccelerationPartial = accelerationPartial;
                    }
                }
            }
        }
    }

    return matchingAccelerationPartial;
}

} // namespace propagators

} // namespace tudat
