#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"

#if( BUILD_WITH_ESTIMATION_TOOLS )
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/thirdBodyGravityPartial.h"
#endif

namespace tudat
{
namespace propagators
{

#if( BUILD_WITH_ESTIMATION_TOOLS )
//! Function to retrieve specific acceleration partial object from list of state derivative partials
std::shared_ptr< acceleration_partials::AccelerationPartial > getAccelerationPartialForBody(
        const orbit_determination::StateDerivativePartialsMap& accelerationPartials,
        const basic_astrodynamics::AvailableAcceleration accelerationType,
        const std::string& bodyExertingAcceleration,
        const std::string& bodyUndergoignAcceleration,
        const std::string& centralBody )
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
                bool partialIdentified = false;
                if( ( accelerationPartial->getAccelerationType( ) == accelerationType ) &&
                        ( accelerationPartial->getAcceleratedBody( ) == bodyExertingAcceleration ) &&
                        ( accelerationPartial->getAcceleratingBody( ) == bodyUndergoignAcceleration ) )
                {
                    if( !basic_astrodynamics::isAccelerationFromThirdBody( accelerationType  ) )
                    {
                        partialIdentified = true;
                    }
                    else if( getCentralBodyNameFromThirdBodyAccelerationPartial( accelerationPartial ) == centralBody )
                    {
                        partialIdentified = true;
                    }
                }

                // Set output; check if multiple models are identified.
                if( partialIdentified )
                {
                    if( matchingAccelerationPartial != nullptr )
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
#endif

template class DynamicsStateDerivativeModel< double, double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class DynamicsStateDerivativeModel< Time, double >;
template class DynamicsStateDerivativeModel< double, long double >;
template class DynamicsStateDerivativeModel< Time, long double >;
#endif


} // namespace propagators

} // namespace tudat
