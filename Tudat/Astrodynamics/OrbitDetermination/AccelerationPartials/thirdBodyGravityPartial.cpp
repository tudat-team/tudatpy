#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/thirdBodyGravityPartial.h"

namespace tudat
{

namespace acceleration_partials
{


std::string getThirdBodyFromAccelerationPartial(
        const boost::shared_ptr< AccelerationPartial > accelerationPartial )
{
    std::string thirdBody;
    if( !basic_astrodynamics::isAccelerationFromThirdBody( accelerationPartial->getAccelerationType( ) ) )
    {
        throw std::runtime_error( "Error, requested third body from acceleration partial, but input is incompatible." );
    }
    else
    {
        if( accelerationPartial->getAccelerationType( ) == basic_astrodynamics::third_body_point_mass_gravity )
        {
            thirdBody = boost::dynamic_pointer_cast< ThirdBodyGravityPartial< CentralGravitationPartial > >(
                        accelerationPartial )->getCentralBodyName( );
        }
        else if( accelerationPartial->getAccelerationType( ) == basic_astrodynamics::third_body_spherical_harmonic_gravity )
        {
            thirdBody = boost::dynamic_pointer_cast< ThirdBodyGravityPartial< SphericalHarmonicsGravityPartial > >(
                        accelerationPartial )->getCentralBodyName( );
        }
        else if( accelerationPartial->getAccelerationType( ) == basic_astrodynamics::third_body_mutual_spherical_harmonic_gravity )
        {
            thirdBody = boost::dynamic_pointer_cast< ThirdBodyGravityPartial< MutualSphericalHarmonicsGravityPartial > >(
                        accelerationPartial )->getCentralBodyName( );
        }
        else
        {
            throw std::runtime_error( "Error, requested third body from acceleration partial, input third-body type not recognized." );
        }
    }

    return thirdBody;
}


} // namespace acceleration_partials

} // namespace tudat
