#include "Tudat/Astrodynamics/Propulsion/thrustFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{

namespace propulsion
{


double computeThrustFromSpecificImpulse(
         const double propellantMassRate, const double specificImpulse )
{
    return  physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * specificImpulse * propellantMassRate;

}

double computePropellantMassRateFromSpecificImpulse(
         const double thrustMagnitude, const double specificImpulse )
{
    return  thrustMagnitude /( physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * specificImpulse );

}


} // namespace propulsion


}
