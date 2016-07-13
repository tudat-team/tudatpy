#ifndef THRUSTFUNCTIONS_H
#define THRUSTFUNCTIONS_H

namespace tudat
{

namespace propulsion
{


double computeThrustFromSpecificImpulse(
         const double propellantMassRate, const double specificImpulse );

double computePropellantMassRateFromSpecificImpulse(
         const double thrustMagnitude, const double specificImpulse );


}

}

#endif // THRUSTFUNCTIONS_H
