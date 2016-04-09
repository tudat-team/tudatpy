#ifndef CREATESTATEDERIVATIVEPARTIALS_H
#define CREATESTATEDERIVATIVEPARTIALS_H

#include "Tudat/Astrodynamics/OrbitDetermination/stateDerivativePartial.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/SimulationSetup/createAccelerationPartials.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{


template< typename StateScalarType, typename TimeType, typename InitialStateParameterType >
std::map< propagators::IntegratedStateType, StateDerivativePartialsMap > createStateDerivativePartials(
        const std::map< propagators::IntegratedStateType,
        std::vector< boost::shared_ptr< propagators::SingleStateTypeDerivative< StateScalarType, TimeType > > > > stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > parametersToEstimate )
{
    std::map< propagators::IntegratedStateType, StateDerivativePartialsMap > stateDerivativePartials;

    for( typename std::map< propagators::IntegratedStateType,
         std::vector< boost::shared_ptr< propagators::SingleStateTypeDerivative< StateScalarType, TimeType > > > >::const_iterator
         stateDerivativeIterator = stateDerivativeModels.begin( ); stateDerivativeIterator != stateDerivativeModels.end( );
         stateDerivativeIterator++ )
    {

        switch( stateDerivativeIterator->first )
        {
        case propagators::transational_state:
        {
            if( stateDerivativeIterator->second.size( ) > 1 )
            {
                std::cerr<<"Error, cannot yet process multiple separate same type propagators when making partial derivatives of translational state."<<std::endl;
            }
            else
            {
                basic_astrodynamics::AccelerationMap accelerationModelList = boost::dynamic_pointer_cast< propagators::NBodyStateDerivative< StateScalarType, TimeType > >(
                            stateDerivativeIterator->second.at( 0 ) )->getAccelerationsMap( );
                stateDerivativePartials[ propagators::transational_state ] = createAccelerationPartialsMap< InitialStateParameterType >(
                            accelerationModelList, bodyMap, parametersToEstimate );
            }
            break;
        }
        default:
            std::cerr<<"Cannot yet create state derivative partial models for type "<<stateDerivativeIterator->first<<std::endl;
        }
    }
    return stateDerivativePartials;
}

}

}

}
#endif // CREATESTATEDERIVATIVEPARTIALS_H
