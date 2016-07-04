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

//! Function to create a set of state derivative partial objects.
/*!
 *  Function to create a set of state derivative partial objects for any propagated state types.
 *  \param stateDerivativeModels List of state derivative models, ordered by state type (key)
 *  \param bodyMap List of boy objects storing environment models of simulation
 *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and
 *  values.
 *  return List partials of state derivative models from. The key is the type of dynamics for which partials are taken,
 *  the values are StateDerivativePartialsMap (see StateDerivativePartialsMap definition for details).
 */
template< typename StateScalarType, typename TimeType, typename InitialStateParameterType >
std::map< propagators::IntegratedStateType, StateDerivativePartialsMap > createStateDerivativePartials(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< boost::shared_ptr< propagators::SingleStateTypeDerivative< StateScalarType, TimeType > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >
        parametersToEstimate )
{
    std::map< propagators::IntegratedStateType, StateDerivativePartialsMap > stateDerivativePartials;

    // Iterate over all state types
    for( typename std::unordered_map< propagators::IntegratedStateType,
         std::vector< boost::shared_ptr< propagators::SingleStateTypeDerivative< StateScalarType, TimeType > > > >::
         const_iterator stateDerivativeIterator = stateDerivativeModels.begin( );
         stateDerivativeIterator != stateDerivativeModels.end( );
         stateDerivativeIterator++ )
    {
        // Identify state types
        switch( stateDerivativeIterator->first )
        {
        case propagators::transational_state:
        {
            if( stateDerivativeIterator->second.size( ) > 1 )
            {
                throw std::runtime_error(
                            "Error, cannot yet process multiple separate same type propagators when making partial derivatives of translational state." );
            }
            else
            {
                // Retrieve acceleration models and create partials
                basic_astrodynamics::AccelerationMap accelerationModelList =
                        boost::dynamic_pointer_cast< propagators::NBodyStateDerivative< StateScalarType, TimeType > >(
                            stateDerivativeIterator->second.at( 0 ) )->getFullAccelerationsMap( );
                stateDerivativePartials[ propagators::transational_state ] =
                        createAccelerationPartialsMap< InitialStateParameterType >(
                            accelerationModelList, bodyMap, parametersToEstimate );
            }
            break;
        }
        default:
            std::string errorMessage = "Cannot yet create state derivative partial models for type " +
                    boost::lexical_cast< std::string >( stateDerivativeIterator->first );
            throw std::runtime_error( errorMessage );
        }
    }
    return stateDerivativePartials;
}

}

}

}
#endif // CREATESTATEDERIVATIVEPARTIALS_H
