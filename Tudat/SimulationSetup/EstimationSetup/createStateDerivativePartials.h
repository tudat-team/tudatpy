/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATESTATEDERIVATIVEPARTIALS_H
#define TUDAT_CREATESTATEDERIVATIVEPARTIALS_H

#include "Tudat/Astrodynamics/OrbitDetermination/stateDerivativePartial.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"
#include "Tudat/SimulationSetup/EstimationSetup/createAccelerationPartials.h"
#include "Tudat/SimulationSetup/EstimationSetup/createTorquePartials.h"

namespace tudat
{

namespace simulation_setup
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
template< typename StateScalarType, typename TimeType >
std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< StateScalarType, TimeType > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > >
        parametersToEstimate )
{
    std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > stateDerivativePartials;

    // Iterate over all state types
    for( typename std::unordered_map< propagators::IntegratedStateType,
         std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< StateScalarType, TimeType > > > >::
         const_iterator stateDerivativeIterator = stateDerivativeModels.begin( );
         stateDerivativeIterator != stateDerivativeModels.end( );
         stateDerivativeIterator++ )
    {
        // Identify state types
        switch( stateDerivativeIterator->first )
        {
        case propagators::translational_state:
        {
            if( stateDerivativeIterator->second.size( ) > 1 )
            {
                throw std::runtime_error(
                            "Error, cannot yet process multiple separate same type propagators when making partial "
                            "derivatives of translational state." );
            }
            else
            {
                // Retrieve acceleration models and create partials
                basic_astrodynamics::AccelerationMap accelerationModelList =
                        std::dynamic_pointer_cast< propagators::NBodyStateDerivative< StateScalarType, TimeType > >(
                            stateDerivativeIterator->second.at( 0 ) )->getFullAccelerationsMap( );
                stateDerivativePartials[ propagators::translational_state ] =
                        createAccelerationPartialsMap< StateScalarType >(
                            accelerationModelList, bodyMap, parametersToEstimate );
            }
            break;
        }
        case propagators::rotational_state:
        {
            if( stateDerivativeIterator->second.size( ) > 1 )
            {
                throw std::runtime_error(
                            "Error, cannot yet process multiple separate same type propagators when making partial derivatives of rotational state." );
            }
            else
            {
                // Retrieve acceleration models and create partials
                basic_astrodynamics::TorqueModelMap torqueModelList =
                        std::dynamic_pointer_cast< propagators::RotationalMotionStateDerivative< StateScalarType, TimeType > >(
                            stateDerivativeIterator->second.at( 0 ) )->getTorquesMap( );
                stateDerivativePartials[ propagators::rotational_state ] =
                        createTorquePartialsMap< StateScalarType >(
                            torqueModelList, bodyMap, parametersToEstimate );
            }
            break;
        }
        default:
            std::string errorMessage = "Cannot yet create state derivative partial models for type " +
                    std::to_string( stateDerivativeIterator->first );
            throw std::runtime_error( errorMessage );
        }
    }

    return stateDerivativePartials;
}

extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< double, double >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, double > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate );

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< long double, double >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, double > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
        parametersToEstimate );
extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< double, Time >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, Time > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate );
extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< long double, Time >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, Time > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
        parametersToEstimate );
#endif

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATESTATEDERIVATIVEPARTIALS_H
