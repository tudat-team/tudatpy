/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/astro/orbit_determination/stateDerivativePartial.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/simulation/estimation_setup/createTorquePartials.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a set of state derivative partial objects.
/*!
 *  Function to create a set of state derivative partial objects for any propagated state types.
 *  \param stateDerivativeModels List of state derivative models, ordered by state type (key)
 *  \param bodies List of boy objects storing environment models of simulation
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
        const simulation_setup::SystemOfBodies& bodies,
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
                            accelerationModelList, bodies, parametersToEstimate );
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
                            torqueModelList, bodies, parametersToEstimate );
            }
            break;
        }
        case propagators::body_mass_state:
        {
            std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > initialMassParameters =
                    getListOfMassStateParametersToEstimate( parametersToEstimate );
            orbit_determination::StateDerivativePartialsMap massPartials;
            for( unsigned int i = 0; i < initialMassParameters.size( ); i++ )
            {
                massPartials.push_back(
                            std::vector< std::shared_ptr< orbit_determination::StateDerivativePartial > >( ) );
            }


            stateDerivativePartials[ propagators::body_mass_state ] = massPartials;
            std::cerr<<"Warning, mass state partials implicitly set to zero - depending on non-conservative force settings"<<
                       " and thrust guidance/mass rate model, this may provide biased results for variational equations "<<
                       "(implicit assumption: mass influence nothing, and nothing influences mass)"<<std::endl;
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
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate );

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< long double, double >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, double > > > >
        stateDerivativeModels,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
        parametersToEstimate );
extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< double, Time >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, Time > > > >
        stateDerivativeModels,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate );
extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< long double, Time >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, Time > > > >
        stateDerivativeModels,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
        parametersToEstimate );
#endif

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATESTATEDERIVATIVEPARTIALS_H
