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
#include "tudat/astro/orbit_determination/massDerivativePartial.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"
#include "tudat/astro/propagators/bodyMassStateDerivative.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/simulation/estimation_setup/createTorquePartials.h"

namespace tudat
{

namespace simulation_setup
{

template< typename InitialStateParameterType = double >
std::shared_ptr< orbit_determination::MassRatePartial > createAnalyticalMassRatePartial(
        std::shared_ptr< basic_astrodynamics::MassRateModel > massRateModel,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > bodyWithMassRate,
        const simulation_setup::SystemOfBodies& bodies = simulation_setup::SystemOfBodies( ),
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >
        parametersToEstimate =
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >( ) )
{
    using namespace gravitation;
    using namespace basic_astrodynamics;
    using namespace electromagnetism;
    using namespace aerodynamics;
    using namespace acceleration_partials;
    using namespace propulsion;

    std::shared_ptr< orbit_determination::MassRatePartial > massRatePartial;

    // Identify current massRate model type
    AvailableMassRateModels massRateType = getMassRateModelType( massRateModel );
    switch( massRateType )
    {
    case from_thrust_mass_rate_model:
    {
        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< FromThrustMassRateModel >( massRateModel ) == nullptr )
        {
            throw std::runtime_error( "Mass rate class type does not match mass rate type (from_thrust_mass_rate_model) when making partial" );
        }
        else
        {
            massRatePartial = std::make_shared< orbit_determination::FromThrustMassRatePartial >
                    ( bodyWithMassRate.first, std::dynamic_pointer_cast< FromThrustMassRateModel >( massRateModel ) );
        }
        break;
    }
    case custom_mass_rate_model:
        throw std::runtime_error( "Custom mass rate model does not yet have an associated partial." );
        break;
    default:
        std::string errorMessage = "Mass rate model " + std::to_string( massRateType ) +
                " not found when making torque partial";
        throw std::runtime_error( errorMessage );
        break;
    }

    return massRatePartial;
}


template< typename InitialStateParameterType >
orbit_determination::StateDerivativePartialsMap createMassRatePartialsMap(
        const basic_astrodynamics::MassRateModelMap& massRateMap,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >
        parametersToEstimate )
{
    // Declare return map.
    orbit_determination::StateDerivativePartialsMap massRatePartialsList;

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatable_parameters::getListOfMassStateParametersToEstimate( parametersToEstimate );
    massRatePartialsList.resize( initialDynamicalParameters.size( ) );

    // Iterate over list of bodies of which the partials of the mass rates acting on them are required.
    for( basic_astrodynamics::MassRateModelMap::const_iterator massRateIterator = massRateMap.begin( );
         massRateIterator != massRateMap.end( ); massRateIterator++ )
    {
        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if( initialDynamicalParameters.at( i )->getParameterName( ).second.first == massRateIterator->first )
            {
                if( ( initialDynamicalParameters.at( i )->getParameterName( ).first ==
                      estimatable_parameters::initial_mass_state ) )
                {
                    // Get object for body undergoing massRate
                    const std::string acceleratedBody = massRateIterator->first;
                    std::shared_ptr< simulation_setup::Body > acceleratedBodyObject = bodies.at( acceleratedBody );

                    // Retrieve list of massRates acting on current body.
                    std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > >  massRateVector =
                            massRateMap.at( acceleratedBody );

                    std::vector< std::shared_ptr< orbit_determination::StateDerivativePartial > > massRatePartialVector;
                    for( unsigned int j = 0; j < massRateVector.size( ); j++ )
                    {
                        // Create single partial object
                        std::shared_ptr< orbit_determination::MassRatePartial > currentMassRatePartial =
                                createAnalyticalMassRatePartial(
                                    massRateVector.at( j ),
                                    std::make_pair( acceleratedBody, acceleratedBodyObject ),
                                    bodies, parametersToEstimate );

                        massRatePartialVector.push_back( currentMassRatePartial );
                    }

                    // Add partials of current body's massRates to vector.
                    massRatePartialsList[ i ] = massRatePartialVector;
                }
            }
        }
    }
    return massRatePartialsList;
}

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
            if( stateDerivativeIterator->second.size( ) > 1 )
            {
                throw std::runtime_error(
                            "Error, cannot yet process multiple separate same type propagators when making partial derivatives of mass state." );
            }
            else
            {
                // Retrieve acceleration models and create partials
                basic_astrodynamics::MassRateModelMap massModelList =
                        std::dynamic_pointer_cast< propagators::BodyMassStateDerivative< StateScalarType, TimeType > >(
                            stateDerivativeIterator->second.at( 0 ) )->getMassRateModels( );
                stateDerivativePartials[ propagators::body_mass_state ] =
                        createMassRatePartialsMap< StateScalarType >(
                            massModelList, bodies, parametersToEstimate );
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
//
//extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< double, double >(
//        const std::unordered_map< propagators::IntegratedStateType,
//        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, double > > > >
//        stateDerivativeModels,
//        const simulation_setup::SystemOfBodies& bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
//        parametersToEstimate );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATESTATEDERIVATIVEPARTIALS_H
