/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/propagation_setup/createMassRateModels.h"

namespace tudat
{

namespace simulation_setup
{


//! Function to create a mass rate model
std::shared_ptr< basic_astrodynamics::MassRateModel >
createMassRateModel(
        const std::string& bodyWithMassRate,
        const std::shared_ptr< MassRateModelSettings > massRateModelSettings,
        const SystemOfBodies& bodies,
        const basic_astrodynamics::AccelerationMap& accelerationModels )
{
    std::shared_ptr< basic_astrodynamics::MassRateModel > massRateModel;

    // Check type of mass rate model.
    switch( massRateModelSettings->massRateType_ )
    {
    case basic_astrodynamics::custom_mass_rate_model:
    {
        // Check input consistency
        std::shared_ptr< CustomMassRateSettings > customMassRateModelSettings =
                std::dynamic_pointer_cast< CustomMassRateSettings >(massRateModelSettings );
        if( customMassRateModelSettings == nullptr )
        {
            throw std::runtime_error( "Error when making cusom mass rate model, input is inconsistent" );
        }
        else
        {
            massRateModel = std::make_shared< basic_astrodynamics::CustomMassRateModel >(
                        customMassRateModelSettings->massRateFunction_ );
        }
        break;
    }
    case basic_astrodynamics::from_thrust_mass_rate_model:
    {
        // Check input consistency
        std::shared_ptr< FromThrustMassRateSettings > fromThrustMassModelSettings =
                std::dynamic_pointer_cast< FromThrustMassRateSettings >(massRateModelSettings );
        if( fromThrustMassModelSettings == nullptr )
        {
            throw std::runtime_error( "Error when making from-engine mass rate model, input is inconsistent" );
        }
        else
        {
            std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > >
                    thrustAccelerations;
            if( accelerationModels.count( bodyWithMassRate ) != 0 )
            {
                if( accelerationModels.at( bodyWithMassRate ).count( bodyWithMassRate ) != 0 )
                {
                    thrustAccelerations = basic_astrodynamics::getAccelerationModelsOfType(
                                accelerationModels.at( bodyWithMassRate ).at( bodyWithMassRate ),
                                basic_astrodynamics::thrust_acceleration );
                }
            }

            if( thrustAccelerations.size( ) == 0 )
            {
                throw std::runtime_error(
                            "Error when making from-thrust mass-rate model for body " + bodyWithMassRate + ", no thrust model is found for this body" );
            }

            std::vector< std::shared_ptr< propulsion::ThrustAcceleration > >
                    explicitThrustAccelerations;

            if( fromThrustMassModelSettings->useAllThrustModels_ == 0 )
            {
                // Retrieve thrust models with the correct id (should be 1)
                for( unsigned int i = 0; i < thrustAccelerations.size( ); i++ )
                {
                    if( std::dynamic_pointer_cast< propulsion::ThrustAcceleration >(
                                thrustAccelerations.at( i ) )->getAssociatedThrustSource( ) ==
                            fromThrustMassModelSettings->associatedThrustSource_ )
                    {
                        explicitThrustAccelerations.push_back( std::dynamic_pointer_cast< propulsion::ThrustAcceleration >(
                                                                   thrustAccelerations.at( i ) ) );
                    }
                }

                if( explicitThrustAccelerations.size( ) != 1 )
                {
                    std::cerr << "Warning when making from-thrust mass-rate model, did not find exactly 1 thrust model with correct identifier" << std::endl;
                }
            }
            else
            {
                // Combine all thrust models into list
                for( unsigned int i = 0; i < thrustAccelerations.size( ); i++ )
                {
                    explicitThrustAccelerations.push_back( std::dynamic_pointer_cast< propulsion::ThrustAcceleration >(
                                                               thrustAccelerations.at( i ) ) );
                }
            }

            // Create mass rate model
            massRateModel = std::make_shared< propulsion::FromThrustMassRateModel >(
                        explicitThrustAccelerations );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error when making mass rate model, type not recognized" );

    }
    return massRateModel;
}


//! Function to create a list of mass rate models for a list of bodies.
basic_astrodynamics::MassRateModelMap createMassRateModelsMap(
        const SystemOfBodies& bodies,
        const SelectedMassRateModelMap& massRateModelSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModels )
{
    // Iterate over all bodies
    std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > > > massRateModels;
    for( std::map< std::string, std::vector< std::shared_ptr< MassRateModelSettings > > >::const_iterator settingsIterator =
         massRateModelSettings.begin( ); settingsIterator != massRateModelSettings.end( ); settingsIterator++)
    {
        // Iterate over all mass model settings for current body.
        for( unsigned int i = 0; i < settingsIterator->second.size( ); i++ )
        {
            massRateModels[ settingsIterator->first ].push_back(
                        createMassRateModel( settingsIterator->first, settingsIterator->second.at( i ), bodies,
                                             accelerationModels ) );
        }
    }
    return massRateModels;

}

} // namespace simulation_setup

} // namespace tudat

