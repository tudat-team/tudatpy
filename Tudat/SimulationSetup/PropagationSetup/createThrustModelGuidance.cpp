/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Propulsion/costateBasedThrustGuidance.h"
#include "Tudat/SimulationSetup/PropagationSetup/createThrustModelGuidance.h"
#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create the object determining the direction of the thrust acceleration.
std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > createThrustGuidanceModel(
        const std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        const std::function< Eigen::Vector3d( ) > bodyFixedThrustOrientation,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
{
    std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustGuidance;

    // Determine thrust direction type
    switch( thrustDirectionGuidanceSettings->thrustDirectionType_ )
    {
    case colinear_with_state_segment_thrust_direction:
    {
        // Check input consistency
        std::shared_ptr< ThrustDirectionFromStateGuidanceSettings > thrustDirectionFromStateGuidanceSettings =
                std::dynamic_pointer_cast< ThrustDirectionFromStateGuidanceSettings >( thrustDirectionGuidanceSettings );
        if( thrustDirectionFromStateGuidanceSettings == nullptr )
        {
            throw std::runtime_error( "Error when thrust colinear with state, input is inconsistent" );
        }
        else
        {
            // Retrieve state function of body for which thrust is to be computed.
            std::function< Eigen::Vector6d( ) > bodyStateFunction =
                    std::bind( &Body::getState, bodyMap.at( nameOfBodyWithGuidance ) );
            std::function< Eigen::Vector6d( ) > centralBodyStateFunction;

            // Retrieve state function of central body (or set to zero if inertial)
            if( thrustDirectionFromStateGuidanceSettings->relativeBody_ != "SSB" &&
                    thrustDirectionFromStateGuidanceSettings->relativeBody_ != "" )
            {
                centralBodyStateFunction = std::bind( &Body::getState, bodyMap.at(
                                                            thrustDirectionFromStateGuidanceSettings->relativeBody_ ) );
                magnitudeUpdateSettings[ propagators::body_translational_state_update ].push_back(
                            thrustDirectionFromStateGuidanceSettings->relativeBody_ );
            }
            else
            {
                centralBodyStateFunction = [ ]( ){ return Eigen::Vector6d::Zero( ); };
            }

            // Define relative state function
            std::function< void( Eigen::Vector6d& ) > stateFunction =
                    std::bind(
                        &ephemerides::getRelativeState, std::placeholders::_1, bodyStateFunction, centralBodyStateFunction );
            std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction;

            // Create force direction function.
            if( thrustDirectionFromStateGuidanceSettings->isColinearWithVelocity_ )
            {
                thrustDirectionFunction =
                        std::bind( &propulsion::getForceDirectionColinearWithVelocity, stateFunction, std::placeholders::_1,
                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
            }
            else
            {
                thrustDirectionFunction =
                        std::bind( &propulsion::getForceDirectionColinearWithPosition, stateFunction, std::placeholders::_1,
                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
            }

            // Create direction guidance
            thrustGuidance =  std::make_shared< propulsion::DirectionBasedForceGuidance >(
                        thrustDirectionFunction, thrustDirectionFromStateGuidanceSettings->relativeBody_,
                        bodyFixedThrustOrientation );
        }
        break;
    }
    case thrust_direction_from_existing_body_orientation:
    {
        std::shared_ptr< Body > bodyWithGuidance = bodyMap.at( nameOfBodyWithGuidance );

        std::function< Eigen::Quaterniond( const double ) > rotationFunction;

        // Retrieve existing body rotation model and set associated update settings.
        if( bodyWithGuidance->getFlightConditions( ) != nullptr )
        {
            rotationFunction = std::bind(
                        &simulation_setup::Body::getCurrentRotationToGlobalFrame,
                        bodyWithGuidance );

            magnitudeUpdateSettings[ propagators::vehicle_flight_conditions_update ].push_back( nameOfBodyWithGuidance );
            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back(
                        thrustDirectionGuidanceSettings->relativeBody_ );
            magnitudeUpdateSettings[ propagators::body_translational_state_update ].push_back(
                        thrustDirectionGuidanceSettings->relativeBody_);
        }
        else if( bodyWithGuidance->getRotationalEphemeris( ) != nullptr )
        {
            rotationFunction = std::bind(
                        &simulation_setup::Body::getCurrentRotationToGlobalFrame,
                        bodyWithGuidance );
            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
        }
        else
        {
            throw std::runtime_error( "Error, requested thrust orientation from existing model, but no such model found" );
        }

        thrustGuidance =  std::make_shared< propulsion::OrientationBasedForceGuidance >(
                    rotationFunction, bodyFixedThrustOrientation );

        break;
    }
    case custom_thrust_direction:
    {
        // Check input consistency
        std::shared_ptr< CustomThrustDirectionSettings > customThrustGuidanceSettings =
                std::dynamic_pointer_cast< CustomThrustDirectionSettings >( thrustDirectionGuidanceSettings );
        if( customThrustGuidanceSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting thrust guidance with custom_thrust_direction, input is inconsistent" );
        }
        else
        {
            // Create direction guidance
            std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
                    customThrustGuidanceSettings->thrustDirectionFunction_;
            thrustGuidance =  std::make_shared< propulsion::DirectionBasedForceGuidance >(
                        thrustDirectionFunction, "", bodyFixedThrustOrientation );
        }

        break;
    }
    case custom_thrust_orientation:
    {
        // Check input consistency
        std::shared_ptr< CustomThrustOrientationSettings > customThrustOrientationSettings =
                std::dynamic_pointer_cast< CustomThrustOrientationSettings >( thrustDirectionGuidanceSettings );
        if( customThrustOrientationSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting thrust guidance with custom_thrust_orientation, input is inconsistent" );
        }
        else
        {
            // Create direction guidance
            thrustGuidance =  std::make_shared< propulsion::OrientationBasedForceGuidance >(
                        customThrustOrientationSettings->thrustOrientationFunction_,
                        bodyFixedThrustOrientation );
            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
        }
        break;
    }
    case mee_costate_based_thrust_direction:
    {
        // Check input consistency
        std::shared_ptr< MeeCostateBasedThrustDirectionSettings > meeCostateBasedThrustSettings =
                std::dynamic_pointer_cast< MeeCostateBasedThrustDirectionSettings >( thrustDirectionGuidanceSettings );

        if( meeCostateBasedThrustSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, input is inconsistent" );
        }
        else
        {
            // Check whether all required environment properties exist
            if( bodyMap.count( meeCostateBasedThrustSettings->relativeBody_ ) == 0 )
            {
                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, central body " +
                                          meeCostateBasedThrustSettings->relativeBody_ + " not found." );
            }
            else if( bodyMap.count( meeCostateBasedThrustSettings->vehicleName_ ) == 0 )
            {
                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, thrusting body " +
                                          meeCostateBasedThrustSettings->vehicleName_ + " not found." );
            }
            else if( bodyMap.at( meeCostateBasedThrustSettings->relativeBody_ )->getGravityFieldModel( ) == nullptr )
            {
                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, central body " +
                                          meeCostateBasedThrustSettings->relativeBody_ + " has no gravity field." );
            }
            else
            {
                // Retrieve required functions and create guidance object
                std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction =
                        std::bind( &simulation_setup::Body::getState,
                                     bodyMap.at( meeCostateBasedThrustSettings->vehicleName_ ) );
                std::function< Eigen::Vector6d( ) > centralBodyStateFunction =
                        std::bind( &simulation_setup::Body::getState,
                                     bodyMap.at( meeCostateBasedThrustSettings->relativeBody_ ) );
                std::function< double( ) > centralBodyGravitationalParameterFunction =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                     bodyMap.at( meeCostateBasedThrustSettings->relativeBody_ )->getGravityFieldModel( ) );

                thrustGuidance =  std::make_shared< propulsion::MeeCostateBasedThrustGuidance >(
                            thrustingBodyStateFunction, centralBodyStateFunction,
                            centralBodyGravitationalParameterFunction,
                            meeCostateBasedThrustSettings->costateFunction_,
                            bodyFixedThrustOrientation );
            }
        }
        break;
    }
    default:
        throw std::runtime_error( "Error, could not find thrust guidance type when creating thrust guidance." );
    }
    return thrustGuidance;
}

//! Function to retrieve the effective thrust direction from a set of thrust sources.
Eigen::Vector3d getCombinedThrustDirection(
        const std::vector< std::function< Eigen::Vector3d( )> >& thrustDirections,
        const std::vector< std::function< double( )> >& thrustMagnitudes )
{
    Eigen::Vector3d thrustDirection = Eigen::Vector3d::Zero( );
    double totalThrust = 0.0;

    for( unsigned int i = 0; i < thrustDirections.size( ); i++ )
    {
        thrustDirection += thrustMagnitudes.at( i )( ) * thrustDirections.at( i )( );
        totalThrust += thrustMagnitudes.at( i )( );
    }
    return thrustDirection / totalThrust;
}

//! Function to create a function that returns the thrust direction in the body-fixed frame.
std::function< Eigen::Vector3d( ) > getBodyFixedThrustDirection(
        const std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string bodyName )
{
    std::function< Eigen::Vector3d( ) > thrustDirectionFunction;

    // Identify magnitude settings type
    switch( thrustMagnitudeSettings->thrustMagnitudeGuidanceType_ )
    {
    case constant_thrust_magnitude:
    {
        // Check input consistency
        std::shared_ptr< ConstantThrustMagnitudeSettings > constantThrustMagnitudeSettings =
                std::dynamic_pointer_cast< ConstantThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( constantThrustMagnitudeSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating body-fixed thrust direction of type constant_thrust_magnitude, input is inconsistent" );
        }
        else
        {
            thrustDirectionFunction = [ = ]( ){ return constantThrustMagnitudeSettings->bodyFixedThrustDirection_; };
        }
        break;
    }
    case from_engine_properties_thrust_magnitude:
    {
        // Check input consistency
        std::shared_ptr< FromBodyThrustMagnitudeSettings > fromEngineThrustMagnitudeSettings =
                std::dynamic_pointer_cast< FromBodyThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( fromEngineThrustMagnitudeSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating body-fixed thrust direction of type from_engine_properties_thrust_magnitude, input is inconsistent" );
        }
        if( bodyMap.at( bodyName )->getVehicleSystems( ) == nullptr )
        {
            throw std::runtime_error( "Error when creating body-fixed thrust direction of type from_engine_properties_thrust_magnitude, no vehicle systems found" );

        }

        // Retrieve single engine
        if( fromEngineThrustMagnitudeSettings->useAllEngines_ == false  )
        {
            // Check if engine model exists
            if( ( bodyMap.at( bodyName )->getVehicleSystems( )->getEngineModels( ).count(
                      thrustMagnitudeSettings->thrustOriginId_ ) == 0 ) )
            {
                throw std::runtime_error( "Error when creating body-fixed thrust direction of type from_engine_properties_thrust_magnitude, no engine of right ID found" );
            }
            else
            {
                thrustDirectionFunction =
                        std::bind( &system_models::EngineModel::getBodyFixedThrustDirection,
                                     bodyMap.at( bodyName )->getVehicleSystems( )->getEngineModels( ).at(
                                         thrustMagnitudeSettings->thrustOriginId_ ) );
            }

        }
        // Retrieve mean thrust direction from all engines
        else
        {
            // Print warning if there are no engines (zero thrust)
            if( ( bodyMap.at( bodyName )->getVehicleSystems( )->getEngineModels( ).size( ) == 0 ) )
            {
                std::cerr << "Warning when creating body-fixed thrust direction of type from_engine_properties_thrust_magnitude; no engines found: returning 0 thrust" << std::endl;
            }

            // Retrieve force directions/magnitudes
            std::vector< std::function< Eigen::Vector3d( )> > thrustDirections;
            std::vector< std::function< double( )> > thrustMagnitudes;

            std::map< std::string, std::shared_ptr< system_models::EngineModel > > engineModels =
                    bodyMap.at( bodyName )->getVehicleSystems( )->getEngineModels( );

            for( std::map< std::string, std::shared_ptr< system_models::EngineModel > >::const_iterator engineIterator =
                 engineModels.begin( ); engineIterator != engineModels.end( ); engineIterator++ )
            {
                thrustDirections.push_back(
                            std::bind( &system_models::EngineModel::getBodyFixedThrustDirection, engineIterator->second ) );
                thrustMagnitudes.push_back(
                            std::bind( &system_models::EngineModel::getCurrentThrust, engineIterator->second ) );
            }

            // Create effective thrust direction function.
            thrustDirectionFunction = std::bind(
                        &getCombinedThrustDirection, thrustDirections, thrustMagnitudes );

        }
        break;

    }
    case thrust_magnitude_from_time_function:
    {
        // Check input consistency
        std::shared_ptr< FromFunctionThrustMagnitudeSettings > fromFunctionThrustMagnitudeSettings =
                std::dynamic_pointer_cast< FromFunctionThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( fromFunctionThrustMagnitudeSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating body-fixed thrust direction of type thrust_magnitude_from_time_function, input is inconsistent" );
        }
        else
        {
            thrustDirectionFunction =  fromFunctionThrustMagnitudeSettings->bodyFixedThrustDirection_;
        }
        break;

    }
    case thrust_magnitude_from_dependent_variables:
    {
        // Check input consistency
        std::shared_ptr< ParameterizedThrustMagnitudeSettings > fromFunctionThrustMagnitudeSettings =
                std::dynamic_pointer_cast< ParameterizedThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( fromFunctionThrustMagnitudeSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating body-fixed thrust direction of type thrust_magnitude_from_dependent_variables, input is inconsistent" );
        }
        else
        {
            thrustDirectionFunction = [ = ]( ){ return fromFunctionThrustMagnitudeSettings->bodyFixedThrustDirection_; };
        }
        break;
    }
    default:
        throw std::runtime_error( "Error when creating body-fixed thrust direction, type not identified" );
    }
    return thrustDirectionFunction;
}

//! Function to create a wrapper object that computes the thrust magnitude
std::shared_ptr< propulsion::ThrustMagnitudeWrapper > createThrustMagnitudeWrapper(
        const std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
{
    std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper;

    // Identify magnitude settings type
    switch( thrustMagnitudeSettings->thrustMagnitudeGuidanceType_ )
    {
    case constant_thrust_magnitude:
    {
        // Check input consistency
        std::shared_ptr< ConstantThrustMagnitudeSettings > constantThrustMagnitudeSettings =
                std::dynamic_pointer_cast< ConstantThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( constantThrustMagnitudeSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating constant thrust magnitude wrapper, input is inconsistent" );
        }

        thrustMagnitudeWrapper = std::make_shared< propulsion::CustomThrustMagnitudeWrapper >(
                    [ = ]( const double ){ return constantThrustMagnitudeSettings->thrustMagnitude_; },
                    [ = ]( const double ){ return constantThrustMagnitudeSettings->specificImpulse_; } );
        break;

    }
    case from_engine_properties_thrust_magnitude:
    {
        // Check input consistency
        std::shared_ptr< FromBodyThrustMagnitudeSettings > fromEngineThrustMagnitudeSettings =
                std::dynamic_pointer_cast< FromBodyThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( fromEngineThrustMagnitudeSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating from-engine thrust magnitude wrapper, input is inconsistent" );
        }
        if( bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( ) == nullptr )
        {
            throw std::runtime_error( "Error when creating from-engine thrust magnitude wrapper, no vehicle systems found" );

        }

        // Retrieve single engine thrust
        if( fromEngineThrustMagnitudeSettings->useAllEngines_ == false  )
        {
            if( ( bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( )->getEngineModels( ).count(
                      thrustMagnitudeSettings->thrustOriginId_ ) == 0 ) )
            {
                throw std::runtime_error( "Error when creating from-engine thrust magnitude wrapper, no engine of right ID found" );
            }
            else
            {
                thrustMagnitudeWrapper = std::make_shared< propulsion::ThrustMagnitudeFromEngineWrapper >(
                            bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( )->getEngineModels( ).at(
                                thrustMagnitudeSettings->thrustOriginId_ ) );
            }
        }
        // Retrieve total engine thrust
        else
        {
            if( ( bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( )->getEngineModels( ).size( ) == 0 ) )
            {
                std::cerr << "Warning when creating from-engine thrust magnitude wrapper for all engines; no engines found: returning 0 thrust" << std::endl;
            }
            thrustMagnitudeWrapper = std::make_shared< propulsion::ThrustMagnitudeFromEngineWrapper >(
                        utilities::createVectorFromMapValues< std::shared_ptr< system_models::EngineModel >, std::string >(
                            bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( )->getEngineModels( ) ));
        }
        break;

    }
    case thrust_magnitude_from_time_function:
    {
        // Check input consistency
        std::shared_ptr< FromFunctionThrustMagnitudeSettings > fromFunctionThrustMagnitudeSettings =
                std::dynamic_pointer_cast< FromFunctionThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( fromFunctionThrustMagnitudeSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating from-function thrust magnitude wrapper, input is inconsistent" );
        }
        thrustMagnitudeWrapper = std::make_shared< propulsion::CustomThrustMagnitudeWrapper >(
                    fromFunctionThrustMagnitudeSettings->thrustMagnitudeFunction_,
                    fromFunctionThrustMagnitudeSettings->specificImpulseFunction_,
                    fromFunctionThrustMagnitudeSettings->isEngineOnFunction_,
                    fromFunctionThrustMagnitudeSettings->customThrustResetFunction_ );
        break;

    }
    case thrust_magnitude_from_dependent_variables:
    {
        // Check input consistency
        std::shared_ptr< ParameterizedThrustMagnitudeSettings > parameterizedThrustMagnitudeSettings =
                std::dynamic_pointer_cast< ParameterizedThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( parameterizedThrustMagnitudeSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating from-function thrust magnitude wrapper, input is inconsistent" );
        }

        // Create indpendent variable functions
        std::vector< std::function< double( ) > > thrustInputVariableFunctions =
                getPropulsionInputVariables(
                    bodyMap.at( nameOfBodyWithGuidance ), parameterizedThrustMagnitudeSettings->thrustIndependentVariables_,
                    parameterizedThrustMagnitudeSettings->thrustGuidanceInputVariables_ );
        std::vector< std::function< double( ) > > specificInputVariableFunctions =
                getPropulsionInputVariables(
                    bodyMap.at( nameOfBodyWithGuidance ), parameterizedThrustMagnitudeSettings->specificImpulseDependentVariables_,
                    parameterizedThrustMagnitudeSettings->specificImpulseGuidanceInputVariables_ );

        // Create thrust magnitude wrapper
        thrustMagnitudeWrapper = std::make_shared< propulsion::ParameterizedThrustMagnitudeWrapper >(
                    parameterizedThrustMagnitudeSettings->thrustMagnitudeFunction_,
                    parameterizedThrustMagnitudeSettings->specificImpulseFunction_,
                    thrustInputVariableFunctions,
                    specificInputVariableFunctions,
                    parameterizedThrustMagnitudeSettings->thrustIndependentVariables_,
                    parameterizedThrustMagnitudeSettings->specificImpulseDependentVariables_,
                    parameterizedThrustMagnitudeSettings->inputUpdateFunction_ );

        break;

    }
    default:
        throw std::runtime_error( "Error when creating thrust magnitude wrapper, type not identified" );
    }

    return thrustMagnitudeWrapper;
}

//! Function to update the thrust magnitude and direction to current time.
void updateThrustMagnitudeAndDirection(
        const std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
        const double currentTime )
{
    thrustMagnitudeWrapper->update( currentTime );
    thrustDirectionGuidance->updateCalculator( currentTime );
}

//! Function to reset the current time variable of the thrust magnitude and direction wrappers
void resetThrustMagnitudeAndDirectionTime(
        const std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
        const double currentTime )
{
    thrustMagnitudeWrapper->resetCurrentTime( currentTime );
    thrustDirectionGuidance->resetCurrentTime( currentTime );
}


}

}

