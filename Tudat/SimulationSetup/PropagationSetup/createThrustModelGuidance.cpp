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
boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > createThrustGuidanceModel(
        const boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        const boost::function< Eigen::Vector3d( ) > bodyFixedThrustOrientation,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
{
    boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustGuidance;

    // Determine thrust direction type
    switch( thrustDirectionGuidanceSettings->thrustDirectionType_ )
    {
    case colinear_with_state_segment_thrust_direction:
    {
        // Check input consistency
        boost::shared_ptr< ThrustDirectionFromStateGuidanceSettings > thrustDirectionFromStateGuidanceSettings =
                boost::dynamic_pointer_cast< ThrustDirectionFromStateGuidanceSettings >( thrustDirectionGuidanceSettings );
        if( thrustDirectionFromStateGuidanceSettings == NULL )
        {
            throw std::runtime_error( "Error when thrust colinear with state, input is inconsistent" );
        }
        else
        {
            // Retrieve state function of body for which thrust is to be computed.
            boost::function< Eigen::Vector6d( ) > bodyStateFunction =
                    boost::bind( &Body::getState, bodyMap.at( nameOfBodyWithGuidance ) );
            boost::function< Eigen::Vector6d( ) > centralBodyStateFunction;

            // Retrieve state function of central body (or set to zero if inertial)
            if( thrustDirectionFromStateGuidanceSettings->relativeBody_ != "SSB" &&
                    thrustDirectionFromStateGuidanceSettings->relativeBody_ != "" )
            {
                centralBodyStateFunction = boost::bind( &Body::getState, bodyMap.at(
                                                            thrustDirectionFromStateGuidanceSettings->relativeBody_ ) );
                magnitudeUpdateSettings[ propagators::body_transational_state_update ].push_back(
                            thrustDirectionFromStateGuidanceSettings->relativeBody_ );
            }
            else
            {
                centralBodyStateFunction = boost::lambda::constant( Eigen::Vector6d::Zero( ) );
            }

            // Define relative state function
            boost::function< void( Eigen::Vector6d& ) > stateFunction =
                    boost::bind(
                        &ephemerides::getRelativeState, _1, bodyStateFunction, centralBodyStateFunction );
            boost::function< Eigen::Vector3d( const double ) > thrustDirectionFunction;

            // Create force direction function.
            if( thrustDirectionFromStateGuidanceSettings->isColinearWithVelocity_ )
            {
                thrustDirectionFunction =
                        boost::bind( &propulsion::getForceDirectionColinearWithVelocity, stateFunction, _1,
                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
            }
            else
            {
                thrustDirectionFunction =
                        boost::bind( &propulsion::getForceDirectionColinearWithPosition, stateFunction, _1,
                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
            }

            // Create direction guidance
            thrustGuidance =  boost::make_shared< propulsion::DirectionBasedForceGuidance >(
                        thrustDirectionFunction, thrustDirectionFromStateGuidanceSettings->relativeBody_,
                        bodyFixedThrustOrientation );
        }
        break;
    }
    case thrust_direction_from_existing_body_orientation:
    {
        boost::shared_ptr< Body > bodyWithGuidance = bodyMap.at( nameOfBodyWithGuidance );

        boost::function< Eigen::Quaterniond( const double ) > rotationFunction;

        // Retrieve existing body rotation model and set associated update settings.
        if( bodyWithGuidance->getFlightConditions( ) != NULL )
        {
            rotationFunction = boost::bind(
                        &simulation_setup::Body::getCurrentRotationToGlobalFrame,
                        bodyWithGuidance );

            magnitudeUpdateSettings[ propagators::vehicle_flight_conditions_update ].push_back( nameOfBodyWithGuidance );
            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back(
                        thrustDirectionGuidanceSettings->relativeBody_ );
            magnitudeUpdateSettings[ propagators::body_transational_state_update ].push_back(
                        thrustDirectionGuidanceSettings->relativeBody_);
        }
        else if( bodyWithGuidance->getRotationalEphemeris( ) != NULL )
        {
            rotationFunction = boost::bind(
                        &simulation_setup::Body::getCurrentRotationToGlobalFrame,
                        bodyWithGuidance );
            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
        }
        else
        {
            throw std::runtime_error( "Error, requested thrust orientation from existing model, but no such model found" );
        }

        thrustGuidance =  boost::make_shared< propulsion::OrientationBasedForceGuidance >(
                    rotationFunction, bodyFixedThrustOrientation );

        break;
    }
    case custom_thrust_direction:
    {
        // Check input consistency
        boost::shared_ptr< CustomThrustDirectionSettings > customThrustGuidanceSettings =
                boost::dynamic_pointer_cast< CustomThrustDirectionSettings >( thrustDirectionGuidanceSettings );
        if( customThrustGuidanceSettings == NULL )
        {
            throw std::runtime_error( "Error when getting thrust guidance with custom_thrust_direction, input is inconsistent" );
        }
        else
        {
            // Create direction guidance
            boost::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
                    customThrustGuidanceSettings->thrustDirectionFunction_;
            thrustGuidance =  boost::make_shared< propulsion::DirectionBasedForceGuidance >(
                        thrustDirectionFunction, "", bodyFixedThrustOrientation );
        }

        break;
    }
    case custom_thrust_orientation:
    {
        // Check input consistency
        boost::shared_ptr< CustomThrustOrientationSettings > customThrustOrientationSettings =
                boost::dynamic_pointer_cast< CustomThrustOrientationSettings >( thrustDirectionGuidanceSettings );
        if( customThrustOrientationSettings == NULL )
        {
            throw std::runtime_error( "Error when getting thrust guidance with custom_thrust_orientation, input is inconsistent" );
        }
        else
        {
            // Create direction guidance
            thrustGuidance =  boost::make_shared< propulsion::OrientationBasedForceGuidance >(
                        customThrustOrientationSettings->thrustOrientationFunction_,
                        bodyFixedThrustOrientation );
            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
        }
        break;
    }
    case mee_costate_based_thrust_direction:
    {
        // Check input consistency
        boost::shared_ptr< MeeCostateBasedThrustDirectionSettings > meeCostateBasedThrustSettings =
                boost::dynamic_pointer_cast< MeeCostateBasedThrustDirectionSettings >( thrustDirectionGuidanceSettings );

        if( meeCostateBasedThrustSettings == NULL )
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
            else if( bodyMap.at( meeCostateBasedThrustSettings->relativeBody_ )->getGravityFieldModel( ) == NULL )
            {
                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, central body " +
                                          meeCostateBasedThrustSettings->relativeBody_ + " has no gravity field." );
            }
            else
            {
                // Retrieve required functions and create guidance object
                boost::function< Eigen::Vector6d( ) > thrustingBodyStateFunction =
                        boost::bind( &simulation_setup::Body::getState,
                                     bodyMap.at( meeCostateBasedThrustSettings->vehicleName_ ) );
                boost::function< Eigen::Vector6d( ) > centralBodyStateFunction =
                        boost::bind( &simulation_setup::Body::getState,
                                     bodyMap.at( meeCostateBasedThrustSettings->relativeBody_ ) );
                boost::function< double( ) > centralBodyGravitationalParameterFunction =
                        boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                     bodyMap.at( meeCostateBasedThrustSettings->relativeBody_ )->getGravityFieldModel( ) );

                thrustGuidance =  boost::make_shared< propulsion::MeeCostateBasedThrustGuidance >(
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
        const std::vector< boost::function< Eigen::Vector3d( )> >& thrustDirections,
        const std::vector< boost::function< double( )> >& thrustMagnitudes )
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
boost::function< Eigen::Vector3d( ) > getBodyFixedThrustDirection(
        const boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string bodyName )
{
    boost::function< Eigen::Vector3d( ) > thrustDirectionFunction;

    // Identify magnitude settings type
    switch( thrustMagnitudeSettings->thrustMagnitudeGuidanceType_ )
    {
    case constant_thrust_magnitude:
    {
        // Check input consistency
        boost::shared_ptr< ConstantThrustEngineSettings > constantThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< ConstantThrustEngineSettings >( thrustMagnitudeSettings );
        if( constantThrustMagnitudeSettings == NULL )
        {
            throw std::runtime_error( "Error when creating body-fixed thrust direction of type constant_thrust_magnitude, input is inconsistent" );
        }
        else
        {
            thrustDirectionFunction = boost::lambda::constant( constantThrustMagnitudeSettings->bodyFixedThrustDirection_ );
        }
        break;
    }
    case from_engine_properties_thrust_magnitude:
    {
        // Check input consistency
        boost::shared_ptr< FromBodyThrustEngineSettings > fromEngineThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< FromBodyThrustEngineSettings >( thrustMagnitudeSettings );
        if( fromEngineThrustMagnitudeSettings == NULL )
        {
            throw std::runtime_error( "Error when creating body-fixed thrust direction of type from_engine_properties_thrust_magnitude, input is inconsistent" );
        }
        if( bodyMap.at( bodyName )->getVehicleSystems( ) == NULL )
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
                        boost::bind( &system_models::EngineModel::getBodyFixedThrustDirection,
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
            std::vector< boost::function< Eigen::Vector3d( )> > thrustDirections;
            std::vector< boost::function< double( )> > thrustMagnitudes;

            std::map< std::string, boost::shared_ptr< system_models::EngineModel > > engineModels =
                    bodyMap.at( bodyName )->getVehicleSystems( )->getEngineModels( );

            for( std::map< std::string, boost::shared_ptr< system_models::EngineModel > >::const_iterator engineIterator =
                 engineModels.begin( ); engineIterator != engineModels.end( ); engineIterator++ )
            {
                thrustDirections.push_back(
                            boost::bind( &system_models::EngineModel::getBodyFixedThrustDirection, engineIterator->second ) );
                thrustMagnitudes.push_back(
                            boost::bind( &system_models::EngineModel::getCurrentThrust, engineIterator->second ) );
            }

            // Create effective thrust direction function.
            thrustDirectionFunction = boost::bind(
                        &getCombinedThrustDirection, thrustDirections, thrustMagnitudes );

        }
        break;

    }
    case thrust_magnitude_from_time_function:
    {
        // Check input consistency
        boost::shared_ptr< FromFunctionThrustEngineSettings > fromFunctionThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< FromFunctionThrustEngineSettings >( thrustMagnitudeSettings );
        if( fromFunctionThrustMagnitudeSettings == NULL )
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
        boost::shared_ptr< ParameterizedThrustMagnitudeSettings > fromFunctionThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< ParameterizedThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( fromFunctionThrustMagnitudeSettings == NULL )
        {
            throw std::runtime_error( "Error when creating body-fixed thrust direction of type thrust_magnitude_from_dependent_variables, input is inconsistent" );
        }
        else
        {
            thrustDirectionFunction = boost::lambda::constant( fromFunctionThrustMagnitudeSettings->bodyFixedThrustDirection_ );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error when creating body-fixed thrust direction, type not identified" );
    }
    return thrustDirectionFunction;
}

//! Function to create a wrapper object that computes the thrust magnitude
boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > createThrustMagnitudeWrapper(
        const boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
{
    boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper;

    // Identify magnitude settings type
    switch( thrustMagnitudeSettings->thrustMagnitudeGuidanceType_ )
    {
    case constant_thrust_magnitude:
    {
        // Check input consistency
        boost::shared_ptr< ConstantThrustEngineSettings > constantThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< ConstantThrustEngineSettings >( thrustMagnitudeSettings );
        if( constantThrustMagnitudeSettings == NULL )
        {
            throw std::runtime_error( "Error when creating constant thrust magnitude wrapper, input is inconsistent" );
        }

        thrustMagnitudeWrapper = boost::make_shared< propulsion::CustomThrustMagnitudeWrapper >(
                    boost::lambda::constant( constantThrustMagnitudeSettings->thrustMagnitude_ ),
                    boost::lambda::constant( constantThrustMagnitudeSettings->specificImpulse_ ) );
        break;

    }
    case from_engine_properties_thrust_magnitude:
    {
        // Check input consistency
        boost::shared_ptr< FromBodyThrustEngineSettings > fromEngineThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< FromBodyThrustEngineSettings >( thrustMagnitudeSettings );
        if( fromEngineThrustMagnitudeSettings == NULL )
        {
            throw std::runtime_error( "Error when creating from-engine thrust magnitude wrapper, input is inconsistent" );
        }
        if( bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( ) == NULL )
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
                thrustMagnitudeWrapper = boost::make_shared< propulsion::ThrustMagnitudeFromEngineWrapper >(
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
            thrustMagnitudeWrapper = boost::make_shared< propulsion::ThrustMagnitudeFromEngineWrapper >(
                        utilities::createVectorFromMapValues< boost::shared_ptr< system_models::EngineModel >, std::string >(
                            bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( )->getEngineModels( ) ));
        }
        break;

    }
    case thrust_magnitude_from_time_function:
    {
        // Check input consistency
        boost::shared_ptr< FromFunctionThrustEngineSettings > fromFunctionThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< FromFunctionThrustEngineSettings >( thrustMagnitudeSettings );
        if( fromFunctionThrustMagnitudeSettings == NULL )
        {
            throw std::runtime_error( "Error when creating from-function thrust magnitude wrapper, input is inconsistent" );
        }
        thrustMagnitudeWrapper = boost::make_shared< propulsion::CustomThrustMagnitudeWrapper >(
                    fromFunctionThrustMagnitudeSettings->thrustMagnitudeFunction_,
                    fromFunctionThrustMagnitudeSettings->specificImpulseFunction_,
                    fromFunctionThrustMagnitudeSettings->isEngineOnFunction_,
                    fromFunctionThrustMagnitudeSettings->customThrustResetFunction_ );
        break;

    }
    case thrust_magnitude_from_dependent_variables:
    {
        // Check input consistency
        boost::shared_ptr< ParameterizedThrustMagnitudeSettings > parameterizedThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< ParameterizedThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( parameterizedThrustMagnitudeSettings == NULL )
        {
            throw std::runtime_error( "Error when creating from-function thrust magnitude wrapper, input is inconsistent" );
        }

        // Create indpendent variable functions
        std::vector< boost::function< double( ) > > thrustInputVariableFunctions =
                getPropulsionInputVariables(
                    bodyMap.at( nameOfBodyWithGuidance ), parameterizedThrustMagnitudeSettings->thrustIndependentVariables_,
                    parameterizedThrustMagnitudeSettings->thrustGuidanceInputVariables_ );
        std::vector< boost::function< double( ) > > specificInputVariableFunctions =
                getPropulsionInputVariables(
                    bodyMap.at( nameOfBodyWithGuidance ), parameterizedThrustMagnitudeSettings->specificImpulseDependentVariables_,
                    parameterizedThrustMagnitudeSettings->specificImpulseGuidanceInputVariables_ );

        // Create thrust magnitude wrapper
        thrustMagnitudeWrapper = boost::make_shared< propulsion::ParameterizedThrustMagnitudeWrapper >(
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
        const boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
        const double currentTime )
{
    thrustMagnitudeWrapper->update( currentTime );
    thrustDirectionGuidance->updateCalculator( currentTime );
}

//! Function to reset the current time variable of the thrust magnitude and direction wrappers
void resetThrustMagnitudeAndDirectionTime(
        const boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
        const double currentTime )
{
    thrustMagnitudeWrapper->resetCurrentTime( currentTime );
    thrustDirectionGuidance->resetCurrentTime( currentTime );
}


}

}

