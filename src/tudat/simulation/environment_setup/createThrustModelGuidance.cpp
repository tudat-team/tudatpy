/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createThrustModelGuidance.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to retrieve the effective thrust direction from a set of thrust sources.
Eigen::Vector3d getCombinedThrustDirection( const std::vector< std::function< Eigen::Vector3d( ) > >& thrustDirections,
                                            const std::vector< std::function< double( ) > >& thrustMagnitudes )
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

//! Function to create a wrapper object that computes the thrust magnitude
std::shared_ptr< propulsion::ThrustMagnitudeWrapper > createThrustMagnitudeWrapper(
        const std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings,
        const SystemOfBodies& bodies,
        const std::string& nameOfBodyWithGuidance,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
{
    std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper;

    // Identify magnitude settings type
    switch( thrustMagnitudeSettings->thrustMagnitudeType_ )
    {
        case constant_thrust_magnitude: {
            // Check input consistency
            std::shared_ptr< ConstantThrustMagnitudeSettings > constantThrustMagnitudeSettings =
                    std::dynamic_pointer_cast< ConstantThrustMagnitudeSettings >( thrustMagnitudeSettings );
            if( constantThrustMagnitudeSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating constant thrust magnitude wrapper, input is inconsistent" );
            }

            thrustMagnitudeWrapper = std::make_shared< propulsion::ConstantThrustMagnitudeWrapper >(
                    constantThrustMagnitudeSettings->thrustMagnitude_, constantThrustMagnitudeSettings->specificImpulse_ );
            break;
        }
        case thrust_magnitude_from_time_function: {
            // Check input consistency
            std::shared_ptr< CustomThrustMagnitudeSettings > fromFunctionThrustMagnitudeSettings =
                    std::dynamic_pointer_cast< CustomThrustMagnitudeSettings >( thrustMagnitudeSettings );
            if( fromFunctionThrustMagnitudeSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating from-function thrust magnitude wrapper, input is inconsistent" );
            }
            if( fromFunctionThrustMagnitudeSettings->inputIsForce_ )
            {
                if( !fromFunctionThrustMagnitudeSettings->specificImpulseIsConstant_ )
                {
                    thrustMagnitudeWrapper = std::make_shared< propulsion::CustomThrustMagnitudeWrapper >(
                            fromFunctionThrustMagnitudeSettings->thrustMagnitudeFunction_,
                            fromFunctionThrustMagnitudeSettings->specificImpulseFunction_ );
                }
                else
                {
                    thrustMagnitudeWrapper = std::make_shared< propulsion::CustomThrustMagnitudeWrapper >(
                            fromFunctionThrustMagnitudeSettings->thrustMagnitudeFunction_,
                            fromFunctionThrustMagnitudeSettings->specificImpulseFunction_( TUDAT_NAN ) );
                }
            }
            else
            {
                if( !fromFunctionThrustMagnitudeSettings->specificImpulseIsConstant_ )
                {
                    thrustMagnitudeWrapper = std::make_shared< propulsion::CustomThrustAccelerationMagnitudeWrapper >(
                            fromFunctionThrustMagnitudeSettings->thrustMagnitudeFunction_,
                            fromFunctionThrustMagnitudeSettings->specificImpulseFunction_ );
                }
                else
                {
                    thrustMagnitudeWrapper = std::make_shared< propulsion::CustomThrustAccelerationMagnitudeWrapper >(
                            fromFunctionThrustMagnitudeSettings->thrustMagnitudeFunction_,
                            fromFunctionThrustMagnitudeSettings->specificImpulseFunction_( TUDAT_NAN ) );
                }
            }
            break;
        }
        case thrust_magnitude_from_dependent_variables: {
            // Check input consistency
            std::shared_ptr< ParameterizedThrustMagnitudeSettings > parameterizedThrustMagnitudeSettings =
                    std::dynamic_pointer_cast< ParameterizedThrustMagnitudeSettings >( thrustMagnitudeSettings );
            if( parameterizedThrustMagnitudeSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating from-function thrust magnitude wrapper, input is inconsistent" );
            }

            // Create indpendent variable functions
            std::vector< std::function< double( ) > > thrustInputVariableFunctions =
                    getPropulsionInputVariables( bodies.at( nameOfBodyWithGuidance ),
                                                 parameterizedThrustMagnitudeSettings->thrustIndependentVariables_,
                                                 parameterizedThrustMagnitudeSettings->thrustGuidanceInputVariables_ );
            std::vector< std::function< double( ) > > specificInputVariableFunctions =
                    getPropulsionInputVariables( bodies.at( nameOfBodyWithGuidance ),
                                                 parameterizedThrustMagnitudeSettings->specificImpulseDependentVariables_,
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

////! Function to update the thrust magnitude and direction to current time.
// void updateThrustSettings(
//         const std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
//         const std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
//         const double currentTime )
//{
//     thrustMagnitudeWrapper->update( currentTime );
//     thrustDirectionGuidance->updateCalculator( currentTime );
// }

////! Function to reset the current time variable of the thrust magnitude and direction wrappers
// void resetThrustSettingsTime(
//         const std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
//         const std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
//         const double currentTime )
//{
//     thrustMagnitudeWrapper->resetCurrentTime( currentTime );
//     thrustDirectionGuidance->resetCurrentTime( currentTime );
// }

}  // namespace simulation_setup

}  // namespace tudat
