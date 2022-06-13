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
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace simulation_setup
{

////! Function to create the object determining the direction of the thrust acceleration.
//std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > createThrustGuidanceModel(
//        const std::shared_ptr< ThrustDirectionSettings > thrustDirectionGuidanceSettings,
//        const SystemOfBodies& bodies,
//        const std::string& nameOfBodyWithGuidance,
//        const std::function< Eigen::Vector3d( ) > bodyFixedThrustOrientation,
//        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
//{
//    std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustGuidance;

//    // Determine thrust direction type
//    switch( thrustDirectionGuidanceSettings->thrustDirectionType_ )
//    {
////    case colinear_with_state_segment_thrust_direction:
////    {
////        // Check input consistency
////        std::shared_ptr< ThrustDirectionFromStateGuidanceSettings > thrustDirectionFromStateGuidanceSettings =
////                std::dynamic_pointer_cast< ThrustDirectionFromStateGuidanceSettings >( thrustDirectionGuidanceSettings );
////        if( thrustDirectionFromStateGuidanceSettings == nullptr )
////        {
////            throw std::runtime_error( "Error when thrust colinear with state, input is inconsistent" );
////        }
////        else
////        {
////            // Retrieve state function of body for which thrust is to be computed.
////            std::function< Eigen::Vector6d( ) > bodyStateFunction =
////                    std::bind( &Body::getState, bodies.at( nameOfBodyWithGuidance ) );
////            std::function< Eigen::Vector6d( ) > centralBodyStateFunction;

////            // Retrieve state function of central body (or set to zero if inertial)
////            if( thrustDirectionFromStateGuidanceSettings->relativeBody_ != "SSB" &&
////                    thrustDirectionFromStateGuidanceSettings->relativeBody_ != "" )
////            {
////                centralBodyStateFunction = std::bind( &Body::getState, bodies.at(
////                                                            thrustDirectionFromStateGuidanceSettings->relativeBody_ ) );
////                magnitudeUpdateSettings[ propagators::body_translational_state_update ].push_back(
////                            thrustDirectionFromStateGuidanceSettings->relativeBody_ );
////            }
////            else
////            {
////                centralBodyStateFunction = [ ]( ){ return Eigen::Vector6d::Zero( ); };
////            }

////            // Define relative state function
////            std::function< void( Eigen::Vector6d& ) > stateFunction =
////                    std::bind(
////                        &ephemerides::getRelativeState, std::placeholders::_1, bodyStateFunction, centralBodyStateFunction );
////            std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction;

////            // Create force direction function.
////            if( thrustDirectionFromStateGuidanceSettings->isColinearWithVelocity_ )
////            {
////                thrustDirectionFunction =
////                        std::bind( &propulsion::getForceDirectionColinearWithVelocity, stateFunction, std::placeholders::_1,
////                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
////            }
////            else
////            {
////                thrustDirectionFunction =
////                        std::bind( &propulsion::getForceDirectionColinearWithPosition, stateFunction, std::placeholders::_1,
////                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
////            }

////            // Create direction guidance
////            thrustGuidance =  std::make_shared< propulsion::DirectionBasedForceGuidance >(
////                        thrustDirectionFunction, thrustDirectionFromStateGuidanceSettings->relativeBody_,
////                        bodyFixedThrustOrientation );
////        }
////        break;
////    }
//    case thrust_direction_from_existing_body_orientation:
//    {
//        std::shared_ptr< Body > bodyWithGuidance = bodies.at( nameOfBodyWithGuidance );

//        std::function< Eigen::Quaterniond( const double ) > rotationFunction;

//        // Retrieve existing body rotation model and set associated update settings.
//        if( bodyWithGuidance->getFlightConditions( ) != nullptr )
//        {
//            rotationFunction = std::bind(
//                        &simulation_setup::Body::getCurrentRotationToGlobalFrame,
//                        bodyWithGuidance );

//            magnitudeUpdateSettings[ propagators::vehicle_flight_conditions_update ].push_back( nameOfBodyWithGuidance );
//            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
//        }
//        else if( bodyWithGuidance->getRotationalEphemeris( ) != nullptr )
//        {
//            rotationFunction = std::bind(
//                        &simulation_setup::Body::getCurrentRotationToGlobalFrame,
//                        bodyWithGuidance );
//            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
//        }
//        else
//        {
//            throw std::runtime_error( "Error, requested thrust orientation from existing model, but no such model found" );
//        }

//        thrustGuidance =  std::make_shared< propulsion::OrientationBasedForceGuidance >(
//                    rotationFunction, bodyFixedThrustOrientation );

//        break;
//    }
//    case custom_thrust_direction:
//    {
//        // Check input consistency
//        std::shared_ptr< CustomThrustDirectionSettings > customThrustGuidanceSettings =
//                std::dynamic_pointer_cast< CustomThrustDirectionSettings >( thrustDirectionGuidanceSettings );
//        if( customThrustGuidanceSettings == nullptr )
//        {
//            throw std::runtime_error( "Error when getting thrust guidance with custom_thrust_direction, input is inconsistent" );
//        }
//        else
//        {
//            // Create direction guidance
//            std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
//                    customThrustGuidanceSettings->thrustDirectionFunction_;
//            thrustGuidance =  std::make_shared< propulsion::DirectionBasedForceGuidance >(
//                        thrustDirectionFunction, "", bodyFixedThrustOrientation );
//        }

//        break;
//    }
//    case custom_thrust_orientation:
//    {
//        // Check input consistency
//        std::shared_ptr< CustomThrustOrientationSettings > customThrustOrientationSettings =
//                std::dynamic_pointer_cast< CustomThrustOrientationSettings >( thrustDirectionGuidanceSettings );
//        if( customThrustOrientationSettings == nullptr )
//        {
//            throw std::runtime_error( "Error when getting thrust guidance with custom_thrust_orientation, input is inconsistent" );
//        }
//        else
//        {
//            // Create direction guidance
//            thrustGuidance =  std::make_shared< propulsion::OrientationBasedForceGuidance >(
//                        customThrustOrientationSettings->thrustOrientationFunction_,
//                        bodyFixedThrustOrientation );
//            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
//        }
//        break;
//    }
////    case mee_costate_based_thrust_direction:
////    {
////        // Check input consistency
////        std::shared_ptr< MeeCostateBasedThrustDirectionSettings > meeCostateBasedThrustSettings =
////                std::dynamic_pointer_cast< MeeCostateBasedThrustDirectionSettings >( thrustDirectionGuidanceSettings );

////        if( meeCostateBasedThrustSettings == nullptr )
////        {
////            throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, input is inconsistent" );
////        }
////        else
////        {
////            // Check whether all required environment properties exist
////            if( bodies.count( meeCostateBasedThrustSettings->relativeBody_ ) == 0 )
////            {
////                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, central body " +
////                                          meeCostateBasedThrustSettings->relativeBody_ + " not found." );
////            }
////            else if( bodies.count( meeCostateBasedThrustSettings->vehicleName_ ) == 0 )
////            {
////                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, thrusting body " +
////                                          meeCostateBasedThrustSettings->vehicleName_ + " not found." );
////            }
////            else if( bodies.at( meeCostateBasedThrustSettings->relativeBody_ )->getGravityFieldModel( ) == nullptr )
////            {
////                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, central body " +
////                                          meeCostateBasedThrustSettings->relativeBody_ + " has no gravity field." );
////            }
////            else
////            {
////                // Retrieve required functions and create guidance object
////                std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction =
////                        std::bind( &simulation_setup::Body::getState,
////                                     bodies.at( meeCostateBasedThrustSettings->vehicleName_ ) );
////                std::function< Eigen::Vector6d( ) > centralBodyStateFunction =
////                        std::bind( &simulation_setup::Body::getState,
////                                     bodies.at( meeCostateBasedThrustSettings->relativeBody_ ) );
////                std::function< double( ) > centralBodyGravitationalParameterFunction =
////                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
////                                     bodies.at( meeCostateBasedThrustSettings->relativeBody_ )->getGravityFieldModel( ) );

////                thrustGuidance =  std::make_shared< propulsion::MeeCostateBasedThrustGuidance >(
////                            thrustingBodyStateFunction, centralBodyStateFunction,
////                            centralBodyGravitationalParameterFunction,
////                            meeCostateBasedThrustSettings->costateFunction_,
////                            bodyFixedThrustOrientation );
////            }
////        }
////        break;
////    }
//    default:
//        throw std::runtime_error( "Error, could not find thrust guidance type when creating thrust guidance." );
//    }
//    return thrustGuidance;
//}

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
    case constant_thrust_magnitude:
    {
        // Check input consistency
        std::shared_ptr< ConstantThrustMagnitudeSettings > constantThrustMagnitudeSettings =
                std::dynamic_pointer_cast< ConstantThrustMagnitudeSettings >( thrustMagnitudeSettings );
        if( constantThrustMagnitudeSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating constant thrust magnitude wrapper, input is inconsistent" );
        }

        thrustMagnitudeWrapper = std::make_shared< propulsion::ConstantThrustMagnitudeWrapper >(
                    constantThrustMagnitudeSettings->thrustMagnitude_,
                    constantThrustMagnitudeSettings->specificImpulse_ );
        break;

    }
    case thrust_magnitude_from_time_function:
    {
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
                    bodies.at( nameOfBodyWithGuidance ), parameterizedThrustMagnitudeSettings->thrustIndependentVariables_,
                    parameterizedThrustMagnitudeSettings->thrustGuidanceInputVariables_ );
        std::vector< std::function< double( ) > > specificInputVariableFunctions =
                getPropulsionInputVariables(
                    bodies.at( nameOfBodyWithGuidance ), parameterizedThrustMagnitudeSettings->specificImpulseDependentVariables_,
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
//    case bang_bang_thrust_magnitude_from_mee_costates:
//    {
//        // Check input consistency
//        std::shared_ptr< FromMeeCostatesBangBangThrustMagnitudeSettings > fromMeeCostatesBangBangThrustMagnitudeSettings =
//                std::dynamic_pointer_cast< FromMeeCostatesBangBangThrustMagnitudeSettings >( thrustMagnitudeSettings );
//        if( fromMeeCostatesBangBangThrustMagnitudeSettings == nullptr )
//        {
//            throw std::runtime_error( "Error when creating bang-bang thrust magnitude wrapper based on Mee costates, input is inconsistent" );
//        }
//        else
//        {
//            // Check whether all required environment properties exist
//            if( bodies.count( fromMeeCostatesBangBangThrustMagnitudeSettings->centralBodyName_ ) == 0 )
//            {
//                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, central body " +
//                                          fromMeeCostatesBangBangThrustMagnitudeSettings->centralBodyName_ + " not found." );
//            }
//            else if( bodies.count( fromMeeCostatesBangBangThrustMagnitudeSettings->vehicleName_ ) == 0 )
//            {
//                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, thrusting body " +
//                                          fromMeeCostatesBangBangThrustMagnitudeSettings->vehicleName_ + " not found." );
//            }
//            else if( bodies.at( fromMeeCostatesBangBangThrustMagnitudeSettings->centralBodyName_ )->getGravityFieldModel( ) == nullptr )
//            {
//                throw std::runtime_error( "Error when getting thrust guidance with mee_costate_based_thrust_direction, central body " +
//                                          fromMeeCostatesBangBangThrustMagnitudeSettings->centralBodyName_ + " has no gravity field." );
//            }
//            else
//            {
//                // Retrieve required functions and create guidance object
//                std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction =
//                        std::bind( &simulation_setup::Body::getState,
//                                     bodies.at( fromMeeCostatesBangBangThrustMagnitudeSettings->vehicleName_ ) );
//                std::function< double( ) > thrustingBodyMassFunction =
//                        std::bind( &simulation_setup::Body::getBodyMass,
//                                     bodies.at( fromMeeCostatesBangBangThrustMagnitudeSettings->vehicleName_ ) );
//                std::function< Eigen::Vector6d( ) > centralBodyStateFunction =
//                        std::bind( &simulation_setup::Body::getState,
//                                     bodies.at( fromMeeCostatesBangBangThrustMagnitudeSettings->centralBodyName_ ) );
//                std::function< double( ) > centralBodyGravitationalParameterFunction =
//                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
//                                     bodies.at( fromMeeCostatesBangBangThrustMagnitudeSettings->centralBodyName_ )->getGravityFieldModel( ) );

//                // Create thrust magnitude wrapper
//                thrustMagnitudeWrapper = std::make_shared< propulsion::MeeCostatesBangBangThrustMagnitudeWrapper >(
//                            thrustingBodyStateFunction,
//                            centralBodyStateFunction,
//                            centralBodyGravitationalParameterFunction,
//                            fromMeeCostatesBangBangThrustMagnitudeSettings->costatesFunction_,
//                            fromMeeCostatesBangBangThrustMagnitudeSettings->maximumThrustMagnitude_,
//                            fromMeeCostatesBangBangThrustMagnitudeSettings->specificImpulseFunction_,
//                            thrustingBodyMassFunction,
//                            fromMeeCostatesBangBangThrustMagnitudeSettings->customThrustResetFunction_ );
//            }
//        }

//        break;
//    }
    default:
        throw std::runtime_error( "Error when creating thrust magnitude wrapper, type not identified" );
    }

    return thrustMagnitudeWrapper;
}

////! Function to update the thrust magnitude and direction to current time.
//void updateThrustSettings(
//        const std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
//        const std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
//        const double currentTime )
//{
//    thrustMagnitudeWrapper->update( currentTime );
//    thrustDirectionGuidance->updateCalculator( currentTime );
//}

////! Function to reset the current time variable of the thrust magnitude and direction wrappers
//void resetThrustSettingsTime(
//        const std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
//        const std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
//        const double currentTime )
//{
//    thrustMagnitudeWrapper->resetCurrentTime( currentTime );
//    thrustDirectionGuidance->resetCurrentTime( currentTime );
//}


}

}

