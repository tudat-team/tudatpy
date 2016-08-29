#include "Tudat/SimulationSetup/createThrustModelGuidance.h"
#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace simulation_setup
{


boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > createThrustGuidanceModel(
        const boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        const boost::function< Eigen::Vector3d( ) > bodyFixedThrustOrientation,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
{
    boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustGuidance;

    switch( thrustDirectionGuidanceSettings->thrustDirectionType_ )
    {
    case colinear_with_state_segment_thrust_direction:
    {
        boost::shared_ptr< ThrustDirectionFromStateGuidanceSettings > thrustDirectionFromStateGuidanceSettings =
                boost::dynamic_pointer_cast< ThrustDirectionFromStateGuidanceSettings >( thrustDirectionGuidanceSettings );
        if( thrustDirectionFromStateGuidanceSettings == NULL )
        {
            throw std::runtime_error( "Error when thrust colinear with state, input is inconsistent" );
        }
        else
        {
            boost::function< basic_mathematics::Vector6d( ) > bodyStateFunction =
                    boost::bind( &Body::getState, bodyMap.at( nameOfBodyWithGuidance ) );
            boost::function< basic_mathematics::Vector6d( ) > centralBodyStateFunction;

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
                centralBodyStateFunction = boost::lambda::constant( basic_mathematics::Vector6d::Zero( ) );
            }

            boost::function< void( basic_mathematics::Vector6d& ) > stateFunction =
                    boost::bind(
                        &ephemerides::getRelativeState, _1, bodyStateFunction, centralBodyStateFunction );
            boost::function< Eigen::Vector3d( const double ) > thrustDirectionFunction;

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

        boost::shared_ptr< CustomThrustDirectionSettings > customThrustGuidanceSettings =
                boost::dynamic_pointer_cast< CustomThrustDirectionSettings >( thrustDirectionGuidanceSettings );
        if( customThrustGuidanceSettings == NULL )
        {
            throw std::runtime_error( "Error when getting thrust guidance with custom_thrust_direction, input is inconsistent" );
        }
        else
        {
            boost::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
                    customThrustGuidanceSettings->thrustDirectionFunction_;


            thrustGuidance =  boost::make_shared< propulsion::DirectionBasedForceGuidance >(
                        thrustDirectionFunction, "", bodyFixedThrustOrientation );
        }

        break;
    }
    case custom_thrust_orientation:
    {
        boost::shared_ptr< CustomThrustOrientationSettings > customThrustOrientationSettings =
                boost::dynamic_pointer_cast< CustomThrustOrientationSettings >( thrustDirectionGuidanceSettings );

        if( customThrustOrientationSettings == NULL )
        {
            throw std::runtime_error( "Error when getting thrust guidance with custom_thrust_orientation, input is inconsistent" );
        }
        else
        {
            thrustGuidance =  boost::make_shared< propulsion::OrientationBasedForceGuidance >(
                        customThrustOrientationSettings->thrustOrientationFunction_,
                        bodyFixedThrustOrientation );
            magnitudeUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyWithGuidance );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error, could not find thrust guidance type when creating thrust guidance." );
    }
    return thrustGuidance;
}

Eigen::Vector3d getCombinedThrustDirection(
        const std::vector< boost::function< Eigen::Vector3d( )> > thrustDirections )
{
    Eigen::Vector3d thrustDirection = thrustDirections.at( 0 )( );
    for( unsigned int i = 1; i < thrustDirections.size( ); i++ )
    {
        if( thrustDirections.at( i )( ) != thrustDirection )
        {
            throw std::runtime_error( "Error, cannot have independently vectored engines in combined thurst model" );
        }
    }
    return thrustDirection;
}


boost::function< Eigen::Vector3d( ) > getBodyFixedThrustDirection(
        const boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string bodyName )
{
    boost::function< Eigen::Vector3d( ) > thrustDirectionFunction;
    switch( thrustMagnitudeSettings->thrustMagnitudeGuidanceType_ )
    {
    case constant_thrust_magnitude:
    {
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


        if( fromEngineThrustMagnitudeSettings->useAllEngines_ == false  )
        {
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
        else
        {
            if( ( bodyMap.at( bodyName )->getVehicleSystems( )->getEngineModels( ).size( ) == 0 ) )
            {
                std::cerr<<"Error when creating body-fixed thrust direction of type from_engine_properties_thrust_magnitude; no engines found: returning 0 thrust"<<std::endl;
            }

            std::vector< boost::function< Eigen::Vector3d( )> > thrustDirections;
            std::map< std::string, boost::shared_ptr< system_models::EngineModel > > engineModels =
                bodyMap.at( bodyName )->getVehicleSystems( )->getEngineModels( );
            for( std::map< std::string, boost::shared_ptr< system_models::EngineModel > >::const_iterator engineIterator =
                 engineModels.begin( ); engineIterator != engineModels.end( ); engineIterator++ )
            {
                thrustDirections.push_back(
                            boost::bind( &system_models::EngineModel::getBodyFixedThrustDirection, engineIterator->second ) );
            }

            thrustDirectionFunction = boost::bind(
                        &getCombinedThrustDirection, thrustDirections );

        }
        break;

    }
    case thrust_magnitude_from_time_function:
    {
        boost::shared_ptr< FromFunctionThrustEngineSettings > fromFunctionThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< FromFunctionThrustEngineSettings >( thrustMagnitudeSettings );
        if( fromFunctionThrustMagnitudeSettings == NULL )
        {
            throw std::runtime_error( "Error when creating body-fixed thrust direction of type thrust_magnitude_from_time_function, input is inconsistent" );
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

boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > createThrustMagnitudeWrapper(
        const boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
{
    boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper;
    switch( thrustMagnitudeSettings->thrustMagnitudeGuidanceType_ )
    {
    case constant_thrust_magnitude:
    {
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
        else
        {
            if( ( bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( )->getEngineModels( ).size( ) == 0 ) )
            {
                std::cerr<<"Error when creating from-engine thrust magnitude wrapper for all engines; no engines found: returning 0 thrust"<<std::endl;
            }
            thrustMagnitudeWrapper = boost::make_shared< propulsion::ThrustMagnitudeFromEngineWrapper >(
                        utilities::createVectorFromMapValues< boost::shared_ptr< system_models::EngineModel >, std::string >(
                                    bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( )->getEngineModels( ) ));
        }
        break;

    }
    case thrust_magnitude_from_time_function:
    {
        boost::shared_ptr< FromFunctionThrustEngineSettings > fromFunctionThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< FromFunctionThrustEngineSettings >( thrustMagnitudeSettings );
        if( fromFunctionThrustMagnitudeSettings == NULL )
        {
            throw std::runtime_error( "Error when creating from-function thrust magnitude wrapper, input is inconsistent" );
        }
        thrustMagnitudeWrapper = boost::make_shared< propulsion::CustomThrustMagnitudeWrapper >(
                    fromFunctionThrustMagnitudeSettings->thrustMagnitudeFunction_,
                    fromFunctionThrustMagnitudeSettings->specificImpulseFunction_,
                    fromFunctionThrustMagnitudeSettings->isEngineOnFunction_ );
        break;

    }
    default:
        throw std::runtime_error( "Error when creating thrust magnitude wrapper, type not identified" );
    }

    return thrustMagnitudeWrapper;
}

void updateThrustMagnitudeAndDirection(
        const boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
        const double currentTime )
{
    thrustMagnitudeWrapper->update( currentTime );
    thrustDirectionGuidance->updateCalculator( currentTime );
}

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

