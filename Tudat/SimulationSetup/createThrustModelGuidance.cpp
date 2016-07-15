#include "Tudat/SimulationSetup/createThrustModelGuidance.h"
#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace simulation_setup
{


boost::shared_ptr< propulsion::ThrustDirectionGuidance > createThrustGuidanceModel(
        const boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
{
    boost::shared_ptr< propulsion::ThrustDirectionGuidance > thrustGuidance;

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

            boost::function< basic_mathematics::Vector6d( ) > stateFunction =
                    boost::bind(
                        &ephemerides::getRelativePosition, bodyStateFunction, centralBodyStateFunction );
            boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction;

            if( thrustDirectionFromStateGuidanceSettings->isColinearWithVelocity_ )
            {
                thrustDirectionFunction =
                        boost::bind( &propulsion::getThrustDirectionColinearWithVelocity, _1, _2,
                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
            }
            else
            {
                thrustDirectionFunction =
                        boost::bind( &propulsion::getThrustDirectionColinearWithPosition, _1, _2,
                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
            }

            thrustGuidance =  boost::make_shared< propulsion::StateBasedThrustGuidance >(
                        thrustDirectionFunction, stateFunction, thrustDirectionFromStateGuidanceSettings->relativeBody_ );
        }
        break;
   }
   case thrust_direction_from_existing_body_orientation:
   {
       boost::shared_ptr< Body > bodyWithGuidance = bodyMap.at( nameOfBodyWithGuidance );
       boost::shared_ptr< Body > relativeBody = bodyMap.at( thrustDirectionGuidanceSettings->relativeBody_ );


       boost::function< Eigen::Quaterniond( const double ) > rotationFunction;
       if( bodyWithGuidance->getFlightConditions( ) != NULL )
       {
           boost::shared_ptr< aerodynamics::FlightConditions > bodyFlightConditions =
                   bodyWithGuidance->getFlightConditions( );
           boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator =
                  bodyFlightConditions->getAerodynamicAngleCalculator( );
           rotationFunction =
                   boost::bind( &reference_frames::AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                                angleCalculator, reference_frames::body_frame, reference_frames::inertial_frame );

           magnitudeUpdateSettings[ propagators::vehicle_flight_conditions_update ].push_back( nameOfBodyWithGuidance );
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

       thrustGuidance =  boost::make_shared< propulsion::DirectOrientationBasedThrustGuidance >(
                   rotationFunction );
       break;
    }
    case custom_thrust_direction:
    {

        boost::shared_ptr< CustomThrustDirectionSettings > customThrustGuidanceSettings =
                boost::dynamic_pointer_cast< CustomThrustDirectionSettings >( thrustDirectionGuidanceSettings );
        if( customThrustGuidanceSettings == NULL )
        {
            throw std::runtime_error( "Error when getting thrust guidance, input is inconsistent" );
        }
        else
        {
            boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction =
                    boost::bind( &propulsion::getThrustDirectionFromTimeOnlyFunction,
                                 _1, _2, customThrustGuidanceSettings->thrustDirectionFunction_ );


            thrustGuidance =  boost::make_shared< propulsion::StateBasedThrustGuidance >(
                        thrustDirectionFunction, boost::lambda::constant( basic_mathematics::Vector6d::Zero( ) ), "" );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error, could not find thrust guidance type when creating thrust guidance." );
    }
    return thrustGuidance;
}

boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > createThrustMagnitudeWrapper(
        const boost::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings )
{
    boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper;
    switch( thrustMagnitudeSettings->thrustMagnitudeGuidanceType_ )
    {
    case constant_thrust_magnitude:
    {
        boost::shared_ptr< ConstantThrustMagnitudeSettings > constantThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< ConstantThrustMagnitudeSettings >( thrustMagnitudeSettings );
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
        boost::shared_ptr< FromEngineThrustMagnitudeSettings > fromEngineThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< FromEngineThrustMagnitudeSettings >( thrustMagnitudeSettings );

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
        boost::shared_ptr< FromFunctionThrustMagnitudeSettings > fromFunctionThrustMagnitudeSettings =
                boost::dynamic_pointer_cast< FromFunctionThrustMagnitudeSettings >( thrustMagnitudeSettings );
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
        const boost::shared_ptr< propulsion::ThrustDirectionGuidance > thrustDirectionGuidance,
        const double currentTime )
{
    thrustMagnitudeWrapper->update( currentTime );
    thrustDirectionGuidance->updateCalculator( currentTime );
}


}

}

