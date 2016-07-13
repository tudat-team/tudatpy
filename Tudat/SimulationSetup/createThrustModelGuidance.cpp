#include "Tudat/SimulationSetup/createThrustModelGuidance.h"

namespace tudat
{

namespace simulation_setup
{


boost::shared_ptr< basic_astrodynamics::ThrustDirectionGuidance > createThrustGuidanceModel(
        const boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance )
{
   boost::shared_ptr< basic_astrodynamics::ThrustDirectionGuidance > thrustGuidance;

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
            boost::function< basic_mathematics::Vector6d( ) > centralBodyStateFunction =
                    boost::bind( &Body::getState, bodyMap.at(
                                     thrustDirectionFromStateGuidanceSettings->relativeBody_ ) );

            boost::function< basic_mathematics::Vector6d( ) > stateFunction =
                    boost::bind(
                        &ephemerides::getRelativePosition, bodyStateFunction, centralBodyStateFunction );
            boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction;

            if( thrustDirectionFromStateGuidanceSettings->isColinearWithVelocity_ )
            {
                thrustDirectionFunction =
                        boost::bind( &basic_astrodynamics::getThrustDirectionColinearWithVelocity, _1, _2,
                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
            }
            else
            {
                thrustDirectionFunction =
                        boost::bind( &basic_astrodynamics::getThrustDirectionColinearWithPosition, _1, _2,
                                     thrustDirectionFromStateGuidanceSettings->directionIsOppositeToVector_ );
            }

            thrustGuidance =  boost::make_shared< basic_astrodynamics::StateBasedThrustGuidance >(
                        thrustDirectionFunction, stateFunction );
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
           bodyWithGuidance->setFlightConditions(
                       createFlightConditions( bodyWithGuidance,
                                               relativeBody,
                                               nameOfBodyWithGuidance,
                                               thrustDirectionGuidanceSettings->relativeBody_ ) );
           bodyFlightConditions = bodyWithGuidance->getFlightConditions( );
           boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator =
                  bodyFlightConditions->getAerodynamicAngleCalculator( );
           rotationFunction =
                   boost::bind( &reference_frames::AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                                angleCalculator, reference_frames::body_frame, reference_frames::inertial_frame );
       }
       else if( bodyWithGuidance->getRotationalEphemeris( ) != NULL )
       {
           rotationFunction = boost::bind(
                       &ephemerides::RotationalEphemeris::getRotationToBaseFrame,
                       bodyWithGuidance->getRotationalEphemeris( ), _1, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
       }
       else
       {
            throw std::runtime_error( "Error, requested thrust orientation from existing model, but no such model found" );
       }

       thrustGuidance =  boost::make_shared< basic_astrodynamics::DirectOrientationBasedThrustGuidance >(
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
                       boost::bind( &basic_astrodynamics::getThrustDirectionFromTimeOnlyFunction,
                                         _1, _2, customThrustGuidanceSettings->thrustDirectionFunction_ );


               thrustGuidance =  boost::make_shared< basic_astrodynamics::StateBasedThrustGuidance >(
                           thrustDirectionFunction, boost::lambda::constant( basic_mathematics::Vector6d::Zero( ) ) );
       }
       break;
   }
   default:
       throw std::runtime_error( "Error, could not find thrust guidance type when creating thrust guidance." );
   }
   return thrustGuidance;
}

boost::shared_ptr< ThrustMagnitudeWrapper > createThrustMagnitudeWrapper(
        const boost::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance )
{
    boost::shared_ptr< ThrustMagnitudeWrapper > thrustMagnitudeWrapper;
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
        thrustMagnitudeWrapper = boost::make_shared< CustomThrustMagnitudeWrapper >(
                    boost::lambda::constant( constantThrustMagnitudeSettings->thrustMagnitude_ ),
                    boost::lambda::constant( constantThrustMagnitudeSettings->specificImpulse_ ) );
        break;

    }
    case from_engine_properties_thrust_magnitude:
    {
        if( bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( ) == NULL )
        {
            throw std::runtime_error( "Error when creating from-engine thrust magnitude wrapper, no vehicle systems found" );

        }

        if( bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( )->getEngineModels( ).count(
                    thrustMagnitudeSettings->thrustOriginId_ ) == 0 )
        {
            throw std::runtime_error( "Error when creating from-engine thrust magnitude wrapper, no engine of right ID found" );
        }

        thrustMagnitudeWrapper = boost::make_shared< ThrustMagnitudeFromEngineWrapper >(
                    bodyMap.at( nameOfBodyWithGuidance )->getVehicleSystems( )->getEngineModels( ).at(
                        thrustMagnitudeSettings->thrustOriginId_ ) );
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
        thrustMagnitudeWrapper = boost::make_shared< CustomThrustMagnitudeWrapper >(
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
        const boost::shared_ptr< ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const boost::shared_ptr< basic_astrodynamics::ThrustDirectionGuidance > thrustDirectionGuidance,
        const double currentTime )
{
    thrustMagnitudeWrapper->update( currentTime );
    thrustDirectionGuidance->updateCalculator( currentTime );
}


}

}

