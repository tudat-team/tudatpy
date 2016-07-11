#ifndef CREATETHRUSTMODELGUIDANCE_H
#define CREATETHRUSTMODELGUIDANCE_H


#include "Tudat/Astrodynamics/BasicAstrodynamics/thrustGuidance.h"
#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

namespace tudat
{

namespace simulation_setup
{

enum ThrustDirectionGuidanceTypes
{
    colinear_with_state_segment_thrust,
    thrust_direction_from_body_orientation
};

class ThrustDirectionGuidanceSettings
{
public:
   ThrustDirectionGuidanceSettings(
           const ThrustDirectionGuidanceTypes thrustDirectionType,
           const std::string relativeBody ):
   thrustDirectionType_( thrustDirectionType ), relativeBody_( relativeBody ){ }

   ThrustDirectionGuidanceTypes thrustDirectionType_;

   std::string relativeBody_;
};

class ThrustDirectionFromStateGuidanceSettings: public ThrustDirectionGuidanceSettings
{
public:
   ThrustDirectionFromStateGuidanceSettings(
           const std::string& centralBody,
           const bool isColinearWithVelocity,
           const bool directionIsOppositeToVector ):
       ThrustDirectionGuidanceSettings( colinear_with_state_segment_thrust, centralBody ),
   isColinearWithVelocity_( isColinearWithVelocity ),
   directionIsOppositeToVector_( directionIsOppositeToVector ){ }

   bool isColinearWithVelocity_;

   bool directionIsOppositeToVector_;

};

boost::shared_ptr< basic_astrodynamics::ThrustGuidance > createThrustGuidanceModel(
        const boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
        const NamedBodyMap& bodyMap,
        const std::string& bodyWithGuidance )
{
   boost::shared_ptr< basic_astrodynamics::ThrustGuidance > thrustGuidance;

   switch( thrustDirectionGuidanceSettings->thrustDirectionType_ )
   {
   case colinear_with_state_segment_thrust:
   {
        boost::shared_ptr< ThrustDirectionFromStateGuidanceSettings > thrustDirectionFromStateGuidanceSettings =
                boost::dynamic_pointer_cast< ThrustDirectionFromStateGuidanceSettings >( thrustDirectionGuidanceSettings );
        if( thrustDirectionFromStateGuidanceSettings == NULL )
        {

        }
        else
        {
            boost::function< basic_mathematics::Vector6d( ) > bodyStateFunction =
                    boost::bind( &Body::getStateInBaseFrameFromEphemeris, bodyMap.at( bodyWithGuidance ) );
            boost::function< basic_mathematics::Vector6d( ) > centralBodyStateFunction =
                    boost::bind( &Body::getStateInBaseFrameFromEphemeris, bodyMap.at(
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
   case thrust_direction_from_body_orientation:
   {
       boost::shared_ptr< Body > bodyWithGuidance = bodyMap.at( bodyWithGuidance );
       boost::shared_ptr< Body > relativeBody = bodyMap.at( thrustDirectionGuidanceSettings->relativeBody_ );


       // Retrieve flight conditions; create object if not yet extant.
       boost::shared_ptr< aerodynamics::FlightConditions > bodyFlightConditions =
               bodyWithGuidance->getFlightConditions( );
       if( bodyFlightConditions == NULL )
       {
           bodyWithGuidance->setFlightConditions(
                       createFlightConditions( bodyWit,hGuidance,
                                               relativeBody,
                                               bodyWithGuidance,
                                               relativeBody ) );
           bodyFlightConditions = bodyWithGuidance->getFlightConditions( );
       }

       boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator =
              bodyFlightConditions->getAerodynamicAngleCalculator( );
       boost::function< Eigen::Quaterniond( ) > rotationFunction =
               boost::bind( &reference_frames::AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                            angleCalculator, reference_frames::body_frame, reference_frames::inertial_frame );
       thrustGuidance =  boost::make_shared< basic_astrodynamics::DirectOrientationBasedThrustGuidance >(
                   rotationFunction );
       break;
   }
   default:
       throw std::runtime_error( "Error, could not find thrust guidance type when creating thrust guidance." );
   }
   return thrustGuidance;
}

}

}

#endif // CREATETHRUSTMODELGUIDANCE_H
