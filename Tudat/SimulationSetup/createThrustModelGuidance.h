#ifndef CREATETHRUSTMODELGUIDANCE_H
#define CREATETHRUSTMODELGUIDANCE_H

#include "Tudat/Astrodynamics/SystemModels/engineModel.h"
#include "Tudat/Astrodynamics/Propulsion/thrustGuidance.h"
#include "Tudat/Astrodynamics/Propagators/environmentUpdateTypes.h"
#include "Tudat/Astrodynamics/Propulsion/thrustMagnitudeWrapper.h"
#include "Tudat/SimulationSetup/body.h"
#include "Tudat/SimulationSetup/createFlightConditions.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/SystemModels/engineModel.h"

namespace tudat
{

namespace simulation_setup
{

enum ThrustDirectionGuidanceTypes
{
    colinear_with_state_segment_thrust_direction,
    thrust_direction_from_existing_body_orientation,
    custom_thrust_direction,
    custom_thrust_orientation

};

class ThrustDirectionGuidanceSettings
{
public:
   ThrustDirectionGuidanceSettings(
           const ThrustDirectionGuidanceTypes thrustDirectionType,
           const std::string relativeBody ):
   thrustDirectionType_( thrustDirectionType ), relativeBody_( relativeBody ){ }

   virtual ~ThrustDirectionGuidanceSettings( ){ }

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
       ThrustDirectionGuidanceSettings( colinear_with_state_segment_thrust_direction, centralBody ),
   isColinearWithVelocity_( isColinearWithVelocity ),
   directionIsOppositeToVector_( directionIsOppositeToVector ){ }

   ~ThrustDirectionFromStateGuidanceSettings( ){ }

   bool isColinearWithVelocity_;

   bool directionIsOppositeToVector_;

};

class CustomThrustDirectionSettings: public ThrustDirectionGuidanceSettings
{
public:
   CustomThrustDirectionSettings(
           const boost::function< Eigen::Vector3d( const double ) > thrustDirectionFunction ):
       ThrustDirectionGuidanceSettings( custom_thrust_direction, "" ),
     thrustDirectionFunction_( thrustDirectionFunction ){ }

   ~CustomThrustDirectionSettings( ){ }

   boost::function< Eigen::Vector3d( const double ) > thrustDirectionFunction_;
};

class CustomThrustOrientationSettings: public ThrustDirectionGuidanceSettings
{
public:
    CustomThrustOrientationSettings(
            const boost::function< Eigen::Quaterniond( const double ) > thrustOrientationFunction,
            const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) ):
        ThrustDirectionGuidanceSettings( custom_thrust_orientation, "" ),
        thrustOrientationFunction_( thrustOrientationFunction ),
        bodyFixedThrustDirection_( bodyFixedThrustDirection ){ }

    ~CustomThrustOrientationSettings( ){ }

    boost::function< Eigen::Quaterniond( const double ) > thrustOrientationFunction_ ;

    Eigen::Vector3d bodyFixedThrustDirection_;
};


boost::shared_ptr< propulsion::ThrustDirectionGuidance > createThrustGuidanceModel(
        const boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings );

enum ThrustMagnitudeTypes
{
    constant_thrust_magnitude,
    from_engine_properties_thrust_magnitude,
    thrust_magnitude_from_time_function
};

class ThrustMagnitudeSettings
{
public:
   ThrustMagnitudeSettings(
           const ThrustMagnitudeTypes thrustMagnitudeGuidanceType,
           const std::string& thrustOriginId ):
   thrustMagnitudeGuidanceType_( thrustMagnitudeGuidanceType ),
   thrustOriginId_( thrustOriginId ){ }

   virtual ~ThrustMagnitudeSettings( ){ }

   ThrustMagnitudeTypes thrustMagnitudeGuidanceType_;

   std::string thrustOriginId_;
};

class ConstantThrustMagnitudeSettings: public ThrustMagnitudeSettings
{
public:
   ConstantThrustMagnitudeSettings(
           const double thrustMagnitude,
           const double specificImpulse ):
       ThrustMagnitudeSettings( constant_thrust_magnitude, "" ),
   thrustMagnitude_( thrustMagnitude ), specificImpulse_( specificImpulse ){ }

   ~ConstantThrustMagnitudeSettings( ){ }

   double thrustMagnitude_;

   double specificImpulse_;
};


class FromEngineThrustMagnitudeSettings: public ThrustMagnitudeSettings
{
public:
    FromEngineThrustMagnitudeSettings(
            const bool useAllEngines = 1,
            const std::string& thrustOrigin = "" ):
        ThrustMagnitudeSettings( from_engine_properties_thrust_magnitude, thrustOrigin ),
        useAllEngines_( useAllEngines ){ }

    bool useAllEngines_;
};

class FromFunctionThrustMagnitudeSettings: public ThrustMagnitudeSettings
{
public:
    FromFunctionThrustMagnitudeSettings(
            const boost::function< double( const double ) > thrustMagnitudeFunction,
            const boost::function< double( const double ) > specificImpulseFunction,
            const boost::function< bool( const double ) > isEngineOnFunction = boost::lambda::constant( true ) ):
        ThrustMagnitudeSettings( thrust_magnitude_from_time_function, "" ),
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
   isEngineOnFunction_( isEngineOnFunction ){ }

   ~FromFunctionThrustMagnitudeSettings( ){ }


   boost::function< double( const double) > thrustMagnitudeFunction_;

   boost::function< double( const double) > specificImpulseFunction_;

   boost::function< bool( const double ) > isEngineOnFunction_;
};


boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > createThrustMagnitudeWrapper(
        const boost::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings );

void updateThrustMagnitudeAndDirection(
        const boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const boost::shared_ptr< propulsion::ThrustDirectionGuidance > thrustDirectionGuidance,
        const double currentTime );


}

}

#endif // CREATETHRUSTMODELGUIDANCE_H
