/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"
#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/flightConditions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/simulation/environment_setup/body.h"

#include <iostream>

namespace tudat
{

namespace simulation_setup
{

RigidBodyProperties::RigidBodyProperties( ):
    currentMass_( TUDAT_NAN ), currentCenterOfMass_( Eigen::Vector3d::Constant( TUDAT_NAN ) ),
    currentInertiaTensor_( Eigen::Matrix3d::Constant( TUDAT_NAN ) ), isBodyInPropagation_( false ), isMassComputed_( false ),
    isComComputed_( false ), isInertiaTensorComputed_( false )
{ }

RigidBodyProperties::~RigidBodyProperties( ) { }

void RigidBodyProperties::update( const double currentTime )
{
    updateMass( currentTime );
    updateMassDistribution( currentTime );
}

void RigidBodyProperties::resetCurrentTime( )
{
    isMassComputed_ = false;
    isComComputed_ = false;
    isInertiaTensorComputed_ = false;
}

double RigidBodyProperties::getCurrentMass( )
{
    if( !isMassComputed_ )
    {
        throw std::runtime_error( "Error when retrieveing mass, mass is not computed/defined." );
    }
    return currentMass_;
}

Eigen::Vector3d RigidBodyProperties::getCurrentCenterOfMass( )
{
    if( !isComComputed_ )
    {
        throw std::runtime_error( "Error when retrieving center of mass, center of mass is not computed/defined." );
    }
    return currentCenterOfMass_;
}

Eigen::Matrix3d RigidBodyProperties::getCurrentInertiaTensor( )
{
    if( !isInertiaTensorComputed_ )
    {
        throw std::runtime_error( "Error when retrieving inertia tensor, inertia tensor is not computed/defined." );
    }
    return currentInertiaTensor_;
}

void RigidBodyProperties::setIsBodyInPropagation( const bool isBodyInPropagation )
{
    isBodyInPropagation_ = isBodyInPropagation;
}

TimeDependentRigidBodyProperties::TimeDependentRigidBodyProperties(
        const std::function< double( const double ) > massFunction,
        const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction,
        const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction ):
    massFunction_( massFunction ), centerOfMassFunction_( centerOfMassFunction ), inertiaTensorFunction_( inertiaTensorFunction )
{ }

TimeDependentRigidBodyProperties::TimeDependentRigidBodyProperties(
        const double mass, const Eigen::Vector3d& centerOfMass, const Eigen::Matrix3d& inertiaTensor ):
    massFunction_( nullptr ), centerOfMassFunction_( nullptr ), inertiaTensorFunction_( nullptr )
{
    currentMass_ = mass;
    isMassComputed_ = true;

    if( !centerOfMass.hasNaN( ) )
    {
        currentCenterOfMass_ = centerOfMass;
        isComComputed_ = true;
    }

    if( !inertiaTensor.hasNaN( ) )
    {
        currentInertiaTensor_ = inertiaTensor;
        isInertiaTensorComputed_ = true;
    }
}

TimeDependentRigidBodyProperties::~TimeDependentRigidBodyProperties( ) { }

void TimeDependentRigidBodyProperties::resetCurrentTime( )
{
    if( massFunction_ != nullptr )
    {
        isMassComputed_ = false;
    }

    if( centerOfMassFunction_ != nullptr )
    {
        isComComputed_ = false;
    }

    if( inertiaTensorFunction_ != nullptr )
    {
        isInertiaTensorComputed_ = false;
    }
}

std::function< double( const double ) > TimeDependentRigidBodyProperties::getMassFunction( )
{
    return massFunction_;
}

void TimeDependentRigidBodyProperties::setMassFunction( const std::function< double( const double ) > massFunction )
{
    massFunction_ = massFunction;
}

void TimeDependentRigidBodyProperties::setCurrentMass( const double currentMass )
{
    currentMass_ = currentMass;
    isMassComputed_ = true;
}

void TimeDependentRigidBodyProperties::setInertiaTensor( const Eigen::Matrix3d& inertiaTensor )
{
    currentInertiaTensor_ = inertiaTensor;
    isInertiaTensorComputed_ = true;
}

void TimeDependentRigidBodyProperties::updateMass( const double currentTime )
{
    if( massFunction_ != nullptr && ( !isMassComputed_ || !isBodyInPropagation_ ) )
    {
        currentMass_ = massFunction_( currentTime );
        isMassComputed_ = true;
    }
}

void TimeDependentRigidBodyProperties::updateMassDistribution( const double currentTime )
{
    if( centerOfMassFunction_ != nullptr && ( !isComComputed_ || !isBodyInPropagation_ ) )
    {
        currentCenterOfMass_ = centerOfMassFunction_( currentTime );
        isComComputed_ = true;
    }

    if( inertiaTensorFunction_ != nullptr && ( !isInertiaTensorComputed_ || !isBodyInPropagation_ ) )
    {
        currentInertiaTensor_ = inertiaTensorFunction_( currentTime );
        isInertiaTensorComputed_ = true;
    }
}

MassDependentRigidBodyProperties::MassDependentRigidBodyProperties(
        const double currentMass,
        const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction,
        const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction ):
    centerOfMassFunction_( centerOfMassFunction ), inertiaTensorFunction_( inertiaTensorFunction )
{
    setCurrentMass( currentMass );
    updateMassDistribution( TUDAT_NAN );
}

MassDependentRigidBodyProperties::~MassDependentRigidBodyProperties( ) { }

void MassDependentRigidBodyProperties::updateMass( const double currentTime ) { }

void MassDependentRigidBodyProperties::updateMassDistribution( const double currentTime )
{
    if( centerOfMassFunction_ != nullptr && ( !isComComputed_ || !isBodyInPropagation_ ) )
    {
        currentCenterOfMass_ = centerOfMassFunction_( currentMass_ );
        isComComputed_ = true;
    }

    if( inertiaTensorFunction_ != nullptr && ( !isInertiaTensorComputed_ || !isBodyInPropagation_ ) )
    {
        currentInertiaTensor_ = inertiaTensorFunction_( currentMass_ );
        isInertiaTensorComputed_ = true;
    }
}

void MassDependentRigidBodyProperties::setCurrentMass( const double currentMass )
{
    currentMass_ = currentMass;
    isMassComputed_ = true;
}

FromGravityFieldRigidBodyProperties::FromGravityFieldRigidBodyProperties(
        const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel ):
    gravityFieldModel_( gravityFieldModel )
{
    currentMass_ = gravityFieldModel_->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT;
    currentCenterOfMass_ = gravityFieldModel_->getCenterOfMass( );
    currentInertiaTensor_ = gravityFieldModel_->getInertiaTensor( );

    isMassComputed_ = true;
    isComComputed_ = true;
    isInertiaTensorComputed_ = true;

    modelIsTimeDependent_ =
            ( std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel ) != nullptr );
}

FromGravityFieldRigidBodyProperties::~FromGravityFieldRigidBodyProperties( ) { }

void FromGravityFieldRigidBodyProperties::resetCurrentTime( )
{
    if( modelIsTimeDependent_ )
    {
        isMassComputed_ = false;
        isComComputed_ = false;
        isInertiaTensorComputed_ = false;
    }
}

void FromGravityFieldRigidBodyProperties::updateMass( const double currentTime )
{
    if( ( modelIsTimeDependent_ && !isMassComputed_ ) || !isBodyInPropagation_ )
    {
        currentMass_ = gravityFieldModel_->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT;
        isMassComputed_ = true;
    }
}

void FromGravityFieldRigidBodyProperties::updateMassDistribution( const double currentTime )
{
    if( ( modelIsTimeDependent_ && ( !isComComputed_ || !isInertiaTensorComputed_ ) ) || !isBodyInPropagation_ )
    {
        currentCenterOfMass_ = gravityFieldModel_->getCenterOfMass( );
        currentInertiaTensor_ = gravityFieldModel_->getInertiaTensor( );
        isComComputed_ = true;
        isInertiaTensorComputed_ = true;
    }
}

void FromGravityFieldRigidBodyProperties::setCurrentMass( const double currentMass )
{
    throw std::runtime_error(
            "Error when resetting body mass; mass cannot be reset for bodies with a gravity field. Reset the gravity field's "
            "gravitational parameter instead." );
}

void FromGravityFieldRigidBodyProperties::setIsBodyInPropagation( const bool isBodyInPropagation )
{
    currentMass_ = gravityFieldModel_->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT;
    currentCenterOfMass_ = gravityFieldModel_->getCenterOfMass( );
    currentInertiaTensor_ = gravityFieldModel_->getInertiaTensor( );

    isBodyInPropagation_ = isBodyInPropagation;
}

Body::Body( const Eigen::Vector6d& state ):
    bodyIsGlobalFrameOrigin_( -1 ), currentState_( state ), timeOfCurrentState_( TUDAT_NAN ),
    ephemerisFrameToBaseFrame_( std::make_shared< BaseStateInterfaceImplementation< double, double > >( "" ) ),
    currentRotationToLocalFrame_( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
    currentRotationToGlobalFrame_( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
    currentRotationToLocalFrameDerivative_( Eigen::Matrix3d::Zero( ) ),
    currentAngularVelocityVectorInGlobalFrame_( Eigen::Vector3d::Zero( ) ),
    currentAngularVelocityVectorInLocalFrame_( Eigen::Vector3d::Zero( ) ), bodyName_( "unnamed_body" )
{
    currentLongState_ = currentState_.cast< long double >( );
    isStateSet_ = false;
    isRotationSet_ = false;
}

void Body::setCurrentRotationToLocalFrameFromEphemeris( const double time )
{
    if( rotationalEphemeris_ != nullptr )
    {
        currentRotationToLocalFrame_ = rotationalEphemeris_->getRotationToTargetFrame( time );
    }
    else
    {
        throw std::runtime_error( "Error, no rotation model found in Body::setCurrentRotationToLocalFrameFromEphemeris" );
    }
    currentRotationToGlobalFrame_ = currentRotationToLocalFrame_.inverse( );
    isRotationSet_ = true;
}

void Body::setCurrentRotationalStateToLocalFrame( const Eigen::Vector7d currentRotationalStateFromLocalToGlobalFrame )
{
    currentRotationToGlobalFrame_ = Eigen::Quaterniond( currentRotationalStateFromLocalToGlobalFrame( 0 ),
                                                        currentRotationalStateFromLocalToGlobalFrame( 1 ),
                                                        currentRotationalStateFromLocalToGlobalFrame( 2 ),
                                                        currentRotationalStateFromLocalToGlobalFrame( 3 ) );

    currentRotationToGlobalFrame_.normalize( );
    currentRotationToLocalFrame_ = currentRotationToGlobalFrame_.inverse( );
    currentAngularVelocityVectorInGlobalFrame_ =
            currentRotationToGlobalFrame_ * currentRotationalStateFromLocalToGlobalFrame.block< 3, 1 >( 4, 0 );

    currentAngularVelocityVectorInLocalFrame_ = currentRotationalStateFromLocalToGlobalFrame.block< 3, 1 >( 4, 0 );

    Eigen::Matrix3d currentRotationMatrixToLocalFrame = currentRotationToLocalFrame_.toRotationMatrix( );
    currentRotationToLocalFrameDerivative_ =
            linear_algebra::getCrossProductMatrix( currentRotationalStateFromLocalToGlobalFrame.block< 3, 1 >( 4, 0 ) ) *
            currentRotationMatrixToLocalFrame;
    isRotationSet_ = true;
}

void Body::setGravityFieldModel( const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel )
{
    gravityFieldModel_ = gravityFieldModel;

    if( massProperties_ != nullptr )
    {
        std::cerr << "Warning when settings gravity field model for body, mass interface already found: overrriding existing mass "
                     "interface"
                  << std::endl;
    }
    else
    {
        massProperties_ = std::make_shared< FromGravityFieldRigidBodyProperties >( gravityFieldModel );
    }
}

void Body::setAerodynamicCoefficientInterface(
        const std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > aerodynamicCoefficientInterface )
{
    aerodynamicCoefficientInterface_ = aerodynamicCoefficientInterface;
    std::vector< std::string > controlSurfaceList = aerodynamicCoefficientInterface->getControlSurfaceNames( );
    for( unsigned int i = 0; i < controlSurfaceList.size( ); i++ )
    {
        if( getVehicleSystems( )->doesControlSurfaceExist( controlSurfaceList.at( i ) ) == 0 )
        {
            vehicleSystems_->setCurrentControlSurfaceDeflection( controlSurfaceList.at( i ), 0.0 );
        }
    }

    if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >( aerodynamicFlightConditions_ ) != nullptr )
    {
        std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >( aerodynamicFlightConditions_ )
                ->resetAerodynamicCoefficientInterface( aerodynamicCoefficientInterface_ );
    }
}

void Body::setFlightConditions( const std::shared_ptr< aerodynamics::FlightConditions > aerodynamicFlightConditions )
{
    aerodynamicFlightConditions_ = aerodynamicFlightConditions;

    if( rotationalEphemeris_ != nullptr &&
        std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >( rotationalEphemeris_ ) == nullptr )
    {
        aerodynamicFlightConditions_->getAerodynamicAngleCalculator( )->setBodyFixedAngleInterface(
                std::make_shared< reference_frames::FromGenericEphemerisAerodynamicAngleInterface >( rotationalEphemeris_ ) );
    }
}

double Body::getGravitationalParameter( )
{
    if( gravityFieldModel_ == nullptr )
    {
        throw std::runtime_error( "Error when retrieveing gravitational parameter from body " + bodyName_ +
                                  ", no gravity field model is defined" );
    }
    return gravityFieldModel_->getGravitationalParameter( );
}

std::shared_ptr< system_models::VehicleSystems > Body::getVehicleSystems( )
{
    if( vehicleSystems_ == nullptr )
    {
        vehicleSystems_ = std::make_shared< system_models::VehicleSystems >( );
    }
    return vehicleSystems_;
}

void Body::setMassProperties( const std::shared_ptr< RigidBodyProperties > massProperties )
{
    if( gravityFieldModel_ != nullptr )
    {
        std::cerr << "Warning, setting body mass distribution, but existing gravity field model and associated mass properties already "
                     "found; overriding existing body mass properties"
                  << std::endl;
    }
    massProperties_ = massProperties;
}

void Body::setBodyMassFunction( const std::function< double( const double ) > bodyMassFunction )
{
    if( massProperties_ == nullptr )
    {
        massProperties_ = std::make_shared< TimeDependentRigidBodyProperties >( bodyMassFunction );
    }
    else if( std::dynamic_pointer_cast< TimeDependentRigidBodyProperties >( massProperties_ ) != nullptr )
    {
        std::dynamic_pointer_cast< TimeDependentRigidBodyProperties >( massProperties_ )->setMassFunction( bodyMassFunction );
    }
    else
    {
        throw std::runtime_error( "Error when resetting body mass function for " + bodyName_ + ", no compatible mass properties found" );
    }
}

void Body::setConstantBodyMass( const double bodyMass )
{
    if( massProperties_ == nullptr )
    {
        massProperties_ = std::make_shared< TimeDependentRigidBodyProperties >( bodyMass );
    }
    else
    {
        massProperties_->setCurrentMass( bodyMass );
    }
}

void Body::setCurrentPropagatedBodyMass( const double bodyMass )
{
    if( massProperties_ == nullptr )
    {
        massProperties_ = std::make_shared< TimeDependentRigidBodyProperties >( bodyMass );
    }
    else
    {
        massProperties_->setCurrentMass( bodyMass );
    }
}

std::function< double( const double ) > Body::getBodyMassFunction( )
{
    if( std::dynamic_pointer_cast< TimeDependentRigidBodyProperties >( massProperties_ ) != nullptr )
    {
        return std::dynamic_pointer_cast< TimeDependentRigidBodyProperties >( massProperties_ )->getMassFunction( );
    }
    throw std::runtime_error( "Error when getting body mass function for " + bodyName_ + ", no compatible mass properties found" );
}

void Body::updateMass( const double time )
{
    if( massProperties_ == nullptr )
    {
        throw std::runtime_error( "Error when updating body mass for " + bodyName_ + ", no mass properties found" );
    }
    massProperties_->updateMass( time );
}

void Body::updateMassDistribution( const double time )
{
    if( massProperties_ == nullptr )
    {
        throw std::runtime_error( "Error when updating body mass for " + bodyName_ + ", no mass properties found" );
    }
    massProperties_->updateMassDistribution( time );
}

double Body::getBodyMass( )
{
    if( massProperties_ == nullptr )
    {
        throw std::runtime_error( "Error when retrieving mass of " + bodyName_ + ", no mass properties found" );
    }
    return massProperties_->getCurrentMass( );
}

Eigen::Matrix3d Body::getBodyInertiaTensor( )
{
    if( massProperties_ == nullptr )
    {
        throw std::runtime_error( "Error when retrieving inertia tensor of " + bodyName_ + ", no mass properties found" );
    }
    return massProperties_->getCurrentInertiaTensor( );
}

void Body::setBodyInertiaTensor( const Eigen::Matrix3d& bodyInertiaTensor )
{
    if( massProperties_ == nullptr )
    {
        throw std::runtime_error( "Error when setting inertia tensor; no body mass properties object found (add a body mass first)" );
    }
    else if( std::dynamic_pointer_cast< TimeDependentRigidBodyProperties >( massProperties_ ) != nullptr )
    {
        std::dynamic_pointer_cast< TimeDependentRigidBodyProperties >( massProperties_ )->setInertiaTensor( bodyInertiaTensor );
    }
    else
    {
        throw std::runtime_error(
                "Error when trying to reset body inertia tensor, rigid body properties already exist, and are not of compatible "
                "type" );
    }
}

void Body::updateConstantEphemerisDependentMemberQuantities( )
{
    if( std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel_ ) != nullptr )
    {
        std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel_ )
                ->updateCorrectionFunctions( );
    }
}

void Body::getPositionByReference( Eigen::Vector3d& position )
{
    position = currentState_.segment( 0, 3 );
}

// template void Body::setStateFromEphemeris< double, double >( const double& time );

//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template<>
Eigen::Matrix< double, 6, 1 > BaseStateInterface::getBaseFrameState( const double time )
{
    return getBaseFrameDoubleState( time );
}

//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template<>
Eigen::Matrix< long double, 6, 1 > BaseStateInterface::getBaseFrameState( const double time )
{
    return getBaseFrameLongDoubleState( time );
}

//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template<>
Eigen::Matrix< double, 6, 1 > BaseStateInterface::getBaseFrameState( const Time time )
{
    return getBaseFrameDoubleState( time );
}

//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template<>
Eigen::Matrix< long double, 6, 1 > BaseStateInterface::getBaseFrameState( const Time time )
{
    return getBaseFrameLongDoubleState( time );
}

template<>
Eigen::Matrix< double, 6, 1 > Body::getTemplatedState( )
{
    return getState( );
}

template<>
Eigen::Matrix< long double, 6, 1 > Body::getTemplatedState( )
{
    return getLongState( );
}

//! Templated function to set the state manually.
template<>
void Body::setTemplatedState( const Eigen::Matrix< double, 6, 1 >& state )
{
    setState( state );
}

//! Templated function to set the state manually.
template<>
void Body::setTemplatedState( const Eigen::Matrix< long double, 6, 1 >& state )
{
    setLongState( state );
}

//! Function to define whether the body is currently being propagated, or not
void Body::setIsBodyInPropagation( const bool isBodyInPropagation )
{
    isBodyInPropagation_ = isBodyInPropagation;

    if( rotationalEphemeris_ != nullptr )
    {
        rotationalEphemeris_->setIsBodyInPropagation( isBodyInPropagation );
    }

    if( massProperties_ != nullptr )
    {
        massProperties_->setIsBodyInPropagation( isBodyInPropagation );
    }

    if( !isBodyInPropagation )
    {
        isStateSet_ = false;
        isRotationSet_ = false;
    }
}

double getBodyGravitationalParameter( const SystemOfBodies& bodies, const std::string bodyName )
{
    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when getting gravitational parameter of body " + bodyName + ", no such body is found" );
    }
    else if( bodies.at( bodyName )->getGravityFieldModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when getting gravitational parameter of body " + bodyName + ", body has not gravity field" );
    }
    return bodies.at( bodyName )->getGravityFieldModel( )->getGravitationalParameter( );
}

//! Function ot retrieve the common global translational state origin of the environment
std::string getGlobalFrameOrigin( const SystemOfBodies& bodies )
{
    std::string globalFrameOrigin = "SSB";

    for( auto bodyIterator: bodies.getMap( ) )
    {
        if( bodyIterator.second->getIsBodyGlobalFrameOrigin( ) == -1 )
        {
            throw std::runtime_error( "Error, body " + bodyIterator.first + " does not have global frame origin set" );
        }
        else if( bodyIterator.second->getIsBodyGlobalFrameOrigin( ) == 1 )
        {
            if( globalFrameOrigin != "SSB" )
            {
                throw std::runtime_error( "Error, body " + bodyIterator.first + " found as global frame origin, but body " +
                                          globalFrameOrigin + " has already been detected as global frame origin." );
            }
            else
            {
                globalFrameOrigin = bodyIterator.first;
            }
        }
    }
    return globalFrameOrigin;
}

std::shared_ptr< ephemerides::ReferenceFrameManager > createFrameManager(
        const std::unordered_map< std::string, std::shared_ptr< Body > > bodies )
{
    // Get ephemerides from bodies
    std::map< std::string, std::shared_ptr< ephemerides::Ephemeris > > ephemerides;
    for( auto bodyIterator: bodies )
    {
        if( bodyIterator.second->getEphemeris( ) != nullptr )
        {
            ephemerides[ bodyIterator.first ] = bodyIterator.second->getEphemeris( );
        }
    }
    return std::make_shared< ephemerides::ReferenceFrameManager >( ephemerides );
}

//! Function to set whether the bodies are currently being propagated, or not
void setAreBodiesInPropagation( const SystemOfBodies& bodies, const bool areBodiesInPropagation )
{
    for( auto bodyIterator: bodies.getMap( ) )
    {
        bodyIterator.second->setIsBodyInPropagation( areBodiesInPropagation );
    }
}

bool isReferencePointGroundStation( const std::shared_ptr< Body > body, const std::string& referencePointName )
{
    bool isReferencePointGroundStation = false;
    if( body->getGroundStationMap( ).count( referencePointName ) > 0 )
    {
        isReferencePointGroundStation = true;
    }
    else
    {
        if( body->getVehicleSystems( ) == nullptr )
        {
            throw std::runtime_error( "Error when finding reference point " + referencePointName + " on " + body->getBodyName( ) +
                                      " , point is not a ground station, and no system models found" );
        }
        else if( !body->getVehicleSystems( )->doesReferencePointExist( referencePointName ) )
        {
            throw std::runtime_error( "Error when finding reference point " + referencePointName + " on " + body->getBodyName( ) +
                                      ", point is not a ground station, and not a system reference point" );
        }
        else
        {
            isReferencePointGroundStation = false;
        }
    }
    return isReferencePointGroundStation;
}

bool isReferencePointGroundStation( const SystemOfBodies& bodies, const std::string& bodyName, const std::string& referencePointName )
{
    return isReferencePointGroundStation( bodies.at( bodyName ), referencePointName );
}

std::shared_ptr< system_models::TimingSystem > getTimingSystem( const std::pair< std::string, std::string > linkEndName,
                                                                const SystemOfBodies& bodyMap )
{
    std::shared_ptr< system_models::TimingSystem > timingSystem = nullptr;

    if( bodyMap.count( linkEndName.first ) > 0 )
    {
        std::shared_ptr< Body > currentBody = bodyMap.at( linkEndName.first );

        if( currentBody->getVehicleSystems( ) != NULL )
        {
            timingSystem = currentBody->getVehicleSystems( )->getTimingSystem( );
        }
        if( timingSystem == nullptr )
        {
            if( currentBody->getGroundStationMap( ).count( linkEndName.second ) > 0 )
            {
                std::shared_ptr< ground_stations::GroundStation > currentGroundStation =
                        currentBody->getGroundStation( linkEndName.second );
                if( currentGroundStation->getTimingSystem( ) != NULL )
                {
                    timingSystem = currentGroundStation->getTimingSystem( );
                }
            }
        }
    }

    if( timingSystem == nullptr )
    {
        throw std::runtime_error( "Error, did not find timing system for +(" + linkEndName.first + "," + linkEndName.second + ")" );
    }

    return timingSystem;
}
}  // namespace simulation_setup

}  // namespace tudat
