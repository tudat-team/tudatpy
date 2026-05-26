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
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/basics/tudatExceptions.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/simulation/environment_setup/body.h"

#include <iostream>

namespace tudat
{

namespace simulation_setup
{

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

std::shared_ptr< BaseStateInterface > Body::getEphemerisFrameToBaseFrame( )
{
    return ephemerisFrameToBaseFrame_;
}

void Body::setEphemerisFrameToBaseFrame( const std::shared_ptr< BaseStateInterface > ephemerisFrameToBaseFrame )
{
    ephemerisFrameToBaseFrame_ = ephemerisFrameToBaseFrame;
}

Eigen::Vector6d Body::getState( )
{
    if( !isStateSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "translational state" );
    }
    return currentState_;
}

void Body::setState( const Eigen::Vector6d& state )
{
    currentState_ = state;
    isStateSet_ = true;
}

void Body::setLongState( const Eigen::Matrix< long double, 6, 1 >& longState )
{
    currentLongState_ = longState;
    currentState_ = longState.cast< double >( );
    isStateSet_ = true;
}

Eigen::Vector7d Body::getRotationalStateVector( )
{
    Eigen::Vector7d rotationalStateVector;
    rotationalStateVector.segment( 0, 4 ) =
            linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( currentRotationToGlobalFrame_ ) );
    rotationalStateVector.segment( 4, 3 ) = currentAngularVelocityVectorInLocalFrame_;
    return rotationalStateVector;
}

Eigen::Vector3d Body::getPosition( )
{
    if( !isStateSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "translational state (position only)" );
    }
    return currentState_.segment( 0, 3 );
}

Eigen::Vector3d Body::getVelocity( )
{
    if( !isStateSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "translational state (velocity only)" );
    }
    return currentState_.segment( 3, 3 );
}

Eigen::Matrix< long double, 6, 1 > Body::getLongState( )
{
    if( !isStateSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "translational state" );
    }
    return currentLongState_;
}

Eigen::Matrix< long double, 3, 1 > Body::getLongPosition( )
{
    if( !isStateSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "translational state (position only)" );
    }
    return currentLongState_.segment( 0, 3 );
}

Eigen::Matrix< long double, 3, 1 > Body::getLongVelocity( )
{
    if( !isStateSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "translational state (velocity only)" );
    }
    return currentLongState_.segment( 3, 3 );
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

Eigen::Quaterniond Body::getCurrentRotationToGlobalFrame( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (rotation body-fixed to global frame)" );
    }
    return currentRotationToGlobalFrame_;
}

Eigen::Quaterniond& Body::getCurrentRotationToGlobalFrameReference( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (rotation body-fixed to global frame)" );
    }
    return currentRotationToGlobalFrame_;
}

Eigen::Quaterniond Body::getCurrentRotationToLocalFrame( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (rotation global to body-fixed frame)" );
    }
    return currentRotationToLocalFrame_;
}

Eigen::Matrix3d Body::getCurrentRotationMatrixToGlobalFrame( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (rotation body-fixed to global frame)" );
    }
    return Eigen::Matrix3d( currentRotationToLocalFrame_.inverse( ) );
}

Eigen::Matrix3d Body::getCurrentRotationMatrixToLocalFrame( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (rotation global to body-fixed frame)" );
    }
    return Eigen::Matrix3d( currentRotationToLocalFrame_ );
}

Eigen::Vector7d Body::getCurrentRotationalState( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (rotation quaternion and angular velocity)" );
    }
    return ( Eigen::VectorXd( 7 ) << linear_algebra::convertQuaternionToVectorFormat( getCurrentRotationToGlobalFrame( ) ),
             getCurrentAngularVelocityVectorInGlobalFrame( ) )
            .finished( );
}

Eigen::Matrix3d Body::getCurrentRotationMatrixDerivativeToGlobalFrame( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (rotation derivative body-fixed to global frame)" );
    }
    if( currentRotationToLocalFrameDerivative_.hasNaN( ) )
    {
        throw std::runtime_error( "Error when retrieving derivative of rotation to global frame from body " + bodyName_ +
                                  ", matrix is undefined" );
    }
    return currentRotationToLocalFrameDerivative_.transpose( );
}

Eigen::Matrix3d Body::getCurrentRotationMatrixDerivativeToLocalFrame( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (rotation derivative global to body-fixed  frame)" );
    }
    if( currentRotationToLocalFrameDerivative_.hasNaN( ) )
    {
        throw std::runtime_error( "Error when retrieving derivative of rotation to local frame from body " + bodyName_ +
                                  ", matrix is undefined" );
    }
    return currentRotationToLocalFrameDerivative_;
}

Eigen::Vector3d Body::getCurrentAngularVelocityVectorInGlobalFrame( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (angular velocity vector in global frame)" );
    }
    return currentAngularVelocityVectorInGlobalFrame_;
}

Eigen::Vector3d Body::getCurrentAngularVelocityVectorInLocalFrame( )
{
    if( !isRotationSet_ )
    {
        throw exceptions::BodyDuringPropagationError( bodyName_, "rotational state (angular velocity vector in body-fixed frame)" );
    }
    return currentAngularVelocityVectorInLocalFrame_;
}

void Body::setEphemeris( const std::shared_ptr< ephemerides::Ephemeris > bodyEphemeris )
{
    bodyEphemeris_ = bodyEphemeris;
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

void Body::setAtmosphereModel( const std::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel )
{
    atmosphereModel_ = atmosphereModel;
}

void Body::setRotationalEphemeris( const std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris )
{
    rotationalEphemeris_ = rotationalEphemeris;
}

void Body::setShapeModel( const std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel )
{
    shapeModel_ = shapeModel;
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

std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > > Body::getBodyDeformationModels( )
{
    return bodyDeformationModels_;
}

std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >& Body::getBodyDeformationModelsReference( )
{
    return bodyDeformationModels_;
}

void Body::addBodyDeformationModel( const std::shared_ptr< basic_astrodynamics::BodyDeformationModel > deformationModel )
{
    bodyDeformationModels_.push_back( deformationModel );
}

void Body::setRadiationPressureInterface( const std::string& radiatingBody,
                                          const std::shared_ptr< electromagnetism::RadiationPressureInterface > radiationPressureInterface )
{
    radiationPressureInterfaces_[ radiatingBody ] = radiationPressureInterface;
}

void Body::setRadiationSourceModel( const std::shared_ptr< electromagnetism::RadiationSourceModel > radiationSourceModel )
{
    radiationSourceModel_ = radiationSourceModel;
}

void Body::setRadiationPressureTargetModels(
        const std::vector< std::shared_ptr< electromagnetism::RadiationPressureTargetModel > > radiationPressureTargetModel )
{
    radiationPressureTargetModels_ = radiationPressureTargetModel;
}

void Body::addRadiationPressureTargetModel(
        const std::shared_ptr< electromagnetism::RadiationPressureTargetModel > radiationPressureTargetModel )
{
    radiationPressureTargetModels_.push_back( radiationPressureTargetModel );
}

void Body::setGravityFieldVariationSet( const std::shared_ptr< gravitation::GravityFieldVariationsSet > gravityFieldVariationSet )
{
    gravityFieldVariationSet_ = gravityFieldVariationSet;
}

std::shared_ptr< gravitation::GravityFieldModel > Body::getGravityFieldModel( )
{
    return gravityFieldModel_;
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

std::shared_ptr< ephemerides::Ephemeris > Body::getEphemeris( )
{
    return bodyEphemeris_;
}

std::shared_ptr< aerodynamics::AtmosphereModel > Body::getAtmosphereModel( )
{
    return atmosphereModel_;
}

std::shared_ptr< ephemerides::RotationalEphemeris > Body::getRotationalEphemeris( )
{
    return rotationalEphemeris_;
}

std::shared_ptr< basic_astrodynamics::BodyShapeModel > Body::getShapeModel( )
{
    return shapeModel_;
}

std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > Body::getAerodynamicCoefficientInterface( )
{
    return aerodynamicCoefficientInterface_;
}

std::shared_ptr< aerodynamics::FlightConditions > Body::getFlightConditions( )
{
    return aerodynamicFlightConditions_;
}

std::map< std::string, std::shared_ptr< electromagnetism::RadiationPressureInterface > > Body::getRadiationPressureInterfaces( )
{
    return radiationPressureInterfaces_;
}

const std::shared_ptr< electromagnetism::RadiationSourceModel > Body::getRadiationSourceModel( ) const
{
    return radiationSourceModel_;
}

const std::vector< std::shared_ptr< electromagnetism::RadiationPressureTargetModel > > Body::getRadiationPressureTargetModels( ) const
{
    return radiationPressureTargetModels_;
}

const std::shared_ptr< electromagnetism::RadiationPressureTargetModel > Body::getRadiationPressureTargetModel( ) const
{
    if( radiationPressureTargetModels_.size( ) == 1 )
    {
        return radiationPressureTargetModels_.at( 0 );
    }
    if( radiationPressureTargetModels_.size( ) == 0 )
    {
        return nullptr;
    }
    throw std::runtime_error( "Error, could not unambiguously retrieve radiation pressure target model, found " +
                              std::to_string( radiationPressureTargetModels_.size( ) ) + " models." );
}

std::pair< bool, std::shared_ptr< gravitation::GravityFieldVariations > > Body::getGravityFieldVariation(
        const gravitation::BodyDeformationTypes& deformationType,
        const std::string identifier )
{
    return gravityFieldVariationSet_->getGravityFieldVariation( deformationType, identifier );
}

std::shared_ptr< gravitation::GravityFieldVariationsSet > Body::getGravityFieldVariationSet( )
{
    return gravityFieldVariationSet_;
}

std::shared_ptr< system_models::VehicleSystems > Body::getVehicleSystems( )
{
    if( vehicleSystems_ == nullptr )
    {
        vehicleSystems_ = std::make_shared< system_models::VehicleSystems >( );
    }
    return vehicleSystems_;
}

void Body::setVehicleSystems( const std::shared_ptr< system_models::VehicleSystems > vehicleSystems )
{
    vehicleSystems_ = vehicleSystems;
}

std::shared_ptr< RigidBodyProperties > Body::getMassProperties( )
{
    return massProperties_;
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

Eigen::Vector3d Body::getBodyFixedCenterOfMass( )
{
    return massProperties_->getCurrentCenterOfMass( );
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

void Body::addGroundStation( const std::string& stationName, const std::shared_ptr< ground_stations::GroundStation >& station )
{
    groundStationMap[ stationName ] = station;
}

std::shared_ptr< ground_stations::GroundStation > Body::getGroundStation( const std::string& stationName ) const
{
    if( groundStationMap.count( stationName ) == 0 )
    {
        throw std::runtime_error( "Error, station " + stationName + " does not exist" );
    }
    return groundStationMap.at( stationName );
}

std::map< std::string, std::shared_ptr< ground_stations::GroundStation > > Body::getGroundStationMap( ) const
{
    return groundStationMap;
}

void Body::recomputeStateOnNextCall( )
{
    timeOfCurrentState_ = Time( TUDAT_NAN );
}

double Body::getDoubleTimeOfCurrentState( )
{
    return static_cast< double >( timeOfCurrentState_ );
}

int Body::getIsBodyGlobalFrameOrigin( )
{
    return bodyIsGlobalFrameOrigin_;
}

void Body::setIsBodyGlobalFrameOrigin( const int bodyIsGlobalFrameOrigin )
{
    bodyIsGlobalFrameOrigin_ = bodyIsGlobalFrameOrigin;
}

void Body::getPositionByReference( Eigen::Vector3d& position )
{
    position = currentState_.segment( 0, 3 );
}

std::string Body::getBodyName( )
{
    return bodyName_;
}

void Body::setBodyName( const std::string bodyName )
{
    bodyName_ = bodyName;
}

void Body::setIonosphereModel( const std::shared_ptr< environment::IonosphereModel >& ionosphereModel )
{
    ionosphereModel_ = ionosphereModel;
}

std::shared_ptr< environment::IonosphereModel > Body::getIonosphereModel( ) const
{
    return ionosphereModel_;
}

// template void Body::setStateFromEphemeris< double, double >( const double& time );

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

    for( auto bodyIterator : bodies.getMap( ) )
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
    for( auto bodyIterator : bodies )
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
    for( auto bodyIterator : bodies.getMap( ) )
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
