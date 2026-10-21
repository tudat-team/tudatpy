/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/rigidBodyProperties.h"

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"

namespace tudat
{

namespace simulation_setup
{

RigidBodyProperties::RigidBodyProperties( ):
    currentMass_( TUDAT_NAN ), currentCenterOfMass_( Eigen::Vector3d::Constant( TUDAT_NAN ) ),
    currentInertiaTensor_( Eigen::Matrix3d::Constant( TUDAT_NAN ) ), isBodyInPropagation_( false ), isMassComputed_( false ),
    isComComputed_( false ), isInertiaTensorComputed_( false )
{}

RigidBodyProperties::~RigidBodyProperties( ) {}

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
{}

TimeDependentRigidBodyProperties::TimeDependentRigidBodyProperties( const double mass,
                                                                    const Eigen::Vector3d& centerOfMass,
                                                                    const Eigen::Matrix3d& inertiaTensor ):
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

TimeDependentRigidBodyProperties::~TimeDependentRigidBodyProperties( ) {}

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

MassDependentRigidBodyProperties::~MassDependentRigidBodyProperties( ) {}

void MassDependentRigidBodyProperties::updateMass( const double currentTime ) {}

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
        const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel ): gravityFieldModel_( gravityFieldModel )
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

FromGravityFieldRigidBodyProperties::~FromGravityFieldRigidBodyProperties( ) {}

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

}  // namespace simulation_setup

}  // namespace tudat
