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
    currentInertiaTensor_( Eigen::Matrix3d::Constant( TUDAT_NAN ) ), currentDerivativeInertiaTensor_( Eigen::Matrix3d::Zero( ) ),
    isBodyInPropagation_( false ), isMassComputed_( false ), isComComputed_( false ), isInertiaTensorComputed_( false ),
    isDerivativeInertiaTensorComputed_( false ), isInertiaTensorAvailable_( false ), isDerivativeInertiaTensorAvailable_( false )
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
    isDerivativeInertiaTensorComputed_ = false;
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
    if( !isInertiaTensorAvailable_ || !isInertiaTensorComputed_ )
    {
        throw std::runtime_error( "Error when retrieving inertia tensor, inertia tensor is not computed/defined." );
    }
    return currentInertiaTensor_;
}

Eigen::Matrix3d RigidBodyProperties::getCurrentDerivativeInertiaTensor( )
{
    if( !isDerivativeInertiaTensorAvailable_ || !isDerivativeInertiaTensorComputed_ )
    {
        throw std::runtime_error( "Error when retrieving derivative of the inertia tensor, not computed/defined." );
    }
    return currentDerivativeInertiaTensor_;
}

bool RigidBodyProperties::isInertiaTensorAvailable( ) const
{
    return isInertiaTensorAvailable_;
}

bool RigidBodyProperties::isInertiaTensorDerivativeAvailable( ) const
{
    return isDerivativeInertiaTensorAvailable_;
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
{
    isInertiaTensorAvailable_ = ( inertiaTensorFunction_ != nullptr );
}

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
        isInertiaTensorAvailable_ = true;
        isDerivativeInertiaTensorAvailable_ = true;
        isDerivativeInertiaTensorComputed_ = true;
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
    isInertiaTensorAvailable_ = true;
    currentDerivativeInertiaTensor_.setZero( );
    isDerivativeInertiaTensorAvailable_ = true;
    isDerivativeInertiaTensorComputed_ = true;
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
    isInertiaTensorAvailable_ = ( inertiaTensorFunction_ != nullptr );
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
    if( gravityFieldModel_ == nullptr )
    {
        throw std::runtime_error( "Error when creating gravity-linked rigid body properties: gravity field is null." );
    }

    modelIsTimeDependent_ =
            ( std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel_ ) != nullptr );

    currentMass_ = gravityFieldModel_->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT;
    currentCenterOfMass_ = gravityFieldModel_->getCenterOfMass( );
    currentDerivativeInertiaTensor_ = Eigen::Matrix3d::Zero( );

    isMassComputed_ = true;
    isComComputed_ = true;
    isInertiaTensorAvailable_ = gravityFieldModel_->hasInertiaTensor( );
    // A static gravity-derived tensor has an identically zero derivative. For a
    // time-dependent field, availability is enabled when a coefficient-rate provider
    // (currently an integrated gravity-field variation) supplies its first value.
    isDerivativeInertiaTensorAvailable_ = isInertiaTensorAvailable_ && !modelIsTimeDependent_;
    if( isInertiaTensorAvailable_ )
    {
        currentInertiaTensor_ = gravityFieldModel_->getInertiaTensor( );
        isInertiaTensorComputed_ = true;
        isDerivativeInertiaTensorComputed_ = !modelIsTimeDependent_;
    }
}

FromGravityFieldRigidBodyProperties::~FromGravityFieldRigidBodyProperties( ) {}

void FromGravityFieldRigidBodyProperties::resetCurrentTime( )
{
    if( modelIsTimeDependent_ )
    {
        isMassComputed_ = false;
        isComComputed_ = false;
        isInertiaTensorComputed_ = false;
        isDerivativeInertiaTensorComputed_ = false;
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
    if( ( modelIsTimeDependent_ && !isComComputed_ ) || !isBodyInPropagation_ )
    {
        currentCenterOfMass_ = gravityFieldModel_->getCenterOfMass( );
        isComComputed_ = true;
    }
    if( isInertiaTensorAvailable_ && ( ( modelIsTimeDependent_ && !isInertiaTensorComputed_ ) || !isBodyInPropagation_ ) )
    {
        currentInertiaTensor_ = gravityFieldModel_->getInertiaTensor( );
        isInertiaTensorComputed_ = true;
    }
}

void FromGravityFieldRigidBodyProperties::updateInertiaTensorDerivative( const Eigen::Vector5d& derivativeDegreeTwoCoefficients )
{
    const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( gravityFieldModel_ );
    if( !isInertiaTensorAvailable_ || sphericalHarmonicsGravityField == nullptr )
    {
        throw std::runtime_error(
                "Error when updating inertia tensor derivative: gravity field does not provide degree-two inertia data." );
    }

    currentDerivativeInertiaTensor_ = gravitation::computeDerivativeInertiaTensor( derivativeDegreeTwoCoefficients[ 0 ],
                                                                                   derivativeDegreeTwoCoefficients[ 1 ],
                                                                                   derivativeDegreeTwoCoefficients[ 2 ],
                                                                                   derivativeDegreeTwoCoefficients[ 3 ],
                                                                                   derivativeDegreeTwoCoefficients[ 4 ],
                                                                                   currentMass_,
                                                                                   sphericalHarmonicsGravityField->getReferenceRadius( ) );
    isDerivativeInertiaTensorAvailable_ = true;
    isDerivativeInertiaTensorComputed_ = true;
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
    isMassComputed_ = true;
    isComComputed_ = true;
    if( isInertiaTensorAvailable_ )
    {
        currentInertiaTensor_ = gravityFieldModel_->getInertiaTensor( );
        isInertiaTensorComputed_ = true;
        if( modelIsTimeDependent_ )
        {
            isDerivativeInertiaTensorComputed_ = false;
        }
    }

    isBodyInPropagation_ = isBodyInPropagation;
}

}  // namespace simulation_setup

}  // namespace tudat
