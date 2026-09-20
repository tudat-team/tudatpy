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

#include <cmath>
#include <typeinfo>

#include "tudat/astro/basic_astro/polyhedronFuntions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"
#include "tudat/astro/gravitation/ringGravityField.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"

namespace tudat
{

namespace simulation_setup
{

namespace
{

void computeSphericalHarmonicMassDistribution(
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField,
        const double scaledMeanMomentOfInertia,
        Eigen::Vector3d& centerOfMass,
        Eigen::Matrix3d& inertiaTensor,
        bool& inertiaTensorAvailable )
{
    centerOfMass.setZero( );
    const Eigen::MatrixXd cosineCoefficients = gravityField->getCosineCoefficients( );
    const Eigen::MatrixXd sineCoefficients = gravityField->getSineCoefficients( );
    if( cosineCoefficients.rows( ) > 1 && cosineCoefficients.cols( ) > 1 && sineCoefficients.rows( ) > 1 &&
        sineCoefficients.cols( ) > 1 )
    {
        centerOfMass =
                ( Eigen::Vector3d( ) << cosineCoefficients( 1, 1 ), sineCoefficients( 1, 1 ), cosineCoefficients( 1, 0 ) ).finished( ) /
                gravityField->getReferenceRadius( ) * std::sqrt( 3.0 );
    }
    inertiaTensorAvailable = cosineCoefficients.rows( ) > 2 && cosineCoefficients.cols( ) > 2 &&
            sineCoefficients.rows( ) > 2 && sineCoefficients.cols( ) > 2 && std::isfinite( scaledMeanMomentOfInertia );
    if( inertiaTensorAvailable )
    {
        inertiaTensor = gravitation::getInertiaTensorFromGravityField( gravityField, scaledMeanMomentOfInertia );
    }
}

void computePolyhedronMassDistribution( const std::shared_ptr< gravitation::PolyhedronGravityField >& gravityField,
                                       Eigen::Vector3d& centerOfMass,
                                       Eigen::Matrix3d& inertiaTensor,
                                       bool& inertiaTensorAvailable )
{
    centerOfMass.setZero( );
    inertiaTensor =
            basic_astrodynamics::computePolyhedronInertiaTensor( gravityField->getVerticesCoordinates( ),
                                                                 gravityField->getVerticesDefiningEachFacet( ),
                                                                 gravityField->getGravitationalParameter( ),
                                                                 physical_constants::GRAVITATIONAL_CONSTANT );
    inertiaTensorAvailable = true;
}

}  // namespace

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
        const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel,
        const double scaledMeanMomentOfInertia ):
    gravityFieldModel_( nullptr ), scaledMeanMomentOfInertia_( scaledMeanMomentOfInertia ), modelIsTimeDependent_( false )
{
    resetGravityFieldModel( gravityFieldModel );
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
        synchronizeMassFromGravityField( );
    }
}

void FromGravityFieldRigidBodyProperties::updateMassDistribution( const double currentTime )
{
    if( ( modelIsTimeDependent_ && ( !isComComputed_ || !isInertiaTensorComputed_ ) ) || !isBodyInPropagation_ )
    {
        synchronizeMassDistributionFromGravityField( );
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

    currentDerivativeInertiaTensor_ = gravitation::getInertiaTensor( derivativeDegreeTwoCoefficients[ 0 ],
                                                                     derivativeDegreeTwoCoefficients[ 1 ],
                                                                     derivativeDegreeTwoCoefficients[ 2 ],
                                                                     derivativeDegreeTwoCoefficients[ 3 ],
                                                                     derivativeDegreeTwoCoefficients[ 4 ],
                                                                     0.0,
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
    synchronizeMassFromGravityField( );
    synchronizeMassDistributionFromGravityField( );
    isBodyInPropagation_ = isBodyInPropagation;
}

void FromGravityFieldRigidBodyProperties::resetGravityFieldModel(
        const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel )
{
    if( gravityFieldModel == nullptr )
    {
        throw std::runtime_error( "Error when creating gravity-linked rigid body properties: gravity field is null." );
    }

    gravityFieldModel_ = gravityFieldModel;
    modelIsTimeDependent_ =
            ( std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel_ ) != nullptr );
    currentDerivativeInertiaTensor_.setZero( );
    synchronizeMassFromGravityField( );
    synchronizeMassDistributionFromGravityField( );
}

void FromGravityFieldRigidBodyProperties::synchronizeMassFromGravityField( )
{
    currentMass_ = gravityFieldModel_->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT;
    isMassComputed_ = true;
}

void FromGravityFieldRigidBodyProperties::synchronizeMassDistributionFromGravityField( )
{
    currentCenterOfMass_.setZero( );
    isInertiaTensorAvailable_ = false;
    const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( gravityFieldModel_ );
    if( sphericalHarmonicsGravityField != nullptr )
    {
        computeSphericalHarmonicMassDistribution( sphericalHarmonicsGravityField,
                                                  scaledMeanMomentOfInertia_,
                                                  currentCenterOfMass_,
                                                  currentInertiaTensor_,
                                                  isInertiaTensorAvailable_ );
    }
    else if( const std::shared_ptr< gravitation::PolyhedronGravityField > polyhedronGravityField =
                     std::dynamic_pointer_cast< gravitation::PolyhedronGravityField >( gravityFieldModel_ ) )
    {
        computePolyhedronMassDistribution(
                polyhedronGravityField, currentCenterOfMass_, currentInertiaTensor_, isInertiaTensorAvailable_ );
    }
    else if( std::dynamic_pointer_cast< gravitation::RingGravityField >( gravityFieldModel_ ) != nullptr )
    {
        // A ring field has its center of mass at the origin; no inertia is inferred here.
    }
    else if( typeid( *gravityFieldModel_ ) == typeid( gravitation::GravityFieldModel ) )
    {
        // A central field gives a center of mass but does not determine inertia.
    }
    else
    {
        throw std::runtime_error( "Error when deriving rigid-body properties: unrecognized gravity field type." );
    }

    isComComputed_ = true;
    isInertiaTensorComputed_ = isInertiaTensorAvailable_;

    // Static gravity-derived inertia has zero derivative. A time-dependent spherical field
    // obtains derivative availability when a coefficient-rate provider supplies a value.
    isDerivativeInertiaTensorAvailable_ = isInertiaTensorAvailable_ && !modelIsTimeDependent_;
    isDerivativeInertiaTensorComputed_ = isDerivativeInertiaTensorAvailable_;
    if( isDerivativeInertiaTensorAvailable_ )
    {
        currentDerivativeInertiaTensor_.setZero( );
    }
}

double FromGravityFieldRigidBodyProperties::getScaledMeanMomentOfInertia( ) const
{
    return scaledMeanMomentOfInertia_;
}

void FromGravityFieldRigidBodyProperties::setScaledMeanMomentOfInertia( const double scaledMeanMomentOfInertia )
{
    scaledMeanMomentOfInertia_ = scaledMeanMomentOfInertia;
    synchronizeMassDistributionFromGravityField( );
}

}  // namespace simulation_setup

}  // namespace tudat
