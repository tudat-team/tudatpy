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
#include <iostream>
#include <typeinfo>

#include "tudat/astro/basic_astro/polyhedronFuntions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"
#include "tudat/astro/gravitation/ringGravityField.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/math/basic/legendrePolynomials.h"

namespace tudat
{

namespace simulation_setup
{

namespace
{

void computeSphericalHarmonicMassDistribution( const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField,
                                               const double scaledMeanMomentOfInertia,
                                               Eigen::Vector3d& centerOfMass,
                                               Eigen::Matrix3d& inertiaTensor,
                                               bool& inertiaTensorAvailable )
{
    centerOfMass.setZero( );
    const Eigen::MatrixXd& cosineCoefficients = gravityField->getCosineCoefficientsReference( );
    const Eigen::MatrixXd& sineCoefficients = gravityField->getSineCoefficientsReference( );
    const double degreeOneScale = gravityField->getReferenceRadius( ) * std::sqrt( 3.0 );
    if( cosineCoefficients.rows( ) > 1 && cosineCoefficients.cols( ) > 0 )
    {
        // For geodesy-normalized degree-one coefficients, r_CM = sqrt(3) * R * (C11, S11, C10).
        centerOfMass.z( ) = degreeOneScale * cosineCoefficients( 1, 0 );
        if( cosineCoefficients.cols( ) > 1 )
        {
            centerOfMass.x( ) = degreeOneScale * cosineCoefficients( 1, 1 );
        }
    }
    if( sineCoefficients.rows( ) > 1 && sineCoefficients.cols( ) > 1 )
    {
        centerOfMass.y( ) = degreeOneScale * sineCoefficients( 1, 1 );
    }
    inertiaTensorAvailable = cosineCoefficients.rows( ) > 2 && cosineCoefficients.cols( ) > 2 && sineCoefficients.rows( ) > 2 &&
            sineCoefficients.cols( ) > 2 && std::isfinite( scaledMeanMomentOfInertia );
    if( inertiaTensorAvailable )
    {
        inertiaTensor = gravitation::getInertiaTensorFromGravityField( gravityField, scaledMeanMomentOfInertia );
    }
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
    if( !isDerivativeInertiaTensorAvailable_ )
    {
        throw std::runtime_error( "Error when retrieving inertia tensor derivative: no derivative model is defined." );
    }
    if( !isDerivativeInertiaTensorComputed_ )
    {
        throw std::runtime_error( "Error when retrieving inertia tensor derivative: not updated for the current epoch." );
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

    if( isDerivativeInertiaTensorAvailable_ && ( !isDerivativeInertiaTensorComputed_ || !isBodyInPropagation_ ) )
    {
        currentDerivativeInertiaTensor_ = inertiaTensorDerivativeFunction_( currentTime );
        isDerivativeInertiaTensorComputed_ = true;
    }
}

void MassDependentRigidBodyProperties::setInertiaTensorDerivativeFunction(
        const std::function< Eigen::Matrix3d( const double ) > inertiaTensorDerivativeFunction )
{
    inertiaTensorDerivativeFunction_ = inertiaTensorDerivativeFunction;
    isDerivativeInertiaTensorAvailable_ = ( inertiaTensorDerivativeFunction_ != nullptr );
    isDerivativeInertiaTensorComputed_ = false;
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

    // Coefficient rates use the same geodesy normalization as the gravity coefficients.
    // The inertia conversion requires unnormalized values; these factors are constant.
    static const double normalization20 = basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    static const double normalization21 = basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    static const double normalization22 = basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
    currentDerivativeInertiaTensor_ = gravitation::getInertiaTensor( derivativeDegreeTwoCoefficients[ 0 ] * normalization20,
                                                                     derivativeDegreeTwoCoefficients[ 1 ] * normalization21,
                                                                     derivativeDegreeTwoCoefficients[ 2 ] * normalization22,
                                                                     derivativeDegreeTwoCoefficients[ 3 ] * normalization21,
                                                                     derivativeDegreeTwoCoefficients[ 4 ] * normalization22,
                                                                     0.0,
                                                                     currentMass_,
                                                                     sphericalHarmonicsGravityField->getReferenceRadius( ) );
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

    const bool gravityFieldChanged = gravityFieldModel_ != gravityFieldModel;
    gravityFieldModel_ = gravityFieldModel;
    modelIsTimeDependent_ =
            ( std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel_ ) != nullptr );
    if( gravityFieldChanged )
    {
        hasWarnedUnrecognizedGravityField_ = false;
        if( const auto polyhedronGravityField = std::dynamic_pointer_cast< gravitation::PolyhedronGravityField >( gravityFieldModel_ ) )
        {
            // The shape is fixed. Integrate over the facets only when linking a new model;
            // subsequent mass changes only scale this tensor per unit mass.
            polyhedronInertiaTensorPerUnitMass_ = basic_astrodynamics::computePolyhedronInertiaTensor(
                    polyhedronGravityField->getVerticesCoordinates( ), polyhedronGravityField->getVerticesDefiningEachFacet( ), 1.0, 1.0 );
        }
    }
    currentDerivativeInertiaTensor_.setZero( );
    isDerivativeInertiaTensorComputed_ = false;
    synchronizeMassFromGravityField( );
    synchronizeMassDistributionFromGravityField( );
    // Every gravity-derived inertia has a rate provider: zero for a static field,
    // or the variation models (including their warning/zero fallback) for a time-dependent field.
    isDerivativeInertiaTensorAvailable_ = isInertiaTensorAvailable_;
}

void FromGravityFieldRigidBodyProperties::synchronizeMassFromGravityField( )
{
    const double mass = gravityFieldModel_->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT;
    if( modelIsTimeDependent_ && mass != currentMass_ )
    {
        // Changing mu changes the inertia-rate scale; the coefficient-rate update must refresh it.
        isDerivativeInertiaTensorComputed_ = false;
    }
    currentMass_ = mass;
    isMassComputed_ = true;
}

void FromGravityFieldRigidBodyProperties::synchronizeMassDistributionFromGravityField( )
{
    currentCenterOfMass_.setZero( );
    isInertiaTensorAvailable_ = false;
    const gravitation::GravityFieldModel& gravityFieldReference = *gravityFieldModel_;
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
        currentInertiaTensor_ = currentMass_ * polyhedronInertiaTensorPerUnitMass_;
        isInertiaTensorAvailable_ = true;
    }
    else if( std::dynamic_pointer_cast< gravitation::RingGravityField >( gravityFieldModel_ ) != nullptr )
    {
        // A ring field has its center of mass at the origin; no inertia is inferred here.
    }
    else if( typeid( gravityFieldReference ) == typeid( gravitation::GravityFieldModel ) )
    {
        // A central field gives a center of mass but does not determine inertia.
    }
    else if( !hasWarnedUnrecognizedGravityField_ )
    {
        std::cerr << "Warning when deriving rigid-body properties: unrecognized gravity field type. "
                     "Mass is derived from the gravitational parameter; the center of mass defaults to the origin "
                     "and the inertia tensor is unavailable. Supply explicit rigid-body properties if needed."
                  << std::endl;
        hasWarnedUnrecognizedGravityField_ = true;
    }

    isComComputed_ = true;
    isInertiaTensorComputed_ = isInertiaTensorAvailable_;

    // A resynchronization does not remove a time-dependent model's rate provider or
    // invalidate an already computed rate. Only an epoch reset invalidates that rate.
    if( !modelIsTimeDependent_ )
    {
        currentDerivativeInertiaTensor_.setZero( );
        isDerivativeInertiaTensorComputed_ = isInertiaTensorAvailable_;
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
    isDerivativeInertiaTensorAvailable_ = isInertiaTensorAvailable_;
    if( !isDerivativeInertiaTensorAvailable_ )
    {
        isDerivativeInertiaTensorComputed_ = false;
    }
}

}  // namespace simulation_setup

}  // namespace tudat
