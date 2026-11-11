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
        const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction,
        const std::function< Eigen::Matrix3d( const double ) > inertiaTensorDerivativeFunction ):
    massFunction_( massFunction ), centerOfMassFunction_( centerOfMassFunction ), inertiaTensorFunction_( inertiaTensorFunction ),
    inertiaTensorDerivativeFunction_( inertiaTensorDerivativeFunction )
{
    isInertiaTensorAvailable_ = ( inertiaTensorFunction_ != nullptr );
    if( inertiaTensorDerivativeFunction_ != nullptr && !isInertiaTensorAvailable_ )
    {
        throw std::runtime_error(
                "Error when creating time-dependent rigid body properties: inertia tensor derivative function is "
                "defined, but inertia tensor function is not." );
    }
    isDerivativeInertiaTensorAvailable_ = ( inertiaTensorDerivativeFunction_ != nullptr );
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

    if( inertiaTensorDerivativeFunction_ != nullptr )
    {
        isDerivativeInertiaTensorComputed_ = false;
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

    // Constant inertia tensor has zero derivative (not if the inertia tensor is recomputed from a function at the next update)
    if( inertiaTensorFunction_ == nullptr )
    {
        currentDerivativeInertiaTensor_.setZero( );
        isDerivativeInertiaTensorAvailable_ = true;
        isDerivativeInertiaTensorComputed_ = true;
    }
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

    if( inertiaTensorDerivativeFunction_ != nullptr && ( !isDerivativeInertiaTensorComputed_ || !isBodyInPropagation_ ) )
    {
        currentDerivativeInertiaTensor_ = inertiaTensorDerivativeFunction_( currentTime );
        isDerivativeInertiaTensorComputed_ = true;
    }
}

MassDependentRigidBodyProperties::MassDependentRigidBodyProperties(
        const double currentMass,
        const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction,
        const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction,
        const std::function< Eigen::Matrix3d( const double ) > inertiaTensorDerivativeFunction ):
    centerOfMassFunction_( centerOfMassFunction ), inertiaTensorFunction_( inertiaTensorFunction ),
    inertiaTensorDerivativeFunction_( inertiaTensorDerivativeFunction )
{
    isInertiaTensorAvailable_ = ( inertiaTensorFunction_ != nullptr );
    if( inertiaTensorDerivativeFunction_ != nullptr && !isInertiaTensorAvailable_ )
    {
        throw std::runtime_error(
                "Error when creating mass-dependent rigid body properties: inertia tensor derivative function is "
                "defined, but inertia tensor function is not." );
    }
    setCurrentMass( currentMass );
    updateMassDistribution( TUDAT_NAN );

    // Set after initial update, since the derivative function requires a time
    isDerivativeInertiaTensorAvailable_ = ( inertiaTensorDerivativeFunction_ != nullptr );
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

void MassDependentRigidBodyProperties::setCurrentMass( const double currentMass )
{
    currentMass_ = currentMass;
    isMassComputed_ = true;
}

FromGravityFieldRigidBodyProperties::FromGravityFieldRigidBodyProperties(
        const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel,
        const double scaledMeanMomentOfInertia ):
    gravityFieldModel_( nullptr ), scaledMeanMomentOfInertia_( scaledMeanMomentOfInertia ), modelIsTimeDependent_( false ),
    massDistributionUpdateTime_( TUDAT_NAN )
{
    resetGravityFieldModel( gravityFieldModel );
}

FromGravityFieldRigidBodyProperties::~FromGravityFieldRigidBodyProperties( ) {}

void FromGravityFieldRigidBodyProperties::resetCurrentTime( )
{
    if( modelIsTimeDependent_ )
    {
        isComComputed_ = false;
        isInertiaTensorComputed_ = false;
        isDerivativeInertiaTensorComputed_ = false;
        massDistributionUpdateTime_ = TUDAT_NAN;
    }
}

void FromGravityFieldRigidBodyProperties::updateMass( const double currentTime )
{
    if( ( modelIsTimeDependent_ && !isMassComputed_ ) || !isBodyInPropagation_ )
    {
        computeMassFromGravityField( );
    }
}

void FromGravityFieldRigidBodyProperties::updateMassDistribution( const double currentTime )
{
    if( ( modelIsTimeDependent_ && ( !isComComputed_ || !isInertiaTensorComputed_ ) ) || !isBodyInPropagation_ )
    {
        computeMassDistributionFromGravityField( currentTime );
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

    // Denormalize coefficient rates, and compute inertia tensor derivative (mean moment of inertia is constant)
    Eigen::MatrixXd normalizedCosineCoefficientDerivatives = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd normalizedSineCoefficientDerivatives = Eigen::Matrix3d::Zero( );
    normalizedCosineCoefficientDerivatives( 2, 0 ) = derivativeDegreeTwoCoefficients( 0 );
    normalizedCosineCoefficientDerivatives( 2, 1 ) = derivativeDegreeTwoCoefficients( 1 );
    normalizedCosineCoefficientDerivatives( 2, 2 ) = derivativeDegreeTwoCoefficients( 2 );
    normalizedSineCoefficientDerivatives( 2, 1 ) = derivativeDegreeTwoCoefficients( 3 );
    normalizedSineCoefficientDerivatives( 2, 2 ) = derivativeDegreeTwoCoefficients( 4 );

    Eigen::MatrixXd unnormalizedCosineCoefficientDerivatives, unnormalizedSineCoefficientDerivatives;
    basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients( normalizedCosineCoefficientDerivatives,
                                                                           normalizedSineCoefficientDerivatives,
                                                                           unnormalizedCosineCoefficientDerivatives,
                                                                           unnormalizedSineCoefficientDerivatives );
    currentDerivativeInertiaTensor_ =
            gravitation::getInertiaTensor( unnormalizedCosineCoefficientDerivatives,
                                           unnormalizedSineCoefficientDerivatives,
                                           0.0,
                                           gravityFieldModel_->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT,
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
    computeMassFromGravityField( );
    computeMassDistributionFromGravityField( TUDAT_NAN );
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

    // Polyhedron inertia tensor per unit mass is recomputed for the new gravity field
    polyhedronInertiaTensorPerUnitMass_.setConstant( TUDAT_NAN );
    currentDerivativeInertiaTensor_.setZero( );
    computeMassFromGravityField( );
    computeMassDistributionFromGravityField( TUDAT_NAN );

    // Inertia tensor derivative is zero for a static field, and computed from gravity field variations for a time-dependent field
    isDerivativeInertiaTensorAvailable_ = isInertiaTensorAvailable_;
}

void FromGravityFieldRigidBodyProperties::computeMassFromGravityField( )
{
    currentMass_ = gravityFieldModel_->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT;
    isMassComputed_ = true;
}

void FromGravityFieldRigidBodyProperties::computeMassDistributionFromGravityField( const double currentTime )
{
    const gravitation::GravityFieldModel& gravityFieldReference = *gravityFieldModel_;
    if( const auto sphericalHarmonicsGravityField =
                std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( gravityFieldModel_ ) )
    {
        computeSphericalHarmonicMassDistribution( sphericalHarmonicsGravityField,
                                                  scaledMeanMomentOfInertia_,
                                                  currentCenterOfMass_,
                                                  currentInertiaTensor_,
                                                  isInertiaTensorAvailable_ );
    }
    else if( const auto polyhedronGravityField = std::dynamic_pointer_cast< gravitation::PolyhedronGravityField >( gravityFieldModel_ ) )
    {
        // The shape is fixed, so the integration over the facets is only done once per gravity field; mass changes only scale this
        // tensor per unit mass.
        if( !polyhedronInertiaTensorPerUnitMass_.allFinite( ) )
        {
            polyhedronInertiaTensorPerUnitMass_ = basic_astrodynamics::computePolyhedronInertiaTensor(
                    polyhedronGravityField->getVerticesCoordinates( ), polyhedronGravityField->getVerticesDefiningEachFacet( ), 1.0, 1.0 );
        }
        currentCenterOfMass_.setZero( );
        currentInertiaTensor_ = gravityFieldModel_->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT *
                polyhedronInertiaTensorPerUnitMass_;
        isInertiaTensorAvailable_ = true;
    }
    else if( std::dynamic_pointer_cast< gravitation::RingGravityField >( gravityFieldModel_ ) != nullptr ||
             typeid( gravityFieldReference ) == typeid( gravitation::GravityFieldModel ) )
    {
        // Point mass and ring gravity fields: center of mass at origin, no inertia tensor
        currentCenterOfMass_.setZero( );
        isInertiaTensorAvailable_ = false;
    }
    else
    {
        throw std::runtime_error( "Error when deriving rigid body properties from gravity field: gravity field type not recognized." );
    }

    // A scaled mean moment of inertia is only used for a spherical harmonic gravity field with degree two coefficients
    if( std::isfinite( scaledMeanMomentOfInertia_ ) &&
        ( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( gravityFieldModel_ ) == nullptr ||
          !isInertiaTensorAvailable_ ) )
    {
        throw std::runtime_error(
                "Error when deriving rigid body properties from gravity field: a scaled mean moment of inertia is defined, but the "
                "gravity field is not a spherical harmonic field with degree two coefficients." );
    }

    isComComputed_ = true;
    isInertiaTensorComputed_ = isInertiaTensorAvailable_;

    // Inertia tensor derivative is zero for a static field. For a time-dependent field, it is computed for the time of this update
    // when it is first requested (see getCurrentDerivativeInertiaTensor)
    massDistributionUpdateTime_ = currentTime;
    if( !modelIsTimeDependent_ )
    {
        currentDerivativeInertiaTensor_.setZero( );
        isDerivativeInertiaTensorComputed_ = isInertiaTensorAvailable_;
    }
    else
    {
        isDerivativeInertiaTensorComputed_ = false;
    }
}

Eigen::Matrix3d FromGravityFieldRigidBodyProperties::getCurrentDerivativeInertiaTensor( )
{
    // Compute inertia tensor derivative from time derivatives of degree two coefficients, once per update of the mass distribution
    if( modelIsTimeDependent_ && isDerivativeInertiaTensorAvailable_ && !isDerivativeInertiaTensorComputed_ &&
        !std::isnan( massDistributionUpdateTime_ ) )
    {
        Eigen::MatrixXd sineCoefficientDerivatives = Eigen::MatrixXd::Zero( 3, 3 );
        Eigen::MatrixXd cosineCoefficientDerivatives = Eigen::MatrixXd::Zero( 3, 3 );
        std::shared_ptr< gravitation::GravityFieldVariationsSet > gravityFieldVariationsSet =
                std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel_ )
                        ->getGravityFieldVariationsSet( );
        if( gravityFieldVariationsSet != nullptr )
        {
            gravityFieldVariationsSet->addSphericalHarmonicsCorrectionTimeDerivatives(
                    massDistributionUpdateTime_, sineCoefficientDerivatives, cosineCoefficientDerivatives );
        }

        Eigen::Vector5d degreeTwoCoefficientDerivatives;
        degreeTwoCoefficientDerivatives << cosineCoefficientDerivatives( 2, 0 ), cosineCoefficientDerivatives( 2, 1 ),
                cosineCoefficientDerivatives( 2, 2 ), sineCoefficientDerivatives( 2, 1 ), sineCoefficientDerivatives( 2, 2 );
        updateInertiaTensorDerivative( degreeTwoCoefficientDerivatives );
    }
    return RigidBodyProperties::getCurrentDerivativeInertiaTensor( );
}

double FromGravityFieldRigidBodyProperties::getScaledMeanMomentOfInertia( ) const
{
    return scaledMeanMomentOfInertia_;
}

void FromGravityFieldRigidBodyProperties::setScaledMeanMomentOfInertia( const double scaledMeanMomentOfInertia )
{
    scaledMeanMomentOfInertia_ = scaledMeanMomentOfInertia;
    computeMassDistributionFromGravityField( TUDAT_NAN );
    isDerivativeInertiaTensorAvailable_ = isInertiaTensorAvailable_;
    if( !isDerivativeInertiaTensorAvailable_ )
    {
        isDerivativeInertiaTensorComputed_ = false;
    }
}

//! Function to retrieve the gravitational parameter function of rigid body properties derived from a gravity field
std::function< double( ) > getGravitationalParameterFunctionOfGravityDerivedProperties(
        const std::shared_ptr< RigidBodyProperties > rigidBodyProperties )
{
    std::function< double( ) > gravitationalParameterFunction;
    const std::shared_ptr< FromGravityFieldRigidBodyProperties > gravityDerivedProperties =
            std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( rigidBodyProperties );
    if( gravityDerivedProperties != nullptr )
    {
        gravitationalParameterFunction =
                std::bind( &gravitation::GravityFieldModel::getGravitationalParameter, gravityDerivedProperties->getGravityFieldModel( ) );
    }
    return gravitationalParameterFunction;
}

//! Function to retrieve the inertia tensor normalization function of rigid body properties derived from a spherical harmonic
//! gravity field
std::function< double( ) > getInertiaTensorNormalizationFunctionOfGravityDerivedProperties(
        const std::shared_ptr< RigidBodyProperties > rigidBodyProperties )
{
    std::function< double( ) > inertiaTensorNormalizationFunction;
    const std::shared_ptr< FromGravityFieldRigidBodyProperties > gravityDerivedProperties =
            std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( rigidBodyProperties );
    if( gravityDerivedProperties != nullptr )
    {
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
                std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        gravityDerivedProperties->getGravityFieldModel( ) );
        if( sphericalHarmonicsGravityField != nullptr )
        {
            inertiaTensorNormalizationFunction = std::bind(
                    &gravitation::SphericalHarmonicsGravityField::getInertiaTensorNormalizationFactor, sphericalHarmonicsGravityField );
        }
    }
    return inertiaTensorNormalizationFunction;
}

}  // namespace simulation_setup

}  // namespace tudat
