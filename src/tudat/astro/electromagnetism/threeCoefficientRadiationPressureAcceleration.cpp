/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved.
 *
 *    This file is part of Tudat. Redistribution and use in source and binary
 *    forms, with or without modification, are permitted exclusively under the
 *    terms of the Modified BSD license.
 */

#include "tudat/astro/electromagnetism/threeCoefficientRadiationPressureAcceleration.h"

#include "tudat/astro/basic_astro/unitConversions.h"

namespace tudat
{
namespace electromagnetism
{

namespace
{

constexpr double earthObliquityDegrees = 23.4;

std::string getValidatedSourceName( const std::shared_ptr< IsotropicPointRadiationSourceModel >& sourceModel )
{
    if( sourceModel == nullptr )
    {
        throw std::runtime_error( "Error when creating three-coefficient radiation-pressure acceleration: source model is null." );
    }
    return sourceModel->getSourceName( );
}

void validateVectorNorm( const Eigen::Vector3d& vector, const std::string& vectorDescription )
{
    if( !( vector.norm( ) > 0.0 ) )
    {
        throw std::runtime_error( "Error when updating three-coefficient radiation-pressure acceleration: " + vectorDescription +
                                  " has zero or invalid norm." );
    }
}

}  // namespace

ThreeCoefficientRadiationPressureAcceleration::ThreeCoefficientRadiationPressureAcceleration(
        const std::shared_ptr< IsotropicPointRadiationSourceModel >& sourceModel,
        const std::shared_ptr< basic_astrodynamics::BodyShapeModel >& sourceBodyShapeModel,
        const std::function< Eigen::Vector3d( ) >& sourcePositionFunction,
        const std::function< Eigen::Vector3d( ) >& sourceVelocityFunction,
        const std::shared_ptr< RadiationPressureTargetModel >& targetModel,
        const std::function< Eigen::Vector3d( ) >& targetPositionFunction,
        const std::function< double( ) >& targetMassFunction,
        const std::shared_ptr< OccultationModel >& sourceToTargetOccultationModel,
        const std::function< Eigen::Vector3d( ) >& referenceBodyPositionFunction,
        const std::function< Eigen::Vector3d( ) >& referenceBodyVelocityFunction,
        const Eigen::Vector3d& coefficients,
        const std::string& referenceBodyName ):
    RadiationPressureAcceleration(
            sourcePositionFunction,
            targetModel,
            targetPositionFunction,
            []( ) { return Eigen::Quaterniond::Identity( ); },
            targetMassFunction,
            sourceToTargetOccultationModel,
            getValidatedSourceName( sourceModel ) ),
    sourceModel_( sourceModel ), sourceBodyShapeModel_( sourceBodyShapeModel ), sourceVelocityFunction_( sourceVelocityFunction ),
    referenceBodyPositionFunction_( referenceBodyPositionFunction ), referenceBodyVelocityFunction_( referenceBodyVelocityFunction ),
    coefficients_( coefficients ), referenceBodyName_( referenceBodyName ), currentBasis_( Eigen::Matrix3d::Constant( TUDAT_NAN ) ),
    currentSourceToReferencePosition_( Eigen::Vector3d::Constant( TUDAT_NAN ) ),
    currentSourceToReferenceVelocity_( Eigen::Vector3d::Constant( TUDAT_NAN ) ), currentAccelerationScalingFactor_( TUDAT_NAN ),
    currentTargetMass_( TUDAT_NAN ), sourceToTargetReceivedFraction_( TUDAT_NAN )
{
    if( sourceModel_ == nullptr || sourceBodyShapeModel_ == nullptr || targetModel_ == nullptr ||
        sourceToTargetOccultationModel_ == nullptr )
    {
        throw std::runtime_error( "Error when creating three-coefficient radiation-pressure acceleration: required model is null." );
    }
    if( referenceBodyName_.empty( ) )
    {
        throw std::runtime_error( "Error when creating three-coefficient radiation-pressure acceleration: reference body is empty." );
    }
    if( !coefficients_.allFinite( ) )
    {
        throw std::runtime_error( "Error when creating three-coefficient radiation-pressure acceleration: coefficients must be finite." );
    }
}

void ThreeCoefficientRadiationPressureAcceleration::resetCoefficients( const Eigen::Vector3d& coefficients )
{
    if( !coefficients.allFinite( ) )
    {
        throw std::runtime_error( "Error when resetting three-coefficient radiation-pressure acceleration: coefficients must be finite." );
    }
    coefficients_ = coefficients;
    resetCurrentTime( );
}

void ThreeCoefficientRadiationPressureAcceleration::computeAcceleration( )
{
    sourceCenterPositionInGlobalFrame_ = sourcePositionFunction_( );
    targetCenterPositionInGlobalFrame_ = targetPositionFunction_( );
    targetCenterPositionInSourceFrame_ = targetCenterPositionInGlobalFrame_ - sourceCenterPositionInGlobalFrame_;
    validateVectorNorm( targetCenterPositionInSourceFrame_, "source-to-target position" );

    currentSourceToReferencePosition_ = referenceBodyPositionFunction_( ) - sourceCenterPositionInGlobalFrame_;
    currentSourceToReferenceVelocity_ = referenceBodyVelocityFunction_( ) - sourceVelocityFunction_( );
    validateVectorNorm( currentSourceToReferencePosition_, "source-to-reference-body position" );

    const Eigen::Vector3d orbitAngularMomentum = currentSourceToReferencePosition_.cross( currentSourceToReferenceVelocity_ );
    validateVectorNorm( orbitAngularMomentum, "reference-body source-centered orbital angular momentum" );

    const Eigen::Vector3d uDirection = -currentSourceToReferencePosition_.normalized( );
    const Eigen::Vector3d orbitNormal = orbitAngularMomentum.normalized( );
    const double obliquity = unit_conversions::convertDegreesToRadians( earthObliquityDegrees );
    const Eigen::Vector3d wDirection = std::cos( obliquity ) * orbitNormal - std::sin( obliquity ) * orbitNormal.cross( uDirection );
    const Eigen::Vector3d vDirection = wDirection.cross( uDirection );
    currentBasis_.col( 0 ) = uDirection;
    currentBasis_.col( 1 ) = vDirection;
    currentBasis_.col( 2 ) = wDirection;

    currentTargetMass_ = targetMassFunction_( );
    if( !( currentTargetMass_ > 0.0 ) )
    {
        throw std::runtime_error(
                "Error when updating three-coefficient radiation-pressure acceleration: target mass is non-positive or invalid." );
    }

    sourceToTargetReceivedFraction_ = sourceToTargetOccultationModel_->evaluateReceivedFractionFromExtendedSource(
            sourceCenterPositionInGlobalFrame_, sourceBodyShapeModel_, targetCenterPositionInGlobalFrame_ );
    receivedIrradiance = sourceModel_->evaluateIrradianceAtPosition( targetCenterPositionInSourceFrame_ ).front( ).first *
            sourceToTargetReceivedFraction_;
    currentAccelerationScalingFactor_ = receivedIrradiance / ( physical_constants::SPEED_OF_LIGHT * currentTargetMass_ );

    if( receivedIrradiance <= 0.0 )
    {
        currentUnscaledAcceleration_.setZero( );
        targetModel_->resetComputations( sourceName_ );
    }
    else
    {
        currentUnscaledAcceleration_ = currentAccelerationScalingFactor_ * currentBasis_ * coefficients_;
    }
    scaleRadiationPressureAcceleration( );
}

}  // namespace electromagnetism
}  // namespace tudat
