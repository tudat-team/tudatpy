/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/simulation/environment_setup/rigidBodyProperties.h"

namespace tudat
{

namespace gravitation
{

void TimeDependentSphericalHarmonicsGravityField::setRigidBodyProperties(
        const std::shared_ptr< simulation_setup::RigidBodyProperties >& rigidBodyProperties )
{
    GravityFieldModel::setRigidBodyProperties( rigidBodyProperties );
    gravityDerivedRigidBodyProperties_ =
            std::dynamic_pointer_cast< simulation_setup::FromGravityFieldRigidBodyProperties >( rigidBodyProperties );
}

void TimeDependentSphericalHarmonicsGravityField::updateInertiaTensorDerivative( const double time )
{
    const std::shared_ptr< simulation_setup::FromGravityFieldRigidBodyProperties > rigidBodyProperties =
            gravityDerivedRigidBodyProperties_.lock( );
    if( rigidBodyProperties == nullptr || !rigidBodyProperties->isInertiaTensorAvailable( ) )
    {
        return;
    }

    Eigen::MatrixXd sineCoefficientDerivatives = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd cosineCoefficientDerivatives = Eigen::MatrixXd::Zero( 3, 3 );
    if( gravityFieldVariationsSet_ != nullptr )
    {
        gravityFieldVariationsSet_->addSphericalHarmonicsCorrectionTimeDerivatives(
                time, sineCoefficientDerivatives, cosineCoefficientDerivatives );
    }

    Eigen::Vector5d degreeTwoCoefficientDerivatives;
    degreeTwoCoefficientDerivatives << cosineCoefficientDerivatives( 2, 0 ), cosineCoefficientDerivatives( 2, 1 ),
            cosineCoefficientDerivatives( 2, 2 ), sineCoefficientDerivatives( 2, 1 ), sineCoefficientDerivatives( 2, 2 );
    rigidBodyProperties->updateInertiaTensorDerivative( degreeTwoCoefficientDerivatives );
}

//! Function to (re)set the gravity field variations
void TimeDependentSphericalHarmonicsGravityField::setFieldVariationSettings(
        const std::shared_ptr< GravityFieldVariationsSet > gravityFieldVariationUpdateSettings,
        const bool updateCorrections )
{
    // Set new variation set.
    gravityFieldVariationsSet_ = gravityFieldVariationUpdateSettings;

    // Update correction functions if necessary.
    if( updateCorrections )
    {
        updateCorrectionFunctions( );
    }
}

//! Function to clear all gravity field variations
void TimeDependentSphericalHarmonicsGravityField::clearVariations( )
{
    gravityFieldVariationsSet_ = std::shared_ptr< GravityFieldVariationsSet >( );
    correctionFunctions_.clear( );
}

//! Update gravity field to current time.
void TimeDependentSphericalHarmonicsGravityField::update( const double time )
{
    // Initialize current coefficients to nominal values.
    sineCoefficients_ = nominalSineCoefficients_;
    cosineCoefficients_ = nominalCosineCoefficients_;

    // Iterate over all corrections.
    for( unsigned int i = 0; i < correctionFunctions_.size( ); i++ )
    {
        // Add correction of this iteration to current coefficients.
        correctionFunctions_[ i ]( time, sineCoefficients_, cosineCoefficients_ );
    }

    notifyMassDistributionUpdate( );
    updateInertiaTensorDerivative( time );
}

}  // namespace gravitation

}  // namespace tudat
