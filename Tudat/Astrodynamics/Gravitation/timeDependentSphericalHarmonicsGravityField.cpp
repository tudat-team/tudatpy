/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"

namespace tudat
{

namespace gravitation
{

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
}

} // namespace gravitation

} // namespace tudat
