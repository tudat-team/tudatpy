#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"

namespace tudat
{

namespace gravitation
{

//! Function to (re)set the gravity field variations
void TimeDependentSphericalHarmonicsGravityField::setFieldVariationSettings(
        const boost::shared_ptr< GravityFieldVariationsSet > gravityFieldVariationUpdateSettings,
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
    gravityFieldVariationsSet_ = boost::shared_ptr< GravityFieldVariationsSet >( );
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

}

}
