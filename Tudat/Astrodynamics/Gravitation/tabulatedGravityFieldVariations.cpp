/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Gravitation/tabulatedGravityFieldVariations.h"

namespace tudat
{

namespace gravitation
{

//! Class constructor
TabulatedGravityFieldVariations::TabulatedGravityFieldVariations(
        const std::map< double, Eigen::MatrixXd >& cosineCoefficientCorrections,
        const std::map< double, Eigen::MatrixXd >& sineCoefficientCorrections,
        const int minimumDegree, const int minimumOrder,
        const boost::shared_ptr< interpolators::InterpolatorSettings >interpolatorType ):
    GravityFieldVariations( minimumDegree, minimumOrder,
                            minimumDegree + cosineCoefficientCorrections.begin( )->second.rows( ) - 1,
                            minimumOrder + cosineCoefficientCorrections.begin( )->second.cols( ) - 1 ),
    interpolatorType_( interpolatorType )
{
    // Create interpolator for tabulated coefficients.
    resetCoefficientInterpolator( cosineCoefficientCorrections, sineCoefficientCorrections );
}

//! Function to (re)set the tabulated spherical harmonic coefficients.
void TabulatedGravityFieldVariations::resetCoefficientInterpolator(
        const std::map< double, Eigen::MatrixXd >& cosineCoefficientCorrections,
        const std::map< double, Eigen::MatrixXd >& sineCoefficientCorrections )
{
    // Set current coefficient tables.
    cosineCoefficientCorrections_ = cosineCoefficientCorrections;
    sineCoefficientCorrections_ = sineCoefficientCorrections;

    // Check consistency of map sizes.
    if( cosineCoefficientCorrections_.size( ) != sineCoefficientCorrections_.size( ) )
    {
        throw std::runtime_error( "Error when resetting tabulated gravity field corrections, sine and cosine data size incompatible" );
    }

    // Create iterators over maps.
    std::map< double, Eigen::MatrixXd >::iterator cosineIterator =
            cosineCoefficientCorrections_.begin( );
    std::map< double, Eigen::MatrixXd >::iterator sineIterator =
            sineCoefficientCorrections_.begin( );

    // Declare map to store concatenated (horizontally) [cosine|sine] coefficients for
    // given discrete times.
    std::map< double, Eigen::MatrixXd > sineCosinePairMap;

    // Declare matrix to set concatenated (horizontally) [cosine|sine] coefficients at
    // current time in loop.
    Eigen::MatrixXd currentCoefficients = Eigen::MatrixXd::Zero( numberOfDegrees_, 2 * numberOfOrders_ );

    // Iterate over all times in correction maps.
    for( unsigned int i = 0; i < cosineCoefficientCorrections_.size( ); i++ )
    {
        // Check whether input times are consistent
        if( cosineIterator->first != sineIterator->first )
        {
            std::string errorMessage = "Error when resetting tabulated gravity field corrections, sine and cosine data time differ by" +
                    std::to_string(  cosineIterator->first - sineIterator->first );
            throw std::runtime_error( errorMessage );
        }
        // Check whether matrix sizes are consistent.
        if( ( cosineIterator->second.rows( ) != numberOfDegrees_ ) ||
                ( sineIterator->second.rows( ) != numberOfDegrees_ ) ||
                ( cosineIterator->second.cols( ) != numberOfOrders_ ) ||
                ( sineIterator->second.cols( ) != numberOfOrders_ ) )
        {
            std::string errorMessage = "Error when resetting tabulated gravity field corrections, sine and cosine blocks of inconsistent size";
            throw std::runtime_error( errorMessage );
        }

        // Concatenate cosine and sine matrices
        currentCoefficients.block( 0, 0, numberOfDegrees_, numberOfOrders_ ) =
                cosineIterator->second;
        currentCoefficients.block( 0, numberOfOrders_, numberOfDegrees_, numberOfOrders_ ) =
                sineIterator->second;

        // Set as entry in map of concatenated matrices
        sineCosinePairMap[ cosineIterator->first ] = currentCoefficients;

        // Increment iterators.
        cosineIterator++;
        sineIterator++;
    }

    // Create interpolator
    variationInterpolator_ =
            interpolators::createOneDimensionalInterpolator< double, Eigen::MatrixXd >(
                sineCosinePairMap, interpolatorType_ );
}

//! Function for calculating corrections by interpolating tabulated corrections.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > TabulatedGravityFieldVariations::
calculateSphericalHarmonicsCorrections(
        const double time )
{
    // Interpolate corrections
    Eigen::MatrixXd cosineSinePair = variationInterpolator_->interpolate( time );

    // Split interpolated concatenated matrix and return.
    return std::make_pair( cosineSinePair.block( 0, 0, numberOfDegrees_, numberOfOrders_ ),
                           cosineSinePair.block( 0, numberOfOrders_, numberOfDegrees_, numberOfOrders_ ) );
}

} // namespace gravitation

} // namespace tudat
