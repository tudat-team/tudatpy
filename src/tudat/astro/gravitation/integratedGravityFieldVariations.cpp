/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#include "tudat/astro/gravitation/integratedGravityFieldVariations.h"

#include <iterator>

#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/math/interpolators/linearInterpolator.h"

namespace tudat
{

namespace gravitation
{

IntegratedGravityFieldVariations::IntegratedGravityFieldVariations( ):
    GravityFieldVariations( 2, 0, 2, 2 ), currentCoefficientCorrections_( Eigen::Vector5d::Zero( ) ),
    currentCoefficientCorrectionDerivative_( Eigen::Vector5d::Zero( ) ), isBodyInPropagation_( false )
{}

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > IntegratedGravityFieldVariations::calculateSphericalHarmonicsCorrections(
        const double time )
{
    const Eigen::Vector5d unnormalisedCorrections = getCoefficientCorrections( time );
    Eigen::MatrixXd cosineCorrections = Eigen::MatrixXd::Zero( 1, 3 );
    Eigen::MatrixXd sineCorrections = Eigen::MatrixXd::Zero( 1, 3 );

    cosineCorrections( 0, 0 ) =
            unnormalisedCorrections( 0 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    cosineCorrections( 0, 1 ) =
            unnormalisedCorrections( 1 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    cosineCorrections( 0, 2 ) =
            unnormalisedCorrections( 2 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
    sineCorrections( 0, 1 ) =
            unnormalisedCorrections( 3 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    sineCorrections( 0, 2 ) =
            unnormalisedCorrections( 4 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

    return std::make_pair( cosineCorrections, sineCorrections );
}

void IntegratedGravityFieldVariations::setCurrentCoefficientCorrections( const Eigen::VectorXd& coefficientCorrections )
{
    if( coefficientCorrections.size( ) != 5 )
    {
        throw std::runtime_error(
                "Error when setting integrated gravity-field variation: expected five degree-two coefficients." );
    }
    currentCoefficientCorrections_ = coefficientCorrections;
}

void IntegratedGravityFieldVariations::setCurrentCoefficientCorrectionDerivative(
        const Eigen::VectorXd& coefficientCorrectionDerivative )
{
    if( coefficientCorrectionDerivative.size( ) != 5 )
    {
        throw std::runtime_error(
                "Error when setting integrated gravity-field variation derivative: expected five degree-two coefficients." );
    }
    currentCoefficientCorrectionDerivative_ = coefficientCorrectionDerivative;
}

void IntegratedGravityFieldVariations::setCoefficientCorrectionHistory(
        const std::map< double, Eigen::Vector5d >& coefficientCorrectionHistory )
{
    if( coefficientCorrectionHistory.empty( ) )
    {
        throw std::runtime_error( "Error when setting integrated gravity-field variation history: input history is empty." );
    }

    coefficientCorrectionHistory_ = coefficientCorrectionHistory;
    currentCoefficientCorrections_ = coefficientCorrectionHistory_.rbegin( )->second;

    if( coefficientCorrectionHistory_.size( ) > 1 )
    {
        coefficientCorrectionInterpolator_ =
                std::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector5d > >( coefficientCorrectionHistory_ );
    }
    else
    {
        coefficientCorrectionInterpolator_.reset( );
    }
}

Eigen::Vector5d IntegratedGravityFieldVariations::getCoefficientCorrections( const double time ) const
{
    if( !isBodyInPropagation_ && coefficientCorrectionInterpolator_ != nullptr )
    {
        return coefficientCorrectionInterpolator_->interpolate( time );
    }
    if( !isBodyInPropagation_ && coefficientCorrectionHistory_.size( ) == 1 )
    {
        return coefficientCorrectionHistory_.begin( )->second;
    }
    return currentCoefficientCorrections_;
}

Eigen::Vector5d IntegratedGravityFieldVariations::getCoefficientCorrectionDerivative( const double time ) const
{
    if( isBodyInPropagation_ || coefficientCorrectionHistory_.size( ) < 2 )
    {
        return currentCoefficientCorrectionDerivative_;
    }

    auto upperEntry = coefficientCorrectionHistory_.upper_bound( time );
    if( upperEntry == coefficientCorrectionHistory_.begin( ) )
    {
        ++upperEntry;
    }
    else if( upperEntry == coefficientCorrectionHistory_.end( ) )
    {
        upperEntry = std::prev( coefficientCorrectionHistory_.end( ) );
    }
    const auto lowerEntry = std::prev( upperEntry );
    return ( upperEntry->second - lowerEntry->second ) / ( upperEntry->first - lowerEntry->first );
}

void IntegratedGravityFieldVariations::setIsBodyInPropagation( const bool isBodyInPropagation )
{
    isBodyInPropagation_ = isBodyInPropagation;
}

bool IntegratedGravityFieldVariations::getIsBodyInPropagation( ) const
{
    return isBodyInPropagation_;
}

}  // namespace gravitation

}  // namespace tudat
