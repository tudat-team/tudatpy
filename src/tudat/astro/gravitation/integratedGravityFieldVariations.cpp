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

namespace tudat
{

namespace gravitation
{

IntegratedGravityFieldVariations::IntegratedGravityFieldVariations( ):
    // The returned correction matrices cover the degree-two, order-zero-to-two block.
    GravityFieldVariations( 2, 0, 2, 2 ), currentCoefficientCorrections_( Eigen::Vector5d::Zero( ) ),
    currentCoefficientCorrectionDerivative_( Eigen::Vector5d::Zero( ) ), isBodyInPropagation_( false )
{}

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > IntegratedGravityFieldVariations::calculateSphericalHarmonicsCorrections( const double time )
{
    // Select either the live propagated state or the installed post-propagation history.
    const Eigen::Vector5d unnormalisedCorrections = getCoefficientCorrections( time );
    Eigen::MatrixXd cosineCorrections = Eigen::MatrixXd::Zero( 1, 3 );
    Eigen::MatrixXd sineCorrections = Eigen::MatrixXd::Zero( 1, 3 );

    // The numerical state is unnormalised, whereas the time-dependent gravity field consumes
    // geodesy-normalized corrections. The zero S20 entry remains untouched by construction.
    cosineCorrections( 0, 0 ) = unnormalisedCorrections( 0 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    cosineCorrections( 0, 1 ) = unnormalisedCorrections( 1 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    cosineCorrections( 0, 2 ) = unnormalisedCorrections( 2 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
    sineCorrections( 0, 1 ) = unnormalisedCorrections( 3 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    sineCorrections( 0, 2 ) = unnormalisedCorrections( 4 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

    return std::make_pair( cosineCorrections, sineCorrections );
}

void IntegratedGravityFieldVariations::setCurrentCoefficientCorrections( const Eigen::VectorXd& coefficientCorrections )
{
    // A fixed size guards the state ordering used by the propagator and inertia coupling.
    if( coefficientCorrections.size( ) != 5 )
    {
        throw std::runtime_error( "Error when setting integrated gravity-field variation: expected five degree-two coefficients." );
    }
    currentCoefficientCorrections_ = coefficientCorrections;
}

void IntegratedGravityFieldVariations::setCurrentCoefficientCorrectionDerivative( const Eigen::VectorXd& coefficientCorrectionDerivative )
{
    // The derivative follows the same unnormalised ordering as the propagated correction state.
    if( coefficientCorrectionDerivative.size( ) != 5 )
    {
        throw std::runtime_error(
                "Error when setting integrated gravity-field variation derivative: expected five degree-two coefficients." );
    }
    currentCoefficientCorrectionDerivative_ = coefficientCorrectionDerivative;
}

void IntegratedGravityFieldVariations::setCoefficientCorrectionHistory(
        const std::map< double, Eigen::Vector5d >& coefficientCorrectionHistory,
        const std::shared_ptr< interpolators::InterpolatorSettings >& interpolatorSettings )
{
    if( coefficientCorrectionHistory.empty( ) )
    {
        throw std::runtime_error( "Error when setting integrated gravity-field variation history: input history is empty." );
    }

    coefficientCorrectionHistory_ = coefficientCorrectionHistory;
    // Keep a meaningful live value as well. This is used if propagation is enabled again before
    // a new state has been injected by the environment updater.
    currentCoefficientCorrections_ = coefficientCorrectionHistory_.rbegin( )->second;

    if( coefficientCorrectionHistory_.size( ) > 1 )
    {
        // Use the result-processing interpolation policy when supplied; otherwise retain a
        // predictable linear default for direct callers.
        const std::shared_ptr< interpolators::InterpolatorSettings > settingsToUse =
                ( interpolatorSettings == nullptr ) ? interpolators::linearInterpolation( ) : interpolatorSettings;
        coefficientCorrectionInterpolator_ =
                interpolators::createOneDimensionalInterpolator< double, Eigen::Vector5d >( coefficientCorrectionHistory_, settingsToUse );
    }
    else
    {
        // A one-point history is handled explicitly as a constant in getCoefficientCorrections.
        coefficientCorrectionInterpolator_.reset( );
    }
}

Eigen::Vector5d IntegratedGravityFieldVariations::getCoefficientCorrections( const double time ) const
{
    // While propagating, the current numerical state is authoritative even when a history from a
    // previous propagation is installed. Outside propagation, the installed history is authoritative.
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
    // The derivative model supplies the live rate during propagation. A history with fewer than
    // two points cannot define a post-propagation slope, so retain the last live rate in that case.
    if( isBodyInPropagation_ || coefficientCorrectionHistory_.size( ) < 2 )
    {
        return currentCoefficientCorrectionDerivative_;
    }

    // Select the adjacent history interval containing time. At either boundary, use the first or
    // last available interval so that derivative evaluation follows the history's endpoint policy.
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
    // A secant is used deliberately: the generic state interpolator interface does not expose an
    // interpolation derivative for every supported interpolator type.
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
