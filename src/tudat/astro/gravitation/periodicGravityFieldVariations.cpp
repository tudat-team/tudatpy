/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/gravitation/periodicGravityFieldVariations.h"

namespace tudat
{

namespace gravitation
{

//! Class constructor
PeriodicGravityFieldVariations::PeriodicGravityFieldVariations( const std::vector< Eigen::MatrixXd >& cosineShAmplitudesCosineTime,
                                                                const std::vector< Eigen::MatrixXd >& cosineShAmplitudesSineTime,
                                                                const std::vector< Eigen::MatrixXd >& sineShAmplitudesCosineTime,
                                                                const std::vector< Eigen::MatrixXd >& sineShAmplitudesSineTime,
                                                                const std::vector< double >& frequencies,
                                                                const double referenceEpoch,
                                                                const int minimumDegree,
                                                                const int minimumOrder ):
    GravityFieldVariations( minimumDegree,
                            minimumOrder,
                            minimumDegree + cosineShAmplitudesCosineTime.at( 0 ).rows( ) - 1,
                            minimumOrder + cosineShAmplitudesCosineTime.at( 0 ).cols( ) - 1 ),
    cosineShAmplitudesCosineTime_( cosineShAmplitudesCosineTime ), cosineShAmplitudesSineTime_( cosineShAmplitudesSineTime ),
    sineShAmplitudesCosineTime_( sineShAmplitudesCosineTime ), sineShAmplitudesSineTime_( sineShAmplitudesSineTime ),
    frequencies_( frequencies ), referenceEpoch_( referenceEpoch )
{
    checkAmplitudes( cosineShAmplitudesCosineTime_ );
    checkAmplitudes( cosineShAmplitudesSineTime_ );
    checkAmplitudes( sineShAmplitudesCosineTime_ );
    checkAmplitudes( sineShAmplitudesSineTime_ );
}

void PeriodicGravityFieldVariations::checkAmplitudes( const std::vector< Eigen::MatrixXd >& amplitudes ) const
{
    if( amplitudes.size( ) != frequencies_.size( ) )
    {
        throw std::runtime_error( "Error configuring periodic gravity variations: amplitude and frequency counts differ." );
    }
    for( const auto& amplitude : amplitudes )
    {
        if( amplitude.rows( ) != numberOfDegrees_ || amplitude.cols( ) != numberOfOrders_ )
        {
            throw std::runtime_error( "Error configuring periodic gravity variations: amplitude block dimensions differ." );
        }
    }
}

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > PeriodicGravityFieldVariations::calculateSphericalHarmonicsCorrections( const double time )
{
    Eigen::MatrixXd cosineCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sineCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );

    for( unsigned int i = 0; i < frequencies_.size( ); i++ )
    {
        cosineCorrections += cosineShAmplitudesCosineTime_.at( i ) * std::cos( frequencies_.at( i ) * ( time - referenceEpoch_ ) ) +
                cosineShAmplitudesSineTime_.at( i ) * std::sin( frequencies_.at( i ) * ( time - referenceEpoch_ ) );
        sineCorrections += sineShAmplitudesCosineTime_.at( i ) * std::cos( frequencies_.at( i ) * ( time - referenceEpoch_ ) ) +
                sineShAmplitudesSineTime_.at( i ) * std::sin( frequencies_.at( i ) * ( time - referenceEpoch_ ) );
    }

    return std::make_pair( cosineCorrections, sineCorrections );
}

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > PeriodicGravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivative(
        const double time )
{
    Eigen::MatrixXd cosineRates = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sineRates = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    const double timeSinceEpoch = time - referenceEpoch_;

    for( unsigned int i = 0; i < frequencies_.size( ); i++ )
    {
        const double argument = frequencies_[ i ] * timeSinceEpoch;
        const double cosineTimeDerivative = -frequencies_[ i ] * std::sin( argument );
        const double sineTimeDerivative = frequencies_[ i ] * std::cos( argument );
        cosineRates += cosineShAmplitudesCosineTime_[ i ] * cosineTimeDerivative + cosineShAmplitudesSineTime_[ i ] * sineTimeDerivative;
        sineRates += sineShAmplitudesCosineTime_[ i ] * cosineTimeDerivative + sineShAmplitudesSineTime_[ i ] * sineTimeDerivative;
    }

    return std::make_pair( cosineRates, sineRates );
}

}  // namespace gravitation

}  // namespace tudat
