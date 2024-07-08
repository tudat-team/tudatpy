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
PeriodicGravityFieldVariations::PeriodicGravityFieldVariations(
        const std::vector< Eigen::MatrixXd >& cosineAmplitudes,
        const std::vector< Eigen::MatrixXd >& sineAmplitudes,
        const std::vector< double >& frequencies,
        const std::vector< double >& phases,
        const double referenceEpoch,
        const int minimumDegree,
        const int minimumOrder ):
    GravityFieldVariations( minimumDegree, minimumOrder, minimumDegree + cosineAmplitudes.at( 0 ).rows( ) - 1,
                            minimumOrder + cosineAmplitudes.at( 0 ).cols( ) - 1 ),
    cosineAmplitudes_( cosineAmplitudes ),
    sineAmplitudes_( sineAmplitudes ),
    frequencies_( frequencies ),
    phases_( phases ),
    referenceEpoch_( referenceEpoch )
{
    if( cosineAmplitudes.size( ) != sineAmplitudes.size( ) )
    {
        throw std::runtime_error( "Error when making periodic gravity field variations, amplitude input sizes inconsistent" );
    }

    if( cosineAmplitudes.size( ) != frequencies.size( ) )
    {
        throw std::runtime_error( "Error when making periodic gravity field variations, frequency input size inconsistent" );
    }

    if( cosineAmplitudes.size( ) != phases.size( ) )
    {
        throw std::runtime_error( "Error when making periodic gravity field variations, phase input size inconsistent" );
    }
}


std::pair< Eigen::MatrixXd, Eigen::MatrixXd > PeriodicGravityFieldVariations::calculateSphericalHarmonicsCorrections(
        const double time )
{
    Eigen::MatrixXd cosineCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sineCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );

    for( unsigned int i = 0; i < frequencies_.size( ); i++ )
    {
        cosineCorrections += cosineAmplitudes_.at( i ) * std::cos(
                    frequencies_.at( i ) * ( time - referenceEpoch_ ) + phases_.at( i ) );
        sineCorrections += sineAmplitudes_.at( i ) * std::cos(
                    frequencies_.at( i ) * ( time - referenceEpoch_ ) + phases_.at( i ) );
    }

    return std::make_pair( cosineCorrections, sineCorrections );
}

} // namespace gravitation

} // namespace tudat
