/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/gravitation/polynomialGravityFieldVariations.h"

namespace tudat
{

namespace gravitation
{

std::pair< int, int > getMaximumDegreeOrderForPolynomialVariations(
    const std::map< int, Eigen::MatrixXd >& cosineAmplitudes,
    const std::map< int, Eigen::MatrixXd >& sineAmplitudes,
    const int minimumDegree,
    const int minimumOrder  )
{
    int maximumOrder = -1;
    int maximumDegree = -1;

    if( cosineAmplitudes.size( ) == 0 && sineAmplitudes.size( ) == 0 )
    {
        throw std::runtime_error( "Error when creating polynomial gravity field variations, no variation blocks defined." );
    }
    for( auto it : cosineAmplitudes )
    {
        if( it.second.rows( ) == 0 || it.second.cols( ) == 0 )
        {
            throw std::runtime_error( "Error when creating polynomial gravity field variations, size of matrix block must be > 0." );
        }
        if( maximumDegree < 0 )
        {
            maximumDegree = minimumDegree + it.second.rows( ) - 1;
            maximumOrder = minimumOrder + it.second.cols( ) - 1;
        }
        else if( ( maximumDegree != minimumDegree + it.second.rows( ) - 1 ) || ( maximumOrder != minimumOrder + it.second.cols( ) - 1 ) )
        {
            throw std::runtime_error( "Error when creating polynomial gravity field variations, size of each matrix block must be equal" );
        }
    }

    for( auto it : sineAmplitudes )
    {
        if( it.second.rows( ) == 0 || it.second.cols( ) == 0 )
        {
            throw std::runtime_error( "Error when creating polynomial gravity field variations, size of matrix block must be > 0." );
        }
        if( maximumDegree < 0 )
        {
            maximumDegree = minimumDegree + it.second.rows( ) - 1;
            maximumOrder = minimumOrder + it.second.cols( ) - 1;
        }
        else if( ( maximumDegree != minimumDegree + it.second.rows( ) - 1 ) || ( maximumOrder != minimumOrder + it.second.cols( ) - 1 ) )
        {
            throw std::runtime_error( "Error when creating polynomial gravity field variations, size of each matrix block must be equal" );
        }
    }

    if( maximumDegree < 0 )
    {
        throw std::runtime_error( "Error when creating polynomial gravity field variations, size of matrix blocks not defined." );
    }

    return std::make_pair( maximumDegree, maximumOrder );
}

//! Class constructor
PolynomialGravityFieldVariations::PolynomialGravityFieldVariations(
        const std::map< int, Eigen::MatrixXd >& cosineAmplitudes,
        const std::map< int, Eigen::MatrixXd >& sineAmplitudes,
        const double referenceEpoch,
        const int minimumDegree,
        const int minimumOrder ):
    GravityFieldVariations( minimumDegree, minimumOrder,
                            getMaximumDegreeOrderForPolynomialVariations( cosineAmplitudes, sineAmplitudes, minimumDegree, minimumOrder ).first,
                            getMaximumDegreeOrderForPolynomialVariations( cosineAmplitudes, sineAmplitudes, minimumDegree, minimumOrder ).second ),
    cosineAmplitudes_( cosineAmplitudes ),
    sineAmplitudes_( sineAmplitudes ),
    referenceEpoch_( referenceEpoch )
{
    for( auto it : cosineAmplitudes_ )
    {
        if( sineAmplitudes_.count( it.first ) != 0 )
        {
            if( ( it.second.rows( ) != sineAmplitudes_.at( it.first ).rows( ) ) || ( it.second.cols( ) != sineAmplitudes_.at( it.first ).cols( ) ) )
            {
                throw std::runtime_error( "Error when creating polynomial gravity field variation, sine and cosine sizes are inconsistent" );
            }
        }
    }

    for( auto it : sineAmplitudes_ )
    {
        if( cosineAmplitudes_.count( it.first ) != 0 )
        {
            if( ( it.second.rows( ) != cosineAmplitudes_.at( it.first ).rows( ) ) || ( it.second.cols( ) != cosineAmplitudes_.at( it.first ).cols( ) ) )
            {
                throw std::runtime_error( "Error when creating polynomial gravity field variation, sine and cosine sizes are inconsistent" );
            }
        }
    }
}


std::pair< Eigen::MatrixXd, Eigen::MatrixXd > PolynomialGravityFieldVariations::calculateSphericalHarmonicsCorrections(
        const double time )
{
    Eigen::MatrixXd cosineCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sineCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );

    for( auto it : cosineAmplitudes_ )
    {
        cosineCorrections += it.second * std::pow( time - referenceEpoch_, it.first );
    }

    for( auto it : sineAmplitudes_ )
    {
        sineCorrections += it.second * std::pow( time - referenceEpoch_, it.first );
    }

    return std::make_pair( cosineCorrections, sineCorrections );
}

} // namespace gravitation

} // namespace tudat
