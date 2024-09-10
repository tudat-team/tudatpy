/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/estimatable_parameters/tidalLoveNumber.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Get value of Love number k_{n}
Eigen::VectorXd FullDegreeTidalLoveNumber::getParameterValue( )
{
    // Retrieve complex Love numbers at required degree
    std::vector< std::complex< double > > fullLoveNumbers = gravityFieldVariationModel_->getLoveNumbersOfDegree( degree_ );

    // Compute mean value across orders
    Eigen::VectorXd meanLoveNumber = Eigen::VectorXd::Zero( parameterSize_ );
    for( int i = 0; i <= degree_; i ++ )
    {
        meanLoveNumber[ 0 ] += fullLoveNumbers[ i ].real( );
        if( useComplexComponents_ )
        {
            meanLoveNumber[ 1 ] += fullLoveNumbers[ i ].imag( );
        }
    }

    // Return mean Love number across orders
    meanLoveNumber = meanLoveNumber / static_cast< double >( degree_ + 1 );
    return meanLoveNumber;
}

//! Reset value of Love number k_{n}
void FullDegreeTidalLoveNumber::setParameterValue( Eigen::VectorXd parameterValue )
{
    // Retrieve complex Love numbers at required degree
    std::vector< std::complex< double > > fullLoveNumbers;
    fullLoveNumbers.resize( degree_ + 1 );

    // Set complex value of Love numbers
    double complexPart = 0.0;
    if( useComplexComponents_ )
    {
        complexPart = parameterValue[ 1 ];
    }
    else
    {
        std::vector< std::complex< double > > currentLoveNumbers = gravityFieldVariationModel_->getLoveNumbersOfDegree(
                    degree_ );
        double meanComplexNumber = 0.0;
        for( int i = 0; i <= degree_; i ++ )
        {
            meanComplexNumber += currentLoveNumbers[ i ].imag( );
        }
        meanComplexNumber /= ( degree_ + 1.0 );
        complexPart = meanComplexNumber;
    }

    // Modify required values of Love numbers
    std::complex< double > complexLoveNumber = std::complex< double >( parameterValue[ 0 ], complexPart );
    for( int i = 0; i <= degree_; i ++ )
    {
        fullLoveNumbers[ i ] = complexLoveNumber;
    }

    // Reset Love numbers
    gravityFieldVariationModel_->resetLoveNumbersOfDegree( fullLoveNumbers, degree_ );
}

//! Get value of Love number k_{n,m}
Eigen::VectorXd SingleDegreeVariableTidalLoveNumber::getParameterValue( )
{
    // Retrieve complex Love numbers at required degree
    std::vector< std::complex< double > > loveNumbers = gravityFieldVariationModel_->getLoveNumbersOfDegree( degree_ );

    // Retrieve required values
    Eigen::VectorXd loveNumberVector = Eigen::VectorXd::Zero( parameterSize_ );
    for( unsigned int i = 0; i < orders_.size( ); i++ )
    {
        if( useComplexComponents_ )
        {
            loveNumberVector[ 2 * i ] = loveNumbers[ orders_[ i ] ].real( );
            loveNumberVector[ 2 * i + 1 ] = loveNumbers[ orders_[ i ] ].imag( );
        }
        else
        {
            loveNumberVector[ i ] = loveNumbers[ orders_[ i ] ].real( );
        }
    }

    return loveNumberVector;
}

//! Reset value of Love number k_{n,m}
void SingleDegreeVariableTidalLoveNumber::setParameterValue( Eigen::VectorXd parameterValue )
{
    // Retrieve current complex Love numbers at required degree
    std::vector< std::complex< double > > fullLoveNumbers =
            gravityFieldVariationModel_->getLoveNumbersOfDegree( degree_ );

    // Modify required values
    for( unsigned int i = 0; i < orders_.size( ); i++ )
    {
        if( useComplexComponents_ )
        {
            fullLoveNumbers[ orders_[ i ] ] = std::complex< double >( parameterValue[ 2 * i ], parameterValue[ 2 * i + 1 ] );
        }
        else
        {
            fullLoveNumbers[ orders_[ i ] ] = std::complex< double >(
                        parameterValue[ i ], fullLoveNumbers[ orders_[ i ] ].imag( ) );
        }
    }

    // Reset values
    gravityFieldVariationModel_->resetLoveNumbersOfDegree( fullLoveNumbers, degree_ );
}

ModeCoupledTidalLoveNumber::ModeCoupledTidalLoveNumber(
    const std::shared_ptr< gravitation::ModeCoupledSolidBodyTideGravityFieldVariations > gravityFieldVariationModel,
    const std::string& associatedBody,
    const std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > loveNumberIndices,
    const bool useComplexComponents ):
    EstimatableParameter< Eigen::VectorXd >( mode_coupled_tidal_love_numbers, associatedBody ),
    gravityFieldVariationModel_( gravityFieldVariationModel ),
    loveNumberIndices_( loveNumberIndices )
{
    if( useComplexComponents )
    {
        throw std::runtime_error( "Error, complex mode-coupled Love numbers not yet supported" );
    }
    std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > loveNumbers = gravityFieldVariationModel->getLoveNumbers( );

    parameterSize_ = 0.0;
    maximumForcingDegree_ = 0;

    int currentResponseIndex = 0;

    std::map< int, std::map< int, int > > orderIndexPerDegree;
    for( auto it: loveNumberIndices )
    {
        if( it.first.first > maximumForcingDegree_ )
        {
            maximumForcingDegree_ = it.first.first;
        }
        if( loveNumbers.count( it.first ) == 0 )
        {
            throw std::runtime_error( "Error when estimating mode-coupled Love number, no number at forcing D/O " +
                                      std::to_string( it.first.first ) + "/" + std::to_string( it.first.second ) + " found ");
        }
        for( unsigned int i = 0; i < it.second.size( ); i++ )
        {
            if( loveNumbers.at( it.first ).count( it.second.at( i ) ) == 0 )
            {
                throw std::runtime_error( "Error when estimating mode-coupled Love number, no number at forcing D/O " +
                                          std::to_string( it.first.first ) + "/" + std::to_string( it.first.second ) +
                                          " and response D/O " +  std::to_string( it.second.at( i ).first ) + "/" + std::to_string( it.second.at( i ).second ) +" found ");
            }
            std::pair< int, int > currentForcingDegreeOrder = std::make_pair(
                it.second.at( i ).first, it.second.at( i ).second );
            if( std::find(responseDegreeOrders_.begin(), responseDegreeOrders_.end(), currentForcingDegreeOrder) == responseDegreeOrders_.end( ) )
            {
                responseDegreeOrders_.push_back( currentForcingDegreeOrder );
                currentResponseIndex = responseDegreeOrders_.size( ) - 1;
            }
            else
            {
                auto findIterator = std::find(responseDegreeOrders_.begin(), responseDegreeOrders_.end(), currentForcingDegreeOrder );
                currentResponseIndex = std::distance(responseDegreeOrders_.begin(), findIterator);
            }
            responseIndices_.push_back( currentResponseIndex );
        }

        int forcingDegree = it.first.first;
        int forcingOrder = it.first.second;
        if( forcingOrdersPerDegree_.count( forcingDegree ) == 0 )
        {
            forcingOrdersPerDegree_[ forcingDegree ].push_back( forcingOrder );
        }
        else
        {
            std::vector<int> ordersInCurrentDegree = forcingOrdersPerDegree_.at( forcingDegree );
            if(std::find(ordersInCurrentDegree.begin(), ordersInCurrentDegree.end(), forcingOrder) == ordersInCurrentDegree.end( ) )
            {
                ordersInCurrentDegree.push_back( forcingOrder );
                forcingOrdersPerDegree_[ forcingDegree ] = ordersInCurrentDegree;
            }
        }
        std::vector<int> ordersInCurrentDegree = forcingOrdersPerDegree_.at( forcingDegree );

        auto findIterator = std::find(ordersInCurrentDegree.begin(), ordersInCurrentDegree.end(), forcingOrder );
        int index = std::distance(ordersInCurrentDegree.begin(), findIterator);
        for( unsigned int i = 0; i < it.second.size( ); i++ )
        {
            parameterForcingDegreeAndOrderIndices_.push_back( std::make_pair( forcingDegree, index ));
        }
        parameterSize_ += it.second.size( );
    }

    if( useComplexComponents )
    {
        parameterSize_ *= 2;
    }
}



}

}


