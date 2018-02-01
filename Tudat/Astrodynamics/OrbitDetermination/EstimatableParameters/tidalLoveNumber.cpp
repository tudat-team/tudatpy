/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/tidalLoveNumber.h"

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

}

}


