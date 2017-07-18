#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/tidalLoveNumber.h"

namespace tudat
{

namespace estimatable_parameters
{

Eigen::VectorXd FullDegreeTidalLoveNumber::getParameterValue( )
{
    std::vector< std::complex< double > > fullLoveNumbers = gravityFieldVariationModel_->getLoveNumbersOfDegree( degree_ );

    Eigen::VectorXd meanLoveNumber = Eigen::VectorXd::Zero( parameterSize_ );
    for( int i = 0; i <= degree_; i ++ )
    {
        meanLoveNumber[ 0 ] += fullLoveNumbers[ i ].real( );
        if( useComplexComponents_ )
        {
            meanLoveNumber[ 1 ] += fullLoveNumbers[ i ].imag( );
        }
    }

    meanLoveNumber = meanLoveNumber / static_cast< double >( degree_ + 1 );

    return meanLoveNumber;
}


void FullDegreeTidalLoveNumber::setParameterValue( Eigen::VectorXd parameterValue )
{
    std::vector< std::complex< double > > fullLoveNumbers;
    fullLoveNumbers.resize( degree_ + 1 );

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

    std::complex< double > complexLoveNumber = std::complex< double >( parameterValue[ 0 ], complexPart );

    for( int i = 0; i <= degree_; i ++ )
    {
        fullLoveNumbers[ i ] = complexLoveNumber;
    }

    gravityFieldVariationModel_->resetLoveNumbersOfDegree( fullLoveNumbers, degree_ );
}


Eigen::VectorXd SingleDegreeVariableTidalLoveNumber::getParameterValue( )
{
    std::vector< std::complex< double > > loveNumbers = gravityFieldVariationModel_->getLoveNumbersOfDegree( degree_ );

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


void SingleDegreeVariableTidalLoveNumber::setParameterValue( Eigen::VectorXd parameterValue )
{
    std::vector< std::complex< double > > fullLoveNumbers =
            gravityFieldVariationModel_->getLoveNumbersOfDegree( degree_ );

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

    gravityFieldVariationModel_->resetLoveNumbersOfDegree( fullLoveNumbers, degree_ );
}

}

}


