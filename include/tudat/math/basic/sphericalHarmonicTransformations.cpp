#include <boost/math/special_functions/binomial.hpp>
#include <iostream>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonicTransformations.h"
#include "Tudat/Mathematics/BasicMathematics/cayleyKleinParameters.h"

namespace tudat
{

namespace basic_mathematics
{

void SphericalHarmonicTransformationCache::updateFromCayleyKleinParameters(
        const std::complex< double > cayleyKleinA,
        const std::complex< double > cayleyKleinB )
{
//    std::cout<<"Updating: "<<cayleyKleinA<<" "<<cayleyKleinB<<std::endl;
    wignerDMatricesCache_->updateMatrices( cayleyKleinA, cayleyKleinB );
}

void SphericalHarmonicTransformationCache::updateFromQuaternion(
        const Eigen::Quaterniond& rotationQuaternion )
{
    std::complex< double > cayleyKleinA;
    std::complex< double > cayleyKleinB;

    convertQuaterionToCayleyKleinParameters( rotationQuaternion, cayleyKleinA, cayleyKleinB );
    updateFromCayleyKleinParameters( cayleyKleinA, cayleyKleinB );
}

void SphericalHarmonicTransformationCache::transformCoefficientsAtDegree(
        const Eigen::MatrixXd& originalCosineCoefficients,
        const Eigen::MatrixXd& originalSineCoefficients,
        Eigen::MatrixXd& currentCosineCoefficients,
        Eigen::MatrixXd& currentSineCoefficients,
        const bool areCoefficientsNormalized )
{
    currentCosineCoefficients.setZero( originalCosineCoefficients.rows( ), originalCosineCoefficients.cols( ) );
    currentSineCoefficients.setZero( originalSineCoefficients.rows( ), originalSineCoefficients.cols( ) );

    double currentMultiplier;
    double orderMMultiplier;

    for( unsigned int l = 0; l < originalCosineCoefficients.rows( ); l++ )
    {
        for( unsigned int m = 0; ( m <= l && m < originalCosineCoefficients.cols( ) ); m++ )
        {
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt( boost::math::factorial< double >( l - m ) / boost::math::factorial< double >( l + m ) );
                currentMultiplier = orderMMultiplier;
            }
            else
            {
                orderMMultiplier = ( m == 0 ? 1.0 : 1.0 / std::sqrt( 2.0 ) );
                currentMultiplier = orderMMultiplier;
            }

            std::complex< double > orderZeroDFunction = getWignerFunctionValue( l, m, 0 );

//            std::cout<<l<<" "<<m<<" "<<0<<" "<<orderZeroDFunction<<std::endl;

            currentCosineCoefficients( l, m ) +=
                    ( currentMultiplier ) * orderZeroDFunction.real( ) * originalCosineCoefficients( l, 0 );
            currentSineCoefficients( l, m ) +=
                    ( currentMultiplier ) * orderZeroDFunction.imag( ) * originalCosineCoefficients( l, 0 );


            for( int k = 1; k <= l; k++ )
            {
                if( !areCoefficientsNormalized )
                {
                    currentMultiplier = ::sqrt( boost::math::factorial< double >( l + k ) / boost::math::factorial< double >( l - k ) ) * orderMMultiplier;
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 ) * orderMMultiplier;
                }


                std::complex< double > orderKDFunction = getWignerFunctionValue( l, m, k );
                std::complex< double > orderMinusKDFunction = getWignerFunctionValue( l, m, -k );

//                std::cout<<l<<" "<<m<<" "<<k<<" "<<orderKDFunction<<std::endl;
//                std::cout<<l<<" "<<m<<" "<<-k<<" "<<orderMinusKDFunction<<std::endl;

                double signMultiplier = ( ( ( k % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );

                currentCosineCoefficients( l, m ) += 0.5 * currentMultiplier * (
                            ( signMultiplier * orderKDFunction.real( ) + orderMinusKDFunction.real( ) ) * originalCosineCoefficients( l, k )  +
                            ( signMultiplier * orderKDFunction.imag( ) - orderMinusKDFunction.imag( ) ) * originalSineCoefficients( l, k ) );

                currentSineCoefficients( l, m ) += 0.5 * currentMultiplier * (
                            ( signMultiplier * orderKDFunction.imag( ) + orderMinusKDFunction.imag( ) ) * originalCosineCoefficients( l, k )  +
                            ( -signMultiplier * orderKDFunction.real( ) + orderMinusKDFunction.real( ) ) * originalSineCoefficients( l, k ) );
            }

            currentCosineCoefficients( l, m ) = currentCosineCoefficients( l, m ) * ( ( ( m % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            currentSineCoefficients( l, m ) = currentSineCoefficients( l, m ) * ( ( ( ( m + 1 ) % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );


            if( m > 0 )
            {
                currentCosineCoefficients( l, m ) = currentCosineCoefficients( l, m ) * 2.0;
                currentSineCoefficients( l, m ) = currentSineCoefficients( l, m ) * 2.0;
            }
        }
    }
}

}

}


