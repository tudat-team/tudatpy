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
    wignerDMatricesCache_->updateMatrices( cayleyKleinA, cayleyKleinB );
}

void updateFromQuaternion(
        const Eigen::Quaterniond& rotationQuaternion )
{
    std::complex< double > cayleyKleinA;
    std::complex< double > cayleyKleinB;

    convertQuaterionToCayleyKleinParameters( rotationQuaternion, cayleyKleinA, cayleyKleinB );
    updateFromCayleyKleinParameters( cayleyKleinA, cayleyKleinB );
}

void SphericalHarmonicTransformationCache::updateFunctionPartials( )
{

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
    for( unsigned int l = 0; l < originalCosineCoefficients.rows( ); l++ )
    {
        for( int m = 0; ( m <= l && m < originalCosineCoefficients.cols( ) ); m++ )
        {
            if( !areCoefficientsNormalized )
            {
                currentMultiplier = factorials_[ l ] / factorials_[ l + m ];
            }
            else
            {
                currentMultiplier = factorials_[ l ]  /
                        std::sqrt( factorials_[ l + m ]  * factorials_[ l - m ] * ( ( m == 0 ) ? 1.0 : 2.0 ) );
            }

            currentCosineCoefficients( l, m ) +=
                    ( currentMultiplier ) * eFunctions_[ 0 + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].real( ) * originalCosineCoefficients( l, 0 );
            currentSineCoefficients( l, m ) +=
                    ( currentMultiplier ) * eFunctions_[ 0 + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].imag( ) * originalCosineCoefficients( l, 0 );


            for( int k = 1; k <= l; k++ )
            {
                if( !areCoefficientsNormalized )
                {
                    currentMultiplier = factorials_[ l + k ] / factorials_[ l + m ];
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 * ( factorials_[ l + k ] * factorials_[ l - k ] ) /
                            ( factorials_[ l + m ]  * factorials_[ l - m ] * ( ( m == 0 ) ? 1.0 : 2.0 ) ) );

                }
                currentCosineCoefficients( l, m ) += 0.5 * currentMultiplier * (
                            ( ( ( ( k % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) ) *eFunctions_[ -k + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].real( ) +
                        eFunctions_[ k + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].real( ) ) * originalCosineCoefficients( l, k )  +
                        ( ( ( ( k % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) ) *eFunctions_[ -k + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].imag( ) -
                        eFunctions_[ k + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].imag( ) ) * originalSineCoefficients( l, k ) );

                currentSineCoefficients( l, m ) += 0.5 * currentMultiplier * (
                            ( ( ( ( k % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) ) * eFunctions_[ -k + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].imag( ) +
                        eFunctions_[ k + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].imag( ) ) * originalCosineCoefficients( l, k )  +
                        ( ( ( ( ( k + 1 ) % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) ) *eFunctions_[ -k + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].real( ) +
                        eFunctions_[ k + maximumDegree_ + kSize * ( -m + maximumDegree_ + mSize * l ) ].real( ) ) * originalSineCoefficients( l, k ) );
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

void SphericalHarmonicTransformationCache::getPartialDerivativesOfTransformedCoefficientsWrtEulerAngles(
        const Eigen::MatrixXd& originalCosineCoefficients,
        const Eigen::MatrixXd& originalSineCoefficients,
        std::vector< Eigen::MatrixXd >& currentCosineCoefficients,
        std::vector< Eigen::MatrixXd >& currentSineCoefficients,
        const bool areCoefficientsNormalized )
{
    currentCosineCoefficients.resize( 3 );
    currentSineCoefficients.resize( 3 );

    for( unsigned int i = 0; i < 3; i++ )
    {
        currentCosineCoefficients[ i ].setZero( originalCosineCoefficients.rows( ), originalCosineCoefficients.cols( ) );
        currentSineCoefficients[ i ].setZero( originalSineCoefficients.rows( ), originalSineCoefficients.cols( ) );
    }

    double currentMultiplier, realScaling, imaginaryScaling, realScalingConj, imaginaryScalingConj;
    std::complex< double > currentEFunctionPartial;


}

}

}


