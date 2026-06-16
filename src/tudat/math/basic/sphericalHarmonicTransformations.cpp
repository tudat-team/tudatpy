#include <boost/math/special_functions/binomial.hpp>
#include <iostream>
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"

#include "tudat/math/basic/sphericalHarmonicTransformations.h"
#include "tudat/math/basic/cayleyKleinParameters.h"

namespace tudat
{

namespace basic_mathematics
{

//! Function to update Wigner D-matrices for current orientation, parameterized by Cayley-Klein parameters
void SphericalHarmonicTransformationCache::updateFromCayleyKleinParameters( const std::complex< double > cayleyKleinA,
                                                                            const std::complex< double > cayleyKleinB )
{
    wignerDMatricesCache_->updateMatrices( cayleyKleinA, cayleyKleinB );
}

//! Function to update Wigner D-matrices for current orientation, parameterized by quaternion
void SphericalHarmonicTransformationCache::updateFromQuaternion( const Eigen::Quaterniond& rotationQuaternion )
{
    std::complex< double > cayleyKleinA;
    std::complex< double > cayleyKleinB;

    convertQuaterionToCayleyKleinParameters( rotationQuaternion, cayleyKleinA, cayleyKleinB );
    updateFromCayleyKleinParameters( cayleyKleinA, cayleyKleinB );
}

void transformSphericalHarmonicCoefficientsWithWignerD( const Eigen::MatrixXd& originalCosineCoefficients,
                                                        const Eigen::MatrixXd& originalSineCoefficients,
                                                        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
                                                        Eigen::MatrixXd& transformedCosineCoefficients,
                                                        Eigen::MatrixXd& transformedSineCoefficients,
                                                        const bool areCoefficientsNormalized )
{
    transformedCosineCoefficients.setZero( originalCosineCoefficients.rows( ), originalCosineCoefficients.cols( ) );
    transformedSineCoefficients.setZero( originalSineCoefficients.rows( ), originalSineCoefficients.cols( ) );

    for( unsigned int l = 0; l < originalCosineCoefficients.rows( ); l++ )
    {
        for( unsigned int m = 0; ( m <= l && m < originalCosineCoefficients.cols( ) ); m++ )
        {
            double orderMMultiplier;
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt( boost::math::factorial< double >( l - m ) / boost::math::factorial< double >( l + m ) );
            }
            else
            {
                orderMMultiplier = ( m == 0 ? 1.0 : 1.0 / std::sqrt( 2.0 ) );
            }

            const std::complex< double > orderZeroDFunction = wignerMatrices.at( l )( m + l, l );
            transformedCosineCoefficients( l, m ) += orderMMultiplier * orderZeroDFunction.real( ) * originalCosineCoefficients( l, 0 );
            transformedSineCoefficients( l, m ) += orderMMultiplier * orderZeroDFunction.imag( ) * originalCosineCoefficients( l, 0 );

            for( int k = 1; k <= static_cast< int >( l ); k++ )
            {
                double currentMultiplier;
                if( !areCoefficientsNormalized )
                {
                    currentMultiplier = std::sqrt( boost::math::factorial< double >( l + k ) / boost::math::factorial< double >( l - k ) ) *
                            orderMMultiplier;
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 ) * orderMMultiplier;
                }

                const double signMultiplier = ( ( ( k % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
                const std::complex< double > orderKDFunction = wignerMatrices.at( l )( m + l, k + l );
                const std::complex< double > orderMinusKDFunction = wignerMatrices.at( l )( m + l, -k + l );

                transformedCosineCoefficients( l, m ) += 0.5 * currentMultiplier *
                        ( ( signMultiplier * orderKDFunction.real( ) + orderMinusKDFunction.real( ) ) * originalCosineCoefficients( l, k ) +
                          ( signMultiplier * orderKDFunction.imag( ) - orderMinusKDFunction.imag( ) ) * originalSineCoefficients( l, k ) );

                transformedSineCoefficients( l, m ) += 0.5 * currentMultiplier *
                        ( ( signMultiplier * orderKDFunction.imag( ) + orderMinusKDFunction.imag( ) ) * originalCosineCoefficients( l, k ) +
                          ( -signMultiplier * orderKDFunction.real( ) + orderMinusKDFunction.real( ) ) * originalSineCoefficients( l, k ) );
            }

            double cosineScaling = ( ( ( m % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            double sineScaling = ( ( ( ( m + 1 ) % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            if( m > 0 )
            {
                cosineScaling *= 2.0;
                sineScaling *= 2.0;
            }

            transformedCosineCoefficients( l, m ) *= cosineScaling;
            transformedSineCoefficients( l, m ) *= sineScaling;
        }
    }
}

//! Function to transform spherical harmonic coefficients using current wignerDMatricesCache_
void SphericalHarmonicTransformationCache::transformCoefficientsAtDegree( const Eigen::MatrixXd& originalCosineCoefficients,
                                                                          const Eigen::MatrixXd& originalSineCoefficients,
                                                                          Eigen::MatrixXd& currentCosineCoefficients,
                                                                          Eigen::MatrixXd& currentSineCoefficients,
                                                                          const bool areCoefficientsNormalized )
{
    transformSphericalHarmonicCoefficientsWithWignerD( originalCosineCoefficients,
                                                       originalSineCoefficients,
                                                       wignerDMatricesCache_->getWignerDMatrices( ),
                                                       currentCosineCoefficients,
                                                       currentSineCoefficients,
                                                       areCoefficientsNormalized );
}

}  // namespace basic_mathematics

}  // namespace tudat
