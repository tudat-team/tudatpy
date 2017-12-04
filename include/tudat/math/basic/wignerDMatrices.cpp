#include <iostream>

#include <Tudat/Mathematics/BasicMathematics/wignerDMatrices.h>

namespace tudat
{

namespace basic_mathematics
{

//! Constructor
WignerDMatricesCache::WignerDMatricesCache( const int maximumDegree ):
    maximumDegree_( maximumDegree )
{
    wignerDMatrices_.resize( maximumDegree + 1 );

    for( int l = 0; l <= maximumDegree; l++ )
    {
        wignerDMatrices_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
    }
    wignerDMatrices_[ 0 ]( 0, 0 ) = std::complex< double >( 1.0, 0.0 );

    computeCoefficients( );
}

//! Function to update contents of this object to new orientation
void WignerDMatricesCache::updateMatrices( const std::complex< double > cayleyKleinA, const std::complex< double > cayleyKleinB )
{
    // Set current orientation
    currentCayleyKleinA_ = cayleyKleinA;
    currentCayleyKleinB_ = cayleyKleinB;
    currentCayleyKleinAConjugate_ = std::conj( currentCayleyKleinA_ );
    currentCayleyKleinBConjugate_ = std::conj( currentCayleyKleinB_ );

    // Explicitly compute coefficients at degree 1
    if( maximumDegree_ > 0 )
    {
        wignerDMatrices_[ 1 ]( 2, 2 ) = currentCayleyKleinA_ * currentCayleyKleinA_;
        wignerDMatrices_[ 1 ]( 2, 1 ) = -std::sqrt( 2.0 ) * currentCayleyKleinA_ * currentCayleyKleinBConjugate_;
        wignerDMatrices_[ 1 ]( 2, 0 ) = currentCayleyKleinBConjugate_ * currentCayleyKleinBConjugate_;
        wignerDMatrices_[ 1 ]( 1, 2 ) = std::sqrt( 2.0 ) * currentCayleyKleinA_ * currentCayleyKleinB_;
        wignerDMatrices_[ 1 ]( 1, 1 ) = std::norm( currentCayleyKleinA_ ) * std::norm( currentCayleyKleinA_ ) -
                std::norm( currentCayleyKleinB_ ) * std::norm( currentCayleyKleinB_ );
        wignerDMatrices_[ 1 ]( 1, 0 ) = -std::sqrt( 2.0 ) * currentCayleyKleinAConjugate_ * currentCayleyKleinBConjugate_;
        wignerDMatrices_[ 1 ]( 0, 2 ) = currentCayleyKleinB_ * currentCayleyKleinB_;
        wignerDMatrices_[ 1 ]( 0, 1 ) = std::sqrt( 2.0 ) * currentCayleyKleinAConjugate_ * currentCayleyKleinB_;
        wignerDMatrices_[ 1 ]( 0, 0 ) = currentCayleyKleinAConjugate_ * currentCayleyKleinAConjugate_;
    }

    // Recursively compute coefficients at degree >1
    for( int l = 2; l <= maximumDegree_; l++ )
    {
        for( int i = l; i <= 2 * l; i++ )
        {
            for( int j = 0; j <= 2 * l; j++ )
            {                
                if( i - 2 >= 0 )
                {
                    wignerDMatrices_[ l ]( i, j ) = 0.0;

                    // For each part in equation, check if contribution is non-zer0
                    if( j > 1 )
                    {
                       wignerDMatrices_[ l ]( i, j ) += coefficientsIndexMinusOne_[ l ]( i, j ) * wignerDMatrices_[ 1 ]( 2, 2 ) *
                                wignerDMatrices_[ l - 1 ]( i - 2, j - 2 );
                    }
                    if( ( j > 0 ) && ( j <= 2 * l - 1 ) )
                    {
                        wignerDMatrices_[ l ]( i, j ) +=
                                coefficientsIndexZero_[ l ]( i, j ) * wignerDMatrices_[ 1 ]( 2, 1 ) *
                                wignerDMatrices_[ l - 1 ]( i - 2, j - 1 );
                    }
                    if( j < 2 * l - 1 )
                    {
                        wignerDMatrices_[ l ]( i, j ) += coefficientsIndexOne_[ l ]( i, j ) * wignerDMatrices_[ 1 ]( 2, 0 ) *
                                wignerDMatrices_[ l - 1 ]( i - 2, j );
                    }
                }
                else
                {
                    wignerDMatrices_[ l ]( i, j ) = std::complex< double >( 0.0, 0.0 );
                }
            }
        }

        // Use symmetry relation to compute coefficients for negative m
        int m = 0, k = 0;
        for( int i = 0; i < l; i++ )
        {
            m = i - l;
            for( int j = 0; j <= 2 * l; j++ )
            {
                k = j - l;
                wignerDMatrices_[ l ]( i, j ) = ( ( ( ( m - k ) % 2 ) == 0 ) ? 1.0 : -1.0 ) *
                        std::conj( wignerDMatrices_[ l ]( -m + l, -k + l ) );
            }

        }
    }
}

//! Function to precompute the coefficients used on the recursive formulation for Wigner D-matrices
void WignerDMatricesCache::computeCoefficients( )
{
    coefficientsIndexMinusOne_.resize( maximumDegree_ + 1 );
    coefficientsIndexZero_.resize( maximumDegree_ + 1 );
    coefficientsIndexOne_.resize( maximumDegree_ + 1 );

    int m = 0, k = 0;
    for( int l = 0; l <= maximumDegree_; l++ )
    {
        // Allocate size for coefficients of current degree
        coefficientsIndexMinusOne_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
        coefficientsIndexZero_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
        coefficientsIndexOne_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );

        // Compute coefficients at current degree.
        for( int i = 0; i <= 2 * l; i++ )
        {
            m = i - l;
            for( int j = 0; j <= 2 * l; j++ )
            {
                k = j - l;
                coefficientsIndexMinusOne_[ l ]( i, j ) = std::sqrt(
                            static_cast< double >( ( l + k ) * ( l + k - 1 ) ) /
                            static_cast< double >( ( l + m ) * ( l + m - 1 ) ) );
                coefficientsIndexZero_[ l ]( i, j ) = std::sqrt(
                            static_cast< double >( 2 * ( l + k ) * ( l - k ) ) /
                            static_cast< double >( ( l + m ) * ( l + m - 1 ) ) );
                coefficientsIndexOne_[ l ]( i, j ) = std::sqrt(
                            static_cast< double >( ( l - k ) * ( l - k - 1 ) ) /
                            static_cast< double >( ( l + m ) * ( l + m - 1 ) ) );
            }
        }
    }

}

} // namespace basic_mathematics

} // namespace tudat
