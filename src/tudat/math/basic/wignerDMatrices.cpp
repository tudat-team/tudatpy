#include <algorithm>
#include <array>
#include <cmath>

#include "tudat/math/basic/wignerDMatrices.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace basic_mathematics
{

namespace
{

std::array< Eigen::MatrixXcd, 4 > computeDerivativeOfDegreeOneWignerDMatrixWrtQuaternion( const Eigen::Vector4d& quaternionVector )
{
    const std::complex< double > imaginaryUnit( 0.0, 1.0 );
    const std::complex< double > cayleyKleinA = std::complex< double >( quaternionVector( 0 ), -quaternionVector( 3 ) );
    const std::complex< double > cayleyKleinB = std::complex< double >( quaternionVector( 2 ), -quaternionVector( 1 ) );
    const std::complex< double > cayleyKleinAConjugate = std::conj( cayleyKleinA );
    const std::complex< double > cayleyKleinBConjugate = std::conj( cayleyKleinB );

    std::array< std::complex< double >, 4 > partialOfCayleyKleinA;
    std::array< std::complex< double >, 4 > partialOfCayleyKleinB;
    std::array< std::complex< double >, 4 > partialOfCayleyKleinAConjugate;
    std::array< std::complex< double >, 4 > partialOfCayleyKleinBConjugate;

    partialOfCayleyKleinA.at( 0 ) = 1.0;
    partialOfCayleyKleinA.at( 1 ) = 0.0;
    partialOfCayleyKleinA.at( 2 ) = 0.0;
    partialOfCayleyKleinA.at( 3 ) = -imaginaryUnit;

    partialOfCayleyKleinB.at( 0 ) = 0.0;
    partialOfCayleyKleinB.at( 1 ) = -imaginaryUnit;
    partialOfCayleyKleinB.at( 2 ) = 1.0;
    partialOfCayleyKleinB.at( 3 ) = 0.0;

    for( int i = 0; i < 4; i++ )
    {
        partialOfCayleyKleinAConjugate.at( i ) = std::conj( partialOfCayleyKleinA.at( i ) );
        partialOfCayleyKleinBConjugate.at( i ) = std::conj( partialOfCayleyKleinB.at( i ) );
    }

    std::array< Eigen::MatrixXcd, 4 > derivatives;
    for( int i = 0; i < 4; i++ )
    {
        derivatives.at( i ) = Eigen::MatrixXcd::Zero( 3, 3 );
        derivatives.at( i )( 2, 2 ) = 2.0 * cayleyKleinA * partialOfCayleyKleinA.at( i );
        derivatives.at( i )( 2, 1 ) = -std::sqrt( 2.0 ) *
                ( partialOfCayleyKleinA.at( i ) * cayleyKleinBConjugate + cayleyKleinA * partialOfCayleyKleinBConjugate.at( i ) );
        derivatives.at( i )( 2, 0 ) = 2.0 * cayleyKleinBConjugate * partialOfCayleyKleinBConjugate.at( i );

        derivatives.at( i )( 1, 2 ) =
                std::sqrt( 2.0 ) * ( partialOfCayleyKleinA.at( i ) * cayleyKleinB + cayleyKleinA * partialOfCayleyKleinB.at( i ) );
        derivatives.at( i )( 1, 1 ) = ( i == 0 ? 2.0 * quaternionVector( 0 )
                                               : ( i == 1 ? -2.0 * quaternionVector( 1 )
                                                          : ( i == 2 ? -2.0 * quaternionVector( 2 ) : 2.0 * quaternionVector( 3 ) ) ) );
        derivatives.at( i )( 1, 0 ) = -std::sqrt( 2.0 ) *
                ( partialOfCayleyKleinAConjugate.at( i ) * cayleyKleinBConjugate +
                  cayleyKleinAConjugate * partialOfCayleyKleinBConjugate.at( i ) );

        derivatives.at( i )( 0, 2 ) = 2.0 * cayleyKleinB * partialOfCayleyKleinB.at( i );
        derivatives.at( i )( 0, 1 ) = std::sqrt( 2.0 ) *
                ( partialOfCayleyKleinAConjugate.at( i ) * cayleyKleinB + cayleyKleinAConjugate * partialOfCayleyKleinB.at( i ) );
        derivatives.at( i )( 0, 0 ) = 2.0 * cayleyKleinAConjugate * partialOfCayleyKleinAConjugate.at( i );
    }
    return derivatives;
}

double getWignerRecursionCoefficientMinusOne( const int degree, const int rowIndex, const int columnIndex )
{
    const int orderM = rowIndex - degree;
    const int orderK = columnIndex - degree;
    return std::sqrt( static_cast< double >( ( degree + orderK ) * ( degree + orderK - 1 ) ) /
                      static_cast< double >( ( degree + orderM ) * ( degree + orderM - 1 ) ) );
}

double getWignerRecursionCoefficientZero( const int degree, const int rowIndex, const int columnIndex )
{
    const int orderM = rowIndex - degree;
    const int orderK = columnIndex - degree;
    return std::sqrt( static_cast< double >( 2 * ( degree + orderK ) * ( degree - orderK ) ) /
                      static_cast< double >( ( degree + orderM ) * ( degree + orderM - 1 ) ) );
}

double getWignerRecursionCoefficientOne( const int degree, const int rowIndex, const int columnIndex )
{
    const int orderM = rowIndex - degree;
    const int orderK = columnIndex - degree;
    return std::sqrt( static_cast< double >( ( degree - orderK ) * ( degree - orderK - 1 ) ) /
                      static_cast< double >( ( degree + orderM ) * ( degree + orderM - 1 ) ) );
}

}  // namespace

//! Constructor
WignerDMatricesCache::WignerDMatricesCache( const int maximumDegree ):
    maximumDegree_( maximumDegree ), areAngularMomentumOperatorsUpdated_( false )
{
    wignerDMatrices_.resize( maximumDegree + 1 );
    angularMomentumOperatorsX_.resize( maximumDegree + 1 );
    angularMomentumOperatorsY_.resize( maximumDegree + 1 );
    angularMomentumOperatorsZ_.resize( maximumDegree + 1 );

    for( int l = 0; l <= maximumDegree; l++ )
    {
        wignerDMatrices_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
        angularMomentumOperatorsX_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
        angularMomentumOperatorsY_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
        angularMomentumOperatorsZ_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
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
        wignerDMatrices_[ 1 ]( 1, 1 ) = std::norm( currentCayleyKleinA_ ) - std::norm( currentCayleyKleinB_ );
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
                        wignerDMatrices_[ l ]( i, j ) += coefficientsIndexZero_[ l ]( i, j ) * wignerDMatrices_[ 1 ]( 2, 1 ) *
                                wignerDMatrices_[ l - 1 ]( i - 2, j - 1 );
                    }
                    if( j < 2 * l - 1 )
                    {
                        wignerDMatrices_[ l ]( i, j ) +=
                                coefficientsIndexOne_[ l ]( i, j ) * wignerDMatrices_[ 1 ]( 2, 0 ) * wignerDMatrices_[ l - 1 ]( i - 2, j );
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
                wignerDMatrices_[ l ]( i, j ) =
                        ( ( ( ( m - k ) % 2 ) == 0 ) ? 1.0 : -1.0 ) * std::conj( wignerDMatrices_[ l ]( -m + l, -k + l ) );
            }
        }
    }

    areAngularMomentumOperatorsUpdated_ = false;
}

void WignerDMatricesCache::computeAngularMomentumOperators( )
{
    const std::complex< double > imaginaryUnit( 0.0, 1.0 );
    const double inverseSquareRootTwo = 1.0 / std::sqrt( 2.0 );

    for( int degree = 0; degree <= maximumDegree_; degree++ )
    {
        for( int rowIndex = 0; rowIndex <= 2 * degree; rowIndex++ )
        {
            const int orderM = rowIndex - degree;

            const double plusScaling =
                    std::sqrt( std::max( 0.0, static_cast< double >( degree * ( degree + 1 ) - orderM * ( orderM - 1 ) ) ) ) /
                    std::sqrt( 2.0 );
            const double minusScaling =
                    std::sqrt( std::max( 0.0, static_cast< double >( degree * ( degree + 1 ) - orderM * ( orderM + 1 ) ) ) ) /
                    std::sqrt( 2.0 );

            for( int columnIndex = 0; columnIndex <= 2 * degree; columnIndex++ )
            {
                const int orderK = columnIndex - degree;

                const auto getWignerCoefficient = [ & ]( const int requestedOrderM, const int requestedOrderK ) {
                    if( std::abs( requestedOrderM ) > degree || std::abs( requestedOrderK ) > degree )
                    {
                        return std::complex< double >( 0.0, 0.0 );
                    }
                    return wignerDMatrices_[ degree ]( requestedOrderM + degree, requestedOrderK + degree );
                };

                const std::complex< double > angularMomentumPlus = imaginaryUnit * plusScaling * getWignerCoefficient( orderM - 1, orderK );
                const std::complex< double > angularMomentumMinus =
                        imaginaryUnit * ( -minusScaling ) * getWignerCoefficient( orderM + 1, orderK );
                const std::complex< double > angularMomentumZero =
                        imaginaryUnit * static_cast< double >( -orderM ) * getWignerCoefficient( orderM, orderK );

                angularMomentumOperatorsX_[ degree ]( rowIndex, columnIndex ) =
                        ( angularMomentumMinus - angularMomentumPlus ) * inverseSquareRootTwo;
                angularMomentumOperatorsY_[ degree ]( rowIndex, columnIndex ) =
                        imaginaryUnit * ( angularMomentumMinus + angularMomentumPlus ) * inverseSquareRootTwo;
                angularMomentumOperatorsZ_[ degree ]( rowIndex, columnIndex ) = angularMomentumZero;
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
                if( ( l + m ) >= 2 )
                {
                    const double denominator = static_cast< double >( ( l + m ) * ( l + m - 1 ) );
                    coefficientsIndexMinusOne_[ l ]( i, j ) = std::sqrt( static_cast< double >( ( l + k ) * ( l + k - 1 ) ) / denominator );
                    coefficientsIndexZero_[ l ]( i, j ) = std::sqrt( static_cast< double >( 2 * ( l + k ) * ( l - k ) ) / denominator );
                    coefficientsIndexOne_[ l ]( i, j ) = std::sqrt( static_cast< double >( ( l - k ) * ( l - k - 1 ) ) / denominator );
                }
            }
        }
    }
}

void computeDerivativeOfWignerDMatricesWrtQuaternion( const Eigen::Quaterniond& rotation,
                                                      const std::shared_ptr< WignerDMatricesCache >& wignerCache,
                                                      std::array< std::vector< Eigen::MatrixXcd >, 4 >& derivatives )
{
    const std::vector< Eigen::MatrixXcd >& wignerMatrices = wignerCache->getWignerDMatrices( );
    const int maximumDegree = static_cast< int >( wignerMatrices.size( ) ) - 1;

    for( int i = 0; i < 4; i++ )
    {
        derivatives.at( i ).resize( maximumDegree + 1 );
        derivatives.at( i ).at( 0 ).setZero( 1, 1 );
    }

    if( maximumDegree == 0 )
    {
        return;
    }

    const Eigen::Vector4d quaternionVector = linear_algebra::convertQuaternionToVectorFormat( rotation );
    const std::array< Eigen::MatrixXcd, 4 > degreeOneDerivatives =
            computeDerivativeOfDegreeOneWignerDMatrixWrtQuaternion( quaternionVector );

    for( int i = 0; i < 4; i++ )
    {
        derivatives.at( i ).at( 1 ) = degreeOneDerivatives.at( i );
    }

    if( maximumDegree == 1 )
    {
        return;
    }

    const Eigen::MatrixXcd& degreeOneWignerMatrix = wignerMatrices.at( 1 );
    for( int degree = 2; degree <= maximumDegree; degree++ )
    {
        for( int i = 0; i < 4; i++ )
        {
            derivatives.at( i ).at( degree ).setZero( 2 * degree + 1, 2 * degree + 1 );
        }

        for( int rowIndex = degree; rowIndex <= 2 * degree; rowIndex++ )
        {
            for( int columnIndex = 0; columnIndex <= 2 * degree; columnIndex++ )
            {
                if( rowIndex - 2 < 0 )
                {
                    continue;
                }

                for( int i = 0; i < 4; i++ )
                {
                    std::complex< double > derivativeEntry = std::complex< double >( 0.0, 0.0 );

                    if( columnIndex > 1 )
                    {
                        const double coefficient = getWignerRecursionCoefficientMinusOne( degree, rowIndex, columnIndex );
                        derivativeEntry += coefficient *
                                ( degreeOneDerivatives.at( i )( 2, 2 ) * wignerMatrices.at( degree - 1 )( rowIndex - 2, columnIndex - 2 ) +
                                  degreeOneWignerMatrix( 2, 2 ) * derivatives.at( i ).at( degree - 1 )( rowIndex - 2, columnIndex - 2 ) );
                    }
                    if( columnIndex > 0 && columnIndex <= 2 * degree - 1 )
                    {
                        const double coefficient = getWignerRecursionCoefficientZero( degree, rowIndex, columnIndex );
                        derivativeEntry += coefficient *
                                ( degreeOneDerivatives.at( i )( 2, 1 ) * wignerMatrices.at( degree - 1 )( rowIndex - 2, columnIndex - 1 ) +
                                  degreeOneWignerMatrix( 2, 1 ) * derivatives.at( i ).at( degree - 1 )( rowIndex - 2, columnIndex - 1 ) );
                    }
                    if( columnIndex < 2 * degree - 1 )
                    {
                        const double coefficient = getWignerRecursionCoefficientOne( degree, rowIndex, columnIndex );
                        derivativeEntry += coefficient *
                                ( degreeOneDerivatives.at( i )( 2, 0 ) * wignerMatrices.at( degree - 1 )( rowIndex - 2, columnIndex ) +
                                  degreeOneWignerMatrix( 2, 0 ) * derivatives.at( i ).at( degree - 1 )( rowIndex - 2, columnIndex ) );
                    }

                    derivatives.at( i ).at( degree )( rowIndex, columnIndex ) = derivativeEntry;
                }
            }
        }

        for( int rowIndex = 0; rowIndex < degree; rowIndex++ )
        {
            const int orderM = rowIndex - degree;
            for( int columnIndex = 0; columnIndex <= 2 * degree; columnIndex++ )
            {
                const int orderK = columnIndex - degree;
                const double signMultiplier = ( ( ( orderM - orderK ) % 2 ) == 0 ) ? 1.0 : -1.0;
                for( int i = 0; i < 4; i++ )
                {
                    derivatives.at( i ).at( degree )( rowIndex, columnIndex ) =
                            signMultiplier * std::conj( derivatives.at( i ).at( degree )( -orderM + degree, -orderK + degree ) );
                }
            }
        }
    }
}

}  // namespace basic_mathematics

}  // namespace tudat
