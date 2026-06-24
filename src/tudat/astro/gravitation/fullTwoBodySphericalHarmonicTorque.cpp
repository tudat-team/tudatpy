#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"

#include <boost/math/special_functions/factorials.hpp>

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace tudat
{

namespace gravitation
{

namespace
{

//! Apply the angular momentum operator to one Wigner-D coefficient entry.
/*!
 * Implements the Cartesian form of \hat{J}(D^l_{m,k}) used in the torque evaluation
 * of Dirkx et al. (2019), Eq. (67), with torque definition from Dirkx et al. (2019), Eq. (60).
 */
Eigen::Vector3cd computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices( const std::vector< Eigen::MatrixXcd >& wignerMatrices,
                                                                                      const int degree,
                                                                                      const int originalOrder,
                                                                                      const int newOrder )
{
    const std::complex< double > imaginaryUnit( 0.0, 1.0 );
    const double inverseSquareRootTwo = 1.0 / std::sqrt( 2.0 );

    const auto getWignerCoefficient = [ & ]( const int requestedOrderM, const int requestedOrderK ) {
        if( std::abs( requestedOrderM ) > degree || std::abs( requestedOrderK ) > degree )
        {
            return std::complex< double >( 0.0, 0.0 );
        }
        return wignerMatrices.at( degree )( requestedOrderM + degree, requestedOrderK + degree );
    };

    const double plusScaling =
            std::sqrt( std::max( 0.0, static_cast< double >( degree * ( degree + 1 ) - originalOrder * ( originalOrder - 1 ) ) ) ) /
            std::sqrt( 2.0 );
    const double minusScaling =
            std::sqrt( std::max( 0.0, static_cast< double >( degree * ( degree + 1 ) - originalOrder * ( originalOrder + 1 ) ) ) ) /
            std::sqrt( 2.0 );

    // Dirkx et al. (2019), Eq. (67): Cartesian components of \hat{J}(D^l_{m,k}) assembled from ladder-operator contributions.
    const std::complex< double > angularMomentumPlus = imaginaryUnit * plusScaling * getWignerCoefficient( originalOrder - 1, newOrder );
    const std::complex< double > angularMomentumMinus =
            imaginaryUnit * ( -minusScaling ) * getWignerCoefficient( originalOrder + 1, newOrder );
    const std::complex< double > angularMomentumZero =
            imaginaryUnit * static_cast< double >( -originalOrder ) * getWignerCoefficient( originalOrder, newOrder );

    Eigen::Vector3cd angularMomentumInCartesianBasis;
    angularMomentumInCartesianBasis( 0 ) = ( angularMomentumMinus - angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 1 ) = imaginaryUnit * ( angularMomentumMinus + angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 2 ) = angularMomentumZero;
    // Returned \hat{J}-mapped term is used directly in the Dirkx et al. (2019), Eq. (60) torque evaluation.
    return angularMomentumInCartesianBasis;
}

}  // namespace

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getBody2TorqueCombinationsToUse(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse )
{
    // Build unique combinations required for body-2 spin torque accumulation (Dirkx et al. (2019), Eq. (67) term set).
    std::set< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > body2TorqueCombinationSet;
    for( const auto& combination : coefficientCombinationsToUse )
    {
        const unsigned int degreeOfBody1ForCombination = std::get< 0 >( combination );
        const unsigned int orderOfBody1ForCombination = std::get< 1 >( combination );
        const unsigned int degreeOfBody2ForCombination = std::get< 2 >( combination );
        for( unsigned int orderOfBody2ForTorque = 0; orderOfBody2ForTorque <= degreeOfBody2ForCombination; orderOfBody2ForTorque++ )
        {
            body2TorqueCombinationSet.insert( std::make_tuple(
                    degreeOfBody1ForCombination, orderOfBody1ForCombination, degreeOfBody2ForCombination, orderOfBody2ForTorque ) );
        }
    }

    return std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >( body2TorqueCombinationSet.begin( ),
                                                                                                body2TorqueCombinationSet.end( ) );
}

void computeTransformedAngularMomentumCoefficientsFromWignerMatrices(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum )
{
    // Explicit-matrix path (used by quaternion-derivative partials).
    // Dirkx et al. (2019), Eq. (67): computes \hat{J}-mapped coefficient fields; these feed
    // the Dirkx et al. (2019), Eq. (60) torque relation.
    for( int i = 0; i < 3; i++ )
    {
        transformedCosineCoefficientsBody2AngularMomentum.at( i ).setZero( cosineCoefficientsBody2.rows( ),
                                                                           cosineCoefficientsBody2.cols( ) );
        transformedSineCoefficientsBody2AngularMomentum.at( i ).setZero( sineCoefficientsBody2.rows( ), sineCoefficientsBody2.cols( ) );
    }

    for( int degree = 0; degree < cosineCoefficientsBody2.rows( ); degree++ )
    {
        for( int orderM = 0; ( orderM <= degree && orderM < cosineCoefficientsBody2.cols( ) ); orderM++ )
        {
            double orderMMultiplier;
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt( boost::math::factorial< double >( degree - orderM ) /
                                              boost::math::factorial< double >( degree + orderM ) );
            }
            else
            {
                orderMMultiplier = ( orderM == 0 ? 1.0 : 1.0 / std::sqrt( 2.0 ) );
            }

            const Eigen::Vector3cd orderZeroAngularMomentumD =
                    computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices( wignerMatrices, degree, orderM, 0 );
            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).real( ) * cosineCoefficientsBody2( degree, 0 );
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).imag( ) * cosineCoefficientsBody2( degree, 0 );
            }

            const int maximumOrderAtDegree = std::min( degree, static_cast< int >( cosineCoefficientsBody2.cols( ) - 1 ) );
            for( int orderK = 1; orderK <= maximumOrderAtDegree; orderK++ )
            {
                double currentMultiplier;
                if( !areCoefficientsNormalized )
                {
                    currentMultiplier = std::sqrt( boost::math::factorial< double >( degree + orderK ) /
                                                   boost::math::factorial< double >( degree - orderK ) ) *
                            orderMMultiplier;
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 ) * orderMMultiplier;
                }

                const double signMultiplier = ( ( ( orderK % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
                const Eigen::Vector3cd positiveOrderAngularMomentumD =
                        computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices( wignerMatrices, degree, orderM, orderK );
                const Eigen::Vector3cd negativeOrderAngularMomentumD =
                        computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices( wignerMatrices, degree, orderM, -orderK );

                for( int i = 0; i < 3; i++ )
                {
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) += 0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) - negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );

                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) += 0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) + negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );
                }
            }

            double cosineFinalScaling = ( ( ( orderM % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            double sineFinalScaling = ( ( ( ( orderM + 1 ) % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            if( orderM > 0 )
            {
                cosineFinalScaling *= 2.0;
                sineFinalScaling *= 2.0;
            }

            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= cosineFinalScaling;
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= sineFinalScaling;
            }
        }
    }
}

void computeTransformedAngularMomentumCoefficientsFromWignerCache(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum )
{
    // Cached path (used during nominal model evaluation).
    // Dirkx et al. (2019), Eq. (67): computes \hat{J}-mapped coefficient fields; these feed
    // the Dirkx et al. (2019), Eq. (60) torque relation.
    for( int i = 0; i < 3; i++ )
    {
        transformedCosineCoefficientsBody2AngularMomentum.at( i ).setZero( cosineCoefficientsBody2.rows( ),
                                                                           cosineCoefficientsBody2.cols( ) );
        transformedSineCoefficientsBody2AngularMomentum.at( i ).setZero( sineCoefficientsBody2.rows( ), sineCoefficientsBody2.cols( ) );
    }

    for( int degree = 0; degree < cosineCoefficientsBody2.rows( ); degree++ )
    {
        for( int orderM = 0; ( orderM <= degree && orderM < cosineCoefficientsBody2.cols( ) ); orderM++ )
        {
            double orderMMultiplier;
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt( boost::math::factorial< double >( degree - orderM ) /
                                              boost::math::factorial< double >( degree + orderM ) );
            }
            else
            {
                orderMMultiplier = ( orderM == 0 ? 1.0 : 1.0 / std::sqrt( 2.0 ) );
            }

            const Eigen::Vector3cd orderZeroAngularMomentumD =
                    wignerCache->getAngularMomentumOperatorOnWignerCoefficient( degree, orderM, 0 );
            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).real( ) * cosineCoefficientsBody2( degree, 0 );
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).imag( ) * cosineCoefficientsBody2( degree, 0 );
            }

            const int maximumOrderAtDegree = std::min( degree, static_cast< int >( cosineCoefficientsBody2.cols( ) - 1 ) );
            for( int orderK = 1; orderK <= maximumOrderAtDegree; orderK++ )
            {
                double currentMultiplier;
                if( !areCoefficientsNormalized )
                {
                    currentMultiplier = std::sqrt( boost::math::factorial< double >( degree + orderK ) /
                                                   boost::math::factorial< double >( degree - orderK ) ) *
                            orderMMultiplier;
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 ) * orderMMultiplier;
                }

                const double signMultiplier = ( ( ( orderK % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
                const Eigen::Vector3cd positiveOrderAngularMomentumD =
                        wignerCache->getAngularMomentumOperatorOnWignerCoefficient( degree, orderM, orderK );
                const Eigen::Vector3cd negativeOrderAngularMomentumD =
                        wignerCache->getAngularMomentumOperatorOnWignerCoefficient( degree, orderM, -orderK );

                for( int i = 0; i < 3; i++ )
                {
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) += 0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) - negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );

                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) += 0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) + negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );
                }
            }

            double cosineFinalScaling = ( ( ( orderM % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            double sineFinalScaling = ( ( ( ( orderM + 1 ) % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            if( orderM > 0 )
            {
                cosineFinalScaling *= 2.0;
                sineFinalScaling *= 2.0;
            }

            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= cosineFinalScaling;
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= sineFinalScaling;
            }
        }
    }
}

void FullTwoBodySphericalHarmonicTorque::computeTransformedAngularMomentumCoefficients(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum )
{
    // Forward to cache-based utility; kept as class entry point for reuse by analytical partials.
    // Dirkx et al. (2019), Eq. (67): this function provides the \hat{J}-mapped coefficients used in Eq. (60).
    computeTransformedAngularMomentumCoefficientsFromWignerCache( cosineCoefficientsBody2,
                                                                  sineCoefficientsBody2,
                                                                  wignerCache,
                                                                  areCoefficientsNormalized,
                                                                  transformedCosineCoefficientsBody2AngularMomentum,
                                                                  transformedSineCoefficientsBody2AngularMomentum );
}

}  // namespace gravitation

}  // namespace tudat
