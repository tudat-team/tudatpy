#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"

#include <boost/math/special_functions/factorials.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

namespace tudat
{

namespace gravitation
{

namespace
{

bool isSingleC21ByC20Combination(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
                coefficientCombinationsToUse )
{
    return coefficientCombinationsToUse.size( ) == 1 &&
            std::get< 0 >( coefficientCombinationsToUse.at( 0 ) ) == 2 &&
            std::get< 1 >( coefficientCombinationsToUse.at( 0 ) ) == 1 &&
            std::get< 2 >( coefficientCombinationsToUse.at( 0 ) ) == 2 &&
            std::get< 3 >( coefficientCombinationsToUse.at( 0 ) ) == 0;
}

bool isUnitC21ByC20DebugCase(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
                coefficientCombinationsToUse,
        const Eigen::Vector3d& bodyFixedRelativePosition,
        const Eigen::MatrixXd& cosineCoefficientsOfBody1,
        const Eigen::MatrixXd& sineCoefficientsOfBody1,
        const Eigen::MatrixXd& cosineCoefficientsOfBody2,
        const Eigen::MatrixXd& sineCoefficientsOfBody2 )
{
    if( !isSingleC21ByC20Combination( coefficientCombinationsToUse ) )
    {
        return false;
    }

    if( cosineCoefficientsOfBody1.rows( ) < 3 || cosineCoefficientsOfBody1.cols( ) < 3 ||
        cosineCoefficientsOfBody2.rows( ) < 3 || cosineCoefficientsOfBody2.cols( ) < 3 )
    {
        return false;
    }

    const double tolerance = 1.0E-14;
    const bool isRequestedRelativePositionOnXAxis =
            std::fabs( std::fabs( bodyFixedRelativePosition.x( ) ) - 1.0 ) < tolerance &&
            std::fabs( bodyFixedRelativePosition.y( ) ) < tolerance &&
            std::fabs( bodyFixedRelativePosition.z( ) ) < tolerance;
    const bool hasBody1OnlyC21 =
            std::fabs( cosineCoefficientsOfBody1( 2, 1 ) - 1.0 ) < tolerance &&
            std::fabs( cosineCoefficientsOfBody1( 2, 0 ) ) < tolerance &&
            std::fabs( cosineCoefficientsOfBody1( 2, 2 ) ) < tolerance &&
            std::fabs( sineCoefficientsOfBody1( 2, 1 ) ) < tolerance &&
            std::fabs( sineCoefficientsOfBody1( 2, 2 ) ) < tolerance;
    const bool hasBody2OnlyC20 =
            std::fabs( cosineCoefficientsOfBody2( 2, 0 ) - 1.0 ) < tolerance &&
            std::fabs( cosineCoefficientsOfBody2( 2, 1 ) ) < tolerance &&
            std::fabs( cosineCoefficientsOfBody2( 2, 2 ) ) < tolerance &&
            std::fabs( sineCoefficientsOfBody2( 2, 1 ) ) < tolerance &&
            std::fabs( sineCoefficientsOfBody2( 2, 2 ) ) < tolerance;
    return isRequestedRelativePositionOnXAxis && hasBody1OnlyC21 && hasBody2OnlyC20;
}

struct DegreeTwoCoefficientDescriptor
{
    bool isValid = false;
    bool isCosine = true;
    int order = -1;
    std::string name = "";
    double value = 0.0;
};

DegreeTwoCoefficientDescriptor getSingleActiveDegreeTwoCoefficient(
        const Eigen::MatrixXd& cosineCoefficients,
        const Eigen::MatrixXd& sineCoefficients )
{
    DegreeTwoCoefficientDescriptor descriptor;
    if( cosineCoefficients.rows( ) < 3 || cosineCoefficients.cols( ) < 3 ||
        sineCoefficients.rows( ) < 3 || sineCoefficients.cols( ) < 3 )
    {
        return descriptor;
    }

    const double tolerance = 1.0E-14;
    int activeCount = 0;
    auto setIfActive = [ & ]( const double value, const bool isCosine, const int order, const std::string& name )
    {
        if( std::fabs( value ) > tolerance )
        {
            activeCount++;
            descriptor.isValid = true;
            descriptor.isCosine = isCosine;
            descriptor.order = order;
            descriptor.name = name;
            descriptor.value = value;
        }
    };

    setIfActive( cosineCoefficients( 2, 0 ), true, 0, "C20" );
    setIfActive( cosineCoefficients( 2, 1 ), true, 1, "C21" );
    setIfActive( sineCoefficients( 2, 1 ), false, 1, "S21" );
    setIfActive( cosineCoefficients( 2, 2 ), true, 2, "C22" );
    setIfActive( sineCoefficients( 2, 2 ), false, 2, "S22" );

    if( activeCount != 1 )
    {
        descriptor = DegreeTwoCoefficientDescriptor( );
    }
    return descriptor;
}

bool isSingleDegreeTwoByDegreeTwoCombination(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
                coefficientCombinationsToUse )
{
    return coefficientCombinationsToUse.size( ) == 1 &&
            std::get< 0 >( coefficientCombinationsToUse.at( 0 ) ) == 2 &&
            std::get< 2 >( coefficientCombinationsToUse.at( 0 ) ) == 2;
}

Eigen::Vector3d computeExpectedJHatFromDerivationDocument( const Eigen::Vector3d& relativePosition )
{
    const double x = relativePosition.x( );
    const double y = relativePosition.y( );
    const double z = relativePosition.z( );
    const double r2 = relativePosition.squaredNorm( );
    const double r = std::sqrt( r2 );
    const double x2 = x * x;
    const double y2 = y * y;
    const double z2 = z * z;
    const double x4 = x2 * x2;
    const double y4 = y2 * y2;
    const double z4 = z2 * z2;

    const double prefactor = 10.0 * std::sqrt( 3.0 ) / std::pow( r, 9.0 );
    Eigen::Vector3d expectedJHat = Eigen::Vector3d::Zero( );
    expectedJHat << -10.0 * x * y * ( r2 - 7.0 * z2 ),
            4.0 * x4 + 3.0 * x2 * y2 - 27.0 * x2 * z2 - y4 + 3.0 * y2 * z2 + 4.0 * z4,
            0.0;
    return prefactor * expectedJHat;
}

//! Apply the angular momentum operator to one Wigner-D coefficient entry.
/*!
 * Implements the Cartesian form of \hat{J}(D^l_{m,k}) used in the torque evaluation
 * of Dirkx et al. (2019), Eq. (67), with torque definition from Eq. (60).
 */
Eigen::Vector3cd computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices(
        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
        const int degree,
        const int originalOrder,
        const int newOrder )
{
    const std::complex< double > imaginaryUnit( 0.0, 1.0 );
    const double inverseSquareRootTwo = 1.0 / std::sqrt( 2.0 );

    const auto getWignerCoefficient = [ & ]( const int requestedOrderM, const int requestedOrderK )
    {
        if( std::abs( requestedOrderM ) > degree || std::abs( requestedOrderK ) > degree )
        {
            return std::complex< double >( 0.0, 0.0 );
        }
        return wignerMatrices.at( degree )( requestedOrderM + degree, requestedOrderK + degree );
    };

    const double plusScaling = std::sqrt(
            std::max( 0.0, static_cast< double >(
                    degree * ( degree + 1 ) - originalOrder * ( originalOrder - 1 ) ) ) ) / std::sqrt( 2.0 );
    const double minusScaling = std::sqrt(
            std::max( 0.0, static_cast< double >(
                    degree * ( degree + 1 ) - originalOrder * ( originalOrder + 1 ) ) ) ) / std::sqrt( 2.0 );

    // Eq. (67): Cartesian components of \hat{J}(D^l_{m,k}) assembled from ladder-operator contributions.
    const std::complex< double > angularMomentumPlus =
            imaginaryUnit * plusScaling * getWignerCoefficient( originalOrder - 1, newOrder );
    const std::complex< double > angularMomentumMinus =
            imaginaryUnit * ( -minusScaling ) * getWignerCoefficient( originalOrder + 1, newOrder );
    const std::complex< double > angularMomentumZero =
            imaginaryUnit * static_cast< double >( -originalOrder ) * getWignerCoefficient( originalOrder, newOrder );

    Eigen::Vector3cd angularMomentumInCartesianBasis;
    angularMomentumInCartesianBasis( 0 ) = ( angularMomentumMinus - angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 1 ) =
            imaginaryUnit * ( angularMomentumMinus + angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 2 ) = angularMomentumZero;
    // Returned \hat{J}-mapped term is used directly in the Eq. (60) torque evaluation.
    return angularMomentumInCartesianBasis;
}

}  // namespace

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getBody2TorqueCombinationsToUse(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
                coefficientCombinationsToUse )
{
    // Build unique combinations required for body-2 spin torque accumulation (Eq. (67) term set).
    std::set< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > body2TorqueCombinationSet;
    for( const auto& combination : coefficientCombinationsToUse )
    {
        const unsigned int degreeOfBody1ForCombination = std::get< 0 >( combination );
        const unsigned int orderOfBody1ForCombination = std::get< 1 >( combination );
        const unsigned int degreeOfBody2ForCombination = std::get< 2 >( combination );
        for( unsigned int orderOfBody2ForTorque = 0; orderOfBody2ForTorque <= degreeOfBody2ForCombination;
             orderOfBody2ForTorque++ )
        {
            body2TorqueCombinationSet.insert( std::make_tuple(
                    degreeOfBody1ForCombination,
                    orderOfBody1ForCombination,
                    degreeOfBody2ForCombination,
                    orderOfBody2ForTorque ) );
        }
    }

    return std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >(
            body2TorqueCombinationSet.begin( ), body2TorqueCombinationSet.end( ) );
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
    // Eq. (67): computes \hat{J}-mapped coefficient fields; these feed the Eq. (60) torque relation.
    for( int i = 0; i < 3; i++ )
    {
        transformedCosineCoefficientsBody2AngularMomentum.at( i ).setZero(
                    cosineCoefficientsBody2.rows( ), cosineCoefficientsBody2.cols( ) );
        transformedSineCoefficientsBody2AngularMomentum.at( i ).setZero(
                    sineCoefficientsBody2.rows( ), sineCoefficientsBody2.cols( ) );
    }

    for( int degree = 0; degree < cosineCoefficientsBody2.rows( ); degree++ )
    {
        for( int orderM = 0; ( orderM <= degree && orderM < cosineCoefficientsBody2.cols( ) ); orderM++ )
        {
            double orderMMultiplier;
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt(
                            boost::math::factorial< double >( degree - orderM ) /
                            boost::math::factorial< double >( degree + orderM ) );
            }
            else
            {
                orderMMultiplier = ( orderM == 0 ? 1.0 : 1.0 / std::sqrt( 2.0 ) );
            }

            const Eigen::Vector3cd orderZeroAngularMomentumD =
                    computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices(
                            wignerMatrices, degree, orderM, 0 );
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
                    currentMultiplier = std::sqrt(
                                boost::math::factorial< double >( degree + orderK ) /
                                boost::math::factorial< double >( degree - orderK ) ) * orderMMultiplier;
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 ) * orderMMultiplier;
                }

                const double signMultiplier = ( ( ( orderK % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
                const Eigen::Vector3cd positiveOrderAngularMomentumD =
                        computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices(
                                wignerMatrices, degree, orderM, orderK );
                const Eigen::Vector3cd negativeOrderAngularMomentumD =
                        computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices(
                                wignerMatrices, degree, orderM, -orderK );

                for( int i = 0; i < 3; i++ )
                {
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) * cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) -
                                negativeOrderAngularMomentumD( i ).imag( ) ) * sineCoefficientsBody2( degree, orderK ) );

                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) +
                                negativeOrderAngularMomentumD( i ).imag( ) ) * cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) * sineCoefficientsBody2( degree, orderK ) );
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
    // Eq. (67): computes \hat{J}-mapped coefficient fields; these feed the Eq. (60) torque relation.
    for( int i = 0; i < 3; i++ )
    {
        transformedCosineCoefficientsBody2AngularMomentum.at( i ).setZero(
                    cosineCoefficientsBody2.rows( ), cosineCoefficientsBody2.cols( ) );
        transformedSineCoefficientsBody2AngularMomentum.at( i ).setZero(
                    sineCoefficientsBody2.rows( ), sineCoefficientsBody2.cols( ) );
    }

    for( int degree = 0; degree < cosineCoefficientsBody2.rows( ); degree++ )
    {
        for( int orderM = 0; ( orderM <= degree && orderM < cosineCoefficientsBody2.cols( ) ); orderM++ )
        {
            double orderMMultiplier;
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt(
                            boost::math::factorial< double >( degree - orderM ) /
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
                    currentMultiplier = std::sqrt(
                                boost::math::factorial< double >( degree + orderK ) /
                                boost::math::factorial< double >( degree - orderK ) ) * orderMMultiplier;
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
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) * cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) -
                                negativeOrderAngularMomentumD( i ).imag( ) ) * sineCoefficientsBody2( degree, orderK ) );

                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) +
                                negativeOrderAngularMomentumD( i ).imag( ) ) * cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) * sineCoefficientsBody2( degree, orderK ) );
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
    // Eq. (67): this function provides the \hat{J}-mapped coefficients used in Eq. (60).
    computeTransformedAngularMomentumCoefficientsFromWignerCache(
            cosineCoefficientsBody2,
            sineCoefficientsBody2,
            wignerCache,
            areCoefficientsNormalized,
            transformedCosineCoefficientsBody2AngularMomentum,
            transformedSineCoefficientsBody2AngularMomentum );
}

void FullTwoBodySphericalHarmonicTorque::updateMembers( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        // Step 1: synchronize the acceleration model state (effective coefficients, SH cache, relative geometry).
        accelerationBetweenBodies_->updateMembers( currentTime );

        // Step 2: retrieve current fields and precompute \hat{J}-transformed body-2 coefficients
        // used in torque Eq. (67) / Eq. (60) evaluation.
        const std::shared_ptr< EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField =
                accelerationBetweenBodies_->getEffectiveMutualPotentialField( );
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache > wignerCache =
                effectiveMutualPotentialField->getTransformationCache( )->getWignerDMatricesCache( );

        const Eigen::MatrixXd& cosineCoefficientsOfBody1 = effectiveMutualPotentialField->getCosineCoefficientsOfBody1( );
        const Eigen::MatrixXd& sineCoefficientsOfBody1 = effectiveMutualPotentialField->getSineCoefficientsOfBody1( );
        const Eigen::MatrixXd& cosineCoefficientsOfBody2 = effectiveMutualPotentialField->getCosineCoefficientsOfBody2( );
        const Eigen::MatrixXd& sineCoefficientsOfBody2 = effectiveMutualPotentialField->getSineCoefficientsOfBody2( );
        computeTransformedAngularMomentumCoefficients(
                    cosineCoefficientsOfBody2,
                    sineCoefficientsOfBody2,
                    wignerCache,
                    accelerationBetweenBodies_->getAreCoefficientsNormalized( ),
                    transformedCosineCoefficientsBody2AngularMomentum_,
                    transformedSineCoefficientsBody2AngularMomentum_ );

        const Eigen::Vector3d bodyFixedRelativePosition = accelerationBetweenBodies_->getCurrentBodyFixedRelativePosition( );
        const DegreeTwoCoefficientDescriptor activeBody1DegreeTwoCoefficient =
                getSingleActiveDegreeTwoCoefficient( cosineCoefficientsOfBody1, sineCoefficientsOfBody1 );
        const DegreeTwoCoefficientDescriptor activeBody2DegreeTwoCoefficient =
                getSingleActiveDegreeTwoCoefficient( cosineCoefficientsOfBody2, sineCoefficientsOfBody2 );
        const bool isExpandedSingleInteractionDebugCase =
                isSingleDegreeTwoByDegreeTwoCombination( coefficientCombinationsToUse_ ) &&
                activeBody1DegreeTwoCoefficient.isValid &&
                activeBody2DegreeTwoCoefficient.isValid &&
                ( activeBody1DegreeTwoCoefficient.name == "C21" ||
                  activeBody1DegreeTwoCoefficient.name == "S21" ||
                  activeBody1DegreeTwoCoefficient.name == "S22" );
        const bool isDebugCase = isUnitC21ByC20DebugCase(
                coefficientCombinationsToUse_,
                bodyFixedRelativePosition,
                cosineCoefficientsOfBody1,
                sineCoefficientsOfBody1,
                cosineCoefficientsOfBody2,
                sineCoefficientsOfBody2 );
        if( isDebugCase )
        {
            const double sqrtThree = std::sqrt( 3.0 );
            std::cout << "[DBG Eq66/67 transformed J*C/J*S for body2 degree-2]" << std::endl;
            std::cout << "  m=0 JC computed="
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 0 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 0 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 0 )
                      << " expected=0 0 0" << std::endl;
            std::cout << "  m=0 JS computed="
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 0 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 0 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 0 )
                      << " expected=0 0 0" << std::endl;
            std::cout << "  m=1 JC computed="
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 1 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 1 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 1 )
                      << " expected=0 " << sqrtThree << " 0" << std::endl;
            std::cout << "  m=1 JS computed="
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 1 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 1 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 1 )
                      << " expected=" << sqrtThree << " 0 0" << std::endl;
            std::cout << "  m=2 JC computed="
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 2 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 2 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 2 )
                      << " expected=0 0 0" << std::endl;
            std::cout << "  m=2 JS computed="
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 2 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 2 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 2 )
                      << " expected=0 0 0" << std::endl;
        }
        const double currentDistance = bodyFixedRelativePosition.norm( );
        const double preMultiplier = accelerationBetweenBodies_->getCurrentGravitationalParameter( ) / currentDistance;
        if( isDebugCase )
        {
            std::cout << "[DBG state] r=" << currentDistance
                      << " preMultiplier(GM/r)=" << preMultiplier
                      << " mu=" << accelerationBetweenBodies_->getCurrentGravitationalParameter( ) << std::endl;
        }

        const std::vector< double >& radius1Powers = accelerationBetweenBodies_->getRadius1Powers( );
        const std::vector< double >& radius2Powers = accelerationBetweenBodies_->getRadius2Powers( );
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
                accelerationBetweenBodies_->getSphericalHarmonicsCache( );

        const Eigen::Vector3d expectedJHatFromDerivation =
                computeExpectedJHatFromDerivationDocument( bodyFixedRelativePosition );
        Eigen::Vector3d debugEq67AlternativeSourceOrder0 = Eigen::Vector3d::Zero( );

        // Step 3: evaluate body-2 spin torque contribution in frame F1 from Eq. (67), then apply Eq. (60).
        Eigen::Vector3d body2TorqueInBodyFixedFrameOfBody1 = Eigen::Vector3d::Zero( );
        for( unsigned int i = 0; i < body2TorqueCombinationsToUse_.size( ); i++ )
        {
            const int degreeOfBody1 = std::get< 0 >( body2TorqueCombinationsToUse_.at( i ) );
            const int orderOfBody1 = std::get< 1 >( body2TorqueCombinationsToUse_.at( i ) );
            const int degreeOfBody2 = std::get< 2 >( body2TorqueCombinationsToUse_.at( i ) );
            const int orderOfBody2 = std::get< 3 >( body2TorqueCombinationsToUse_.at( i ) );

            const double equatorialRadiusRatioPower =
                    radius1Powers.at( degreeOfBody1 ) * radius2Powers.at( degreeOfBody2 );
            const int totalDegree = degreeOfBody1 + degreeOfBody2;

            // Expand each stored non-negative-order tuple into all signed-order combinations
            // consistent with Eq. (49)/(67) real-coefficient convention.
            for( int j = 0; j < 4; j++ )
            {
                int signedOrderOfBody1 = 0;
                int signedOrderOfBody2 = 0;
                const bool computeTerm = getSignedOrdersForCombinationCase(
                        j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 );
                if( !computeTerm )
                {
                    continue;
                }

                const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                const double legendrePolynomial = sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                            totalDegree, totalOrder );
                const double cosineOfMultipleLongitude = sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder );
                const double sineOfMultipleLongitude = sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder );

                const double signOrderBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
                const double signOrderBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
                const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;

                const double body1CosineCoefficient = cosineCoefficientsOfBody1(
                            degreeOfBody1, std::abs( signedOrderOfBody1 ) );
                const double body1SineCoefficient = sineCoefficientsOfBody1(
                            degreeOfBody1, std::abs( signedOrderOfBody1 ) );
                const double multiplier = getMutualPotentialEffectiveCoefficientMultiplier(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2,
                            accelerationBetweenBodies_->getAreCoefficientsNormalized( ) );

                Eigen::Vector3d angularMomentumTransformedCosineCoefficientsBody2;
                Eigen::Vector3d angularMomentumTransformedSineCoefficientsBody2;
                for( int k = 0; k < 3; k++ )
                {
                    angularMomentumTransformedCosineCoefficientsBody2( k ) =
                            transformedCosineCoefficientsBody2AngularMomentum_.at( k )(
                                degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                    angularMomentumTransformedSineCoefficientsBody2( k ) =
                            transformedSineCoefficientsBody2AngularMomentum_.at( k )(
                                degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                }

                const Eigen::Vector3d effectiveAngularMomentumCosineCoefficients =
                        ( body1CosineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 -
                          signOrderBody1 * signOrderBody2 * body1SineCoefficient * angularMomentumTransformedSineCoefficientsBody2 ) *
                        multiplier;
                const Eigen::Vector3d effectiveAngularMomentumSineCoefficients =
                        ( signOrderBody2 * body1CosineCoefficient * angularMomentumTransformedSineCoefficientsBody2 +
                          signOrderBody1 * body1SineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 ) *
                        signTotalOrder * multiplier;

                // Scalar basis and effective angular-momentum coefficients form one Eq. (67) term.
                const Eigen::Vector3d currentEq67Contribution =
                        equatorialRadiusRatioPower * legendrePolynomial *
                        ( effectiveAngularMomentumCosineCoefficients * cosineOfMultipleLongitude +
                          effectiveAngularMomentumSineCoefficients * sineOfMultipleLongitude );
                body2TorqueInBodyFixedFrameOfBody1 += currentEq67Contribution;

                Eigen::Vector3d currentEq67ContributionSourceOrder0 = Eigen::Vector3d::Zero( );
                if( isDebugCase && coefficientCombinationsToUse_.size( ) == 1 )
                {
                    const int selectedOrderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( 0 ) );
                    const int signedSelectedOrderOfBody2 = ( selectedOrderOfBody2 == 0 ) ? 0 : signedOrderOfBody2;
                    const double signSelectedOrderBody2 = ( signedSelectedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
                    const double signSelectedTotalOrder =
                            ( ( signedOrderOfBody1 + signedSelectedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;
                    const int selectedTotalOrder = std::abs( signedOrderOfBody1 + signedSelectedOrderOfBody2 );
                    const double selectedLegendrePolynomial =
                            sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial( totalDegree, selectedTotalOrder );
                    const double selectedCosineOfMultipleLongitude =
                            sphericalHarmonicsCache->getCosineOfMultipleLongitude( selectedTotalOrder );
                    const double selectedSineOfMultipleLongitude =
                            sphericalHarmonicsCache->getSineOfMultipleLongitude( selectedTotalOrder );
                    const double selectedOrderMultiplier = getMutualPotentialEffectiveCoefficientMultiplier(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedSelectedOrderOfBody2,
                            accelerationBetweenBodies_->getAreCoefficientsNormalized( ) );
                    const Eigen::Vector3d effectiveAngularMomentumCosineCoefficientsSourceOrder0 =
                            ( body1CosineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 -
                              signOrderBody1 * signSelectedOrderBody2 * body1SineCoefficient *
                              angularMomentumTransformedSineCoefficientsBody2 ) * selectedOrderMultiplier;
                    const Eigen::Vector3d effectiveAngularMomentumSineCoefficientsSourceOrder0 =
                            ( signSelectedOrderBody2 * body1CosineCoefficient * angularMomentumTransformedSineCoefficientsBody2 +
                              signOrderBody1 * body1SineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 ) *
                            signSelectedTotalOrder * selectedOrderMultiplier;
                    currentEq67ContributionSourceOrder0 =
                            equatorialRadiusRatioPower * selectedLegendrePolynomial *
                            ( effectiveAngularMomentumCosineCoefficientsSourceOrder0 * selectedCosineOfMultipleLongitude +
                              effectiveAngularMomentumSineCoefficientsSourceOrder0 * selectedSineOfMultipleLongitude );
                    debugEq67AlternativeSourceOrder0 += currentEq67ContributionSourceOrder0;
                }

                if( isDebugCase )
                {
                    std::cout << "[DBG Eq67 term] (l1,m1,l2,m2)=("
                              << degreeOfBody1 << "," << signedOrderOfBody1 << ","
                              << degreeOfBody2 << "," << signedOrderOfBody2 << ") "
                              << "term_computed=" << currentEq67Contribution.transpose( )
                              << " term_alt_source_order0=" << currentEq67ContributionSourceOrder0.transpose( )
                              << " running_sum_computed=" << body2TorqueInBodyFixedFrameOfBody1.transpose( )
                              << " running_sum_expected_final="
                              << expectedJHatFromDerivation.transpose( ) << std::endl;
                    std::cout << "  [DBG Eq67 term factors] P_lm=" << legendrePolynomial
                              << " cos(m*lon)=" << cosineOfMultipleLongitude
                              << " sin(m*lon)=" << sineOfMultipleLongitude
                              << " multiplier=" << multiplier
                              << " body1(C,S)=(" << body1CosineCoefficient << "," << body1SineCoefficient << ")"
                              << " JC=" << angularMomentumTransformedCosineCoefficientsBody2.transpose( )
                              << " JS=" << angularMomentumTransformedSineCoefficientsBody2.transpose( ) << std::endl;
                }

            }
        }
        if( isDebugCase )
        {
            std::cout << "[DBG Eq67 sum Jhat(V1-2)] computed="
                      << body2TorqueInBodyFixedFrameOfBody1.transpose( )
                      << " alt_source_order0=" << debugEq67AlternativeSourceOrder0.transpose( )
                      << " expected=" << expectedJHatFromDerivation.transpose( ) << std::endl;
        }
        const Eigen::Vector3d body2TorqueEq67SumBeforePremult = body2TorqueInBodyFixedFrameOfBody1;
        body2TorqueInBodyFixedFrameOfBody1 *= -preMultiplier;
        // Eq. (60): M_2 = -\hat{J}(V_1-2), with preMultiplier carrying the common -GM/r factor.
        if( isDebugCase )
        {
            const Eigen::Vector3d expectedBody2TorqueFromDerivation = -preMultiplier * expectedJHatFromDerivation;
            std::cout << "[DBG Eq60 M2=-Jhat(V1-2)] computed="
                      << body2TorqueInBodyFixedFrameOfBody1.transpose( )
                      << " expected=" << expectedBody2TorqueFromDerivation.transpose( ) << std::endl;
        }

        // Step 4: compute total relative-frame torque from translational side using Eq. (68),
        // then isolate body-1 torque by subtracting body-2 contribution.
        const Eigen::Vector3d totalTorqueInBodyFixedFrameOfBody1 =
                bodyFixedRelativePosition.cross( accelerationBetweenBodies_->getMutualPotentialGradient( ) );
        const Eigen::Vector3d body1TorqueInBodyFixedFrameOfBody1 =
                totalTorqueInBodyFixedFrameOfBody1 - body2TorqueInBodyFixedFrameOfBody1;
        if( isExpandedSingleInteractionDebugCase )
        {
            std::cout << "[DBG degree2 Eq67/Eq68 breakdown] case="
                      << activeBody1DegreeTwoCoefficient.name << "x"
                      << activeBody2DegreeTwoCoefficient.name
                      << " r=" << bodyFixedRelativePosition.transpose( ) << std::endl;
            std::cout << "  Eq67_sum_before_premult="
                      << body2TorqueEq67SumBeforePremult.transpose( )
                      << " preMultiplier(GM/r)=" << preMultiplier << std::endl;
            std::cout << "  body2_torque_from_Eq60="
                      << body2TorqueInBodyFixedFrameOfBody1.transpose( ) << std::endl;
            std::cout << "  total_torque_from_Eq68=rxdUdr="
                      << totalTorqueInBodyFixedFrameOfBody1.transpose( ) << std::endl;
            std::cout << "  body1_torque_from_Eq68_minus_Eq60="
                      << body1TorqueInBodyFixedFrameOfBody1.transpose( ) << std::endl;
            std::cout << "  transformed_JC(m=0,1,2)_x="
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 0 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 1 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 2 ) << std::endl;
            std::cout << "  transformed_JC(m=0,1,2)_y="
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 0 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 1 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 2 ) << std::endl;
            std::cout << "  transformed_JC(m=0,1,2)_z="
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 0 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 1 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 2 ) << std::endl;
            std::cout << "  transformed_JS(m=0,1,2)_x="
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 0 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 1 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 2 ) << std::endl;
            std::cout << "  transformed_JS(m=0,1,2)_y="
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 0 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 1 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 2 ) << std::endl;
            std::cout << "  transformed_JS(m=0,1,2)_z="
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 0 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 1 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 2 ) << std::endl;
        }
        if( isDebugCase )
        {
            // For this exact debug setup (unit C21 x C20, r=[+/-1,0,0]), Eq. (8) analytical body-1 torque
            // in the test-side normalized convention is [0, -35*sqrt(3)/2, 0].
            const Eigen::Vector3d expectedCurrentTorqueFromDerivation(
                    0.0, -35.0 * std::sqrt( 3.0 ) / 2.0, 0.0 );
            const Eigen::Vector3d expectedBody1TorqueFromDerivation = -expectedCurrentTorqueFromDerivation;
            const Eigen::Vector3d expectedBody2TorqueFromDerivationClosure =
                    totalTorqueInBodyFixedFrameOfBody1 - expectedBody1TorqueFromDerivation;
            const Eigen::Vector3d expectedEq67SumFromDerivationClosure =
                    -expectedBody2TorqueFromDerivationClosure / preMultiplier;
            const double requiredEq67ScalingFromClosure =
                    ( std::fabs( body2TorqueEq67SumBeforePremult.y( ) ) > 1.0E-15 ) ?
                    ( expectedEq67SumFromDerivationClosure.y( ) / body2TorqueEq67SumBeforePremult.y( ) ) : TUDAT_NAN;
            const double computedJCyDegreeTwoOrderOne =
                    transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 1 );
            const double expectedJCyFromClosure =
                    computedJCyDegreeTwoOrderOne * requiredEq67ScalingFromClosure;

            const Eigen::Vector3d expectedTotalTorqueFromEquation68 =
                    body1TorqueInBodyFixedFrameOfBody1 + body2TorqueInBodyFixedFrameOfBody1;
            std::cout << "[DBG Eq68 torque balance] r x dU/dr="
                      << totalTorqueInBodyFixedFrameOfBody1.transpose( )
                      << " expected_eq67_from_closure="
                      << expectedEq67SumFromDerivationClosure.transpose( )
                      << " eq67_required_scaling_from_closure=" << requiredEq67ScalingFromClosure
                      << " JC_y_computed=" << computedJCyDegreeTwoOrderOne
                      << " JC_y_implied_from_closure=" << expectedJCyFromClosure
                      << " expected(body1+body2)="
                      << expectedTotalTorqueFromEquation68.transpose( ) << std::endl;

        }

        // Step 5: return requested body's torque, applying frame mapping for body 2 using Eq. (69).
        if( acceleratedBodyIsBody1_ )
        {
            currentTorque_ = -body1TorqueInBodyFixedFrameOfBody1;
        }
        else
        {
            currentTorque_ = -( accelerationBetweenBodies_->getCurrentRotationFromBody2ToBody1( ).inverse( ) *
                    body2TorqueInBodyFixedFrameOfBody1 );
        }

        currentTime_ = currentTime;
    }
}

}

}
