#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"

#include <boost/math/special_functions/factorials.hpp>

#include <algorithm>
#include <cmath>

namespace tudat
{

namespace gravitation
{

Eigen::Vector3cd FullTwoBodySphericalHarmonicTorque::computeAngularMomentumOperatorOnWignerCoefficient(
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const int degree,
        const int orderM,
        const int orderK )
{
    const std::complex< double > imaginaryUnit( 0.0, 1.0 );
    const double inverseSquareRootTwo = 1.0 / std::sqrt( 2.0 );

    const auto getWignerCoefficient = [ & ]( const int requestedOrderM, const int requestedOrderK )
    {
        if( std::abs( requestedOrderM ) > degree || std::abs( requestedOrderK ) > degree )
        {
            return std::complex< double >( 0.0, 0.0 );
        }
        return wignerCache->getWignerDCoefficient( degree, requestedOrderM, requestedOrderK );
    };

    const double plusScaling = std::sqrt(
                std::max( 0.0, static_cast< double >(
                              degree * ( degree + 1 ) - orderK * ( orderK - 1 ) ) ) ) / 2.0;
    const double minusScaling = std::sqrt(
                std::max( 0.0, static_cast< double >(
                              degree * ( degree + 1 ) - orderK * ( orderK + 1 ) ) ) ) / 2.0;

    const std::complex< double > angularMomentumPlus =
            imaginaryUnit * plusScaling * getWignerCoefficient( orderM, orderK - 1 );
    const std::complex< double > angularMomentumMinus =
            imaginaryUnit * ( -minusScaling ) * getWignerCoefficient( orderM, orderK + 1 );
    const std::complex< double > angularMomentumZero =
            imaginaryUnit * static_cast< double >( -orderK ) * getWignerCoefficient( orderM, orderK );

    Eigen::Vector3cd angularMomentumInCartesianBasis;
    angularMomentumInCartesianBasis( 0 ) = ( angularMomentumMinus - angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 1 ) = imaginaryUnit * ( angularMomentumMinus + angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 2 ) = angularMomentumZero;
    return angularMomentumInCartesianBasis;
}

void FullTwoBodySphericalHarmonicTorque::computeTransformedAngularMomentumCoefficients(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum )
{
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

            const Eigen::Vector3cd orderZeroAngularMomentumD = computeAngularMomentumOperatorOnWignerCoefficient(
                        wignerCache, degree, orderM, 0 );
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
                const Eigen::Vector3cd positiveOrderAngularMomentumD = computeAngularMomentumOperatorOnWignerCoefficient(
                            wignerCache, degree, orderM, orderK );
                const Eigen::Vector3cd negativeOrderAngularMomentumD = computeAngularMomentumOperatorOnWignerCoefficient(
                            wignerCache, degree, orderM, -orderK );

                for( int i = 0; i < 3; i++ )
                {
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
                              cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) - negativeOrderAngularMomentumD( i ).imag( ) ) *
                              sineCoefficientsBody2( degree, orderK ) );

                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) + negativeOrderAngularMomentumD( i ).imag( ) ) *
                              cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
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

void FullTwoBodySphericalHarmonicTorque::updateMembers( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        accelerationBetweenBodies_->updateMembers( currentTime );

        const std::shared_ptr< EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField =
                accelerationBetweenBodies_->getEffectiveMutualPotentialField( );
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache > wignerCache =
                effectiveMutualPotentialField->getTransformationCache( )->getWignerDMatricesCache( );

        const Eigen::MatrixXd cosineCoefficientsOfBody1 = effectiveMutualPotentialField->getCosineCoefficientsOfBody1( );
        const Eigen::MatrixXd sineCoefficientsOfBody1 = effectiveMutualPotentialField->getSineCoefficientsOfBody1( );
        const Eigen::MatrixXd cosineCoefficientsOfBody2 = effectiveMutualPotentialField->getCosineCoefficientsOfBody2( );
        const Eigen::MatrixXd sineCoefficientsOfBody2 = effectiveMutualPotentialField->getSineCoefficientsOfBody2( );

        std::array< Eigen::MatrixXd, 3 > transformedCosineCoefficientsBody2AngularMomentum;
        std::array< Eigen::MatrixXd, 3 > transformedSineCoefficientsBody2AngularMomentum;
        computeTransformedAngularMomentumCoefficients(
                    cosineCoefficientsOfBody2,
                    sineCoefficientsOfBody2,
                    wignerCache,
                    accelerationBetweenBodies_->getAreCoefficientsNormalized( ),
                    transformedCosineCoefficientsBody2AngularMomentum,
                    transformedSineCoefficientsBody2AngularMomentum );

        const Eigen::Vector3d bodyFixedRelativePosition = accelerationBetweenBodies_->getCurrentBodyFixedRelativePosition( );
        const double currentDistance = bodyFixedRelativePosition.norm( );
        const double preMultiplier = accelerationBetweenBodies_->getCurrentGravitationalParameter( ) / currentDistance;

        const std::vector< double > radius1Powers = accelerationBetweenBodies_->getRadius1Powers( );
        const std::vector< double > radius2Powers = accelerationBetweenBodies_->getRadius2Powers( );
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
                accelerationBetweenBodies_->getSphericalHarmonicsCache( );

        Eigen::Vector3d body2TorqueInBodyFixedFrameOfBody1 = Eigen::Vector3d::Zero( );
        for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
        {
            const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
            const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );

            const double equatorialRadiusRatioPower =
                    radius1Powers.at( degreeOfBody1 ) * radius2Powers.at( degreeOfBody2 );
            const int totalDegree = degreeOfBody1 + degreeOfBody2;

            for( unsigned int j = 0; j < 4; j++ )
            {
                int signedOrderOfBody1 = orderOfBody1;
                int signedOrderOfBody2 = orderOfBody2;
                bool computeTerm = true;

                switch( j )
                {
                case 0:
                    signedOrderOfBody1 = orderOfBody1;
                    signedOrderOfBody2 = orderOfBody2;
                    break;
                case 1:
                    signedOrderOfBody1 = -orderOfBody1;
                    signedOrderOfBody2 = orderOfBody2;
                    computeTerm = ( orderOfBody1 != 0 );
                    break;
                case 2:
                    signedOrderOfBody1 = orderOfBody1;
                    signedOrderOfBody2 = -orderOfBody2;
                    computeTerm = ( orderOfBody2 != 0 );
                    break;
                case 3:
                    signedOrderOfBody1 = -orderOfBody1;
                    signedOrderOfBody2 = -orderOfBody2;
                    computeTerm = ( orderOfBody1 != 0 && orderOfBody2 != 0 );
                    break;
                default:
                    break;
                }

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
                const double multiplier = effectiveMutualPotentialField->getMultiplier(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

                Eigen::Vector3d angularMomentumTransformedCosineCoefficientsBody2;
                Eigen::Vector3d angularMomentumTransformedSineCoefficientsBody2;
                for( int k = 0; k < 3; k++ )
                {
                    angularMomentumTransformedCosineCoefficientsBody2( k ) =
                            transformedCosineCoefficientsBody2AngularMomentum.at( k )( degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                    angularMomentumTransformedSineCoefficientsBody2( k ) =
                            transformedSineCoefficientsBody2AngularMomentum.at( k )( degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                }

                const Eigen::Vector3d effectiveAngularMomentumCosineCoefficients =
                        ( body1CosineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 -
                          signOrderBody1 * signOrderBody2 * body1SineCoefficient * angularMomentumTransformedSineCoefficientsBody2 ) *
                        multiplier;
                const Eigen::Vector3d effectiveAngularMomentumSineCoefficients =
                        ( signOrderBody2 * body1CosineCoefficient * angularMomentumTransformedSineCoefficientsBody2 +
                          signOrderBody1 * body1SineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 ) *
                        signTotalOrder * multiplier;

                body2TorqueInBodyFixedFrameOfBody1 +=
                        equatorialRadiusRatioPower * legendrePolynomial *
                        ( effectiveAngularMomentumCosineCoefficients * cosineOfMultipleLongitude +
                          effectiveAngularMomentumSineCoefficients * sineOfMultipleLongitude );
            }
        }
        body2TorqueInBodyFixedFrameOfBody1 *= preMultiplier;

        const Eigen::Vector3d totalTorqueInBodyFixedFrameOfBody1 =
                bodyFixedRelativePosition.cross( accelerationBetweenBodies_->getMutualPotentialGradient( ) );
        const Eigen::Vector3d body1TorqueInBodyFixedFrameOfBody1 =
                totalTorqueInBodyFixedFrameOfBody1 - body2TorqueInBodyFixedFrameOfBody1;

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
