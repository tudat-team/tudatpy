#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"

#include <boost/math/special_functions/factorials.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
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
                    degree * ( degree + 1 ) - originalOrder * ( originalOrder - 1 ) ) ) ) / 2.0;
    const double minusScaling = std::sqrt(
            std::max( 0.0, static_cast< double >(
                    degree * ( degree + 1 ) - originalOrder * ( originalOrder + 1 ) ) ) ) / 2.0;

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

template< typename AngularMomentumOperatorGetter >
//! Shared implementation of transformed angular-momentum coefficient evaluation.
/*!
 * Computes \hat{J}\bar{C}^{2,F1}_{l,m} and \hat{J}\bar{S}^{2,F1}_{l,m} from body-2 coefficients and an
 * angular-momentum-on-Wigner accessor. The resulting fields are the coefficient part of Eq. (67).
 */
void computeTransformedAngularMomentumCoefficientsImpl(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        AngularMomentumOperatorGetter getAngularMomentumOperator,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum )
{
    // Initialize output fields for the 3 Cartesian components of \hat{J}.
    for( int i = 0; i < 3; i++ )
    {
        transformedCosineCoefficientsBody2AngularMomentum.at( i ).setZero(
                    cosineCoefficientsBody2.rows( ), cosineCoefficientsBody2.cols( ) );
        transformedSineCoefficientsBody2AngularMomentum.at( i ).setZero(
                    sineCoefficientsBody2.rows( ), sineCoefficientsBody2.cols( ) );
    }

    // Loop over (l,m) entries of the transformed coefficient field.
    for( int degree = 0; degree < cosineCoefficientsBody2.rows( ); degree++ )
    {
        for( int orderM = 0; ( orderM <= degree && orderM < cosineCoefficientsBody2.cols( ) ); orderM++ )
        {
            // Normalization-dependent prefactor relating real SH coefficients to Wigner-D terms.
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

            // Contribution from k=0.
            const Eigen::Vector3cd orderZeroAngularMomentumD = getAngularMomentumOperator( degree, orderM, 0 );
            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).real( ) * cosineCoefficientsBody2( degree, 0 );
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).imag( ) * cosineCoefficientsBody2( degree, 0 );
            }

            // Contributions from k>0, combining +/-k real-coefficient terms.
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
                const Eigen::Vector3cd positiveOrderAngularMomentumD = getAngularMomentumOperator( degree, orderM, orderK );
                const Eigen::Vector3cd negativeOrderAngularMomentumD = getAngularMomentumOperator( degree, orderM, -orderK );

                for( int i = 0; i < 3; i++ )
                {
                    // Eq. (67): real-coefficient combination for \hat{J}\bar{C}^{2,F1}_{l,m}.
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
                              cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) - negativeOrderAngularMomentumD( i ).imag( ) ) *
                              sineCoefficientsBody2( degree, orderK ) );

                    // Eq. (67): real-coefficient combination for \hat{J}\bar{S}^{2,F1}_{l,m}.
                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) + negativeOrderAngularMomentumD( i ).imag( ) ) *
                              cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
                              sineCoefficientsBody2( degree, orderK ) );
                }
            }

            // Final scaling to convert to real cosine/sine coefficient convention.
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
    computeTransformedAngularMomentumCoefficientsImpl(
            cosineCoefficientsBody2,
            sineCoefficientsBody2,
            std::bind(
                    &computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices,
                    std::cref( wignerMatrices ),
                    std::placeholders::_1,
                    std::placeholders::_2,
                    std::placeholders::_3 ),
            areCoefficientsNormalized,
            transformedCosineCoefficientsBody2AngularMomentum,
            transformedSineCoefficientsBody2AngularMomentum );
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
    computeTransformedAngularMomentumCoefficientsImpl(
            cosineCoefficientsBody2,
            sineCoefficientsBody2,
            std::bind(
                    &basic_mathematics::WignerDMatricesCache::getAngularMomentumOperatorOnWignerCoefficient,
                    wignerCache,
                    std::placeholders::_1,
                    std::placeholders::_2,
                    std::placeholders::_3 ),
            areCoefficientsNormalized,
            transformedCosineCoefficientsBody2AngularMomentum,
            transformedSineCoefficientsBody2AngularMomentum );
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
        const double currentDistance = bodyFixedRelativePosition.norm( );
        const double preMultiplier = accelerationBetweenBodies_->getCurrentGravitationalParameter( ) / currentDistance;

        const std::vector< double >& radius1Powers = accelerationBetweenBodies_->getRadius1Powers( );
        const std::vector< double >& radius2Powers = accelerationBetweenBodies_->getRadius2Powers( );
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
                accelerationBetweenBodies_->getSphericalHarmonicsCache( );

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
                body2TorqueInBodyFixedFrameOfBody1 +=
                        equatorialRadiusRatioPower * legendrePolynomial *
                        ( effectiveAngularMomentumCosineCoefficients * cosineOfMultipleLongitude +
                          effectiveAngularMomentumSineCoefficients * sineOfMultipleLongitude );

            }
        }
        body2TorqueInBodyFixedFrameOfBody1 *= -preMultiplier;
        // Eq. (60): M_2 = -\hat{J}(V_1-2), with preMultiplier carrying the common -GM/r factor.

        // Step 4: compute total relative-frame torque from translational side using Eq. (68),
        // then isolate body-1 torque by subtracting body-2 contribution.
        const Eigen::Vector3d totalTorqueInBodyFixedFrameOfBody1 =
                bodyFixedRelativePosition.cross( accelerationBetweenBodies_->getMutualPotentialGradient( ) );
        const Eigen::Vector3d body1TorqueInBodyFixedFrameOfBody1 =
                totalTorqueInBodyFixedFrameOfBody1 - body2TorqueInBodyFixedFrameOfBody1;

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
