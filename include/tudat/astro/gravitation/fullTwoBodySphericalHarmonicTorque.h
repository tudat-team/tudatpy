#ifndef MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
#define MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP

#include <array>
#include <complex>
#include <tuple>
#include <memory>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicAcceleration.h"
#include "tudat/math/basic/wignerDMatrices.h"

namespace tudat
{

namespace gravitation
{

//! Build combination list used for body-2 spin torque evaluation.
/*!
 * Expands each selected (l1,m1,l2,m2) interaction tuple to all body-2 orders needed by the
 * angular-momentum-operator torque summation (Dirkx et al. (2019), Eq. (67)).
 */
std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getBody2TorqueCombinationsToUse(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse );

//! Compute transformed angular-momentum coefficient fields from explicit Wigner D-matrices.
/*!
 * Evaluates \hat{J}\bar{C}_{l,m}^{2,F1} and \hat{J}\bar{S}_{l,m}^{2,F1} needed for the body-2 torque
 * contribution in Dirkx et al. (2019), Eqs. (60) and (67).
 */
void computeTransformedAngularMomentumCoefficientsFromWignerMatrices(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum );

//! Compute transformed angular-momentum coefficient fields from a Wigner D cache.
/*!
 * Cache-based variant of computeTransformedAngularMomentumCoefficientsFromWignerMatrices.
 */
void computeTransformedAngularMomentumCoefficientsFromWignerCache(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum );

//! Full two-body mutual spherical-harmonic torque model (independent of the fourth-degree tensor torque model).
class FullTwoBodySphericalHarmonicTorque : public basic_astrodynamics::TorqueModel
{
public:
    //! Constructor.
    FullTwoBodySphericalHarmonicTorque( const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > accelerationBetweenBodies,
                                        const bool acceleratedBodyIsBody1,
                                        const std::function< double( ) > bodyUndergoingTorqueMassFunction ):
        accelerationBetweenBodies_( accelerationBetweenBodies ), bodyUndergoingTorqueMassFunction_( bodyUndergoingTorqueMassFunction ),
        acceleratedBodyIsBody1_( acceleratedBodyIsBody1 )
    {
        coefficientCombinationsToUse_ = accelerationBetweenBodies_->getEffectiveMutualPotentialField( )->getCoefficientCombinationsToUse( );
        body2TorqueCombinationsToUse_ = gravitation::getBody2TorqueCombinationsToUse( coefficientCombinationsToUse_ );
    }

    //! Update torque to current epoch.
    void updateMembers( const double currentTime = TUDAT_NAN ) override
    {
        if( !( currentTime_ == currentTime ) )
        {
            accelerationBetweenBodies_->updateMembers( currentTime );

            const std::shared_ptr< EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField =
                    accelerationBetweenBodies_->getEffectiveMutualPotentialField( );
            const std::shared_ptr< basic_mathematics::WignerDMatricesCache > wignerCache =
                    effectiveMutualPotentialField->getTransformationCache( )->getWignerDMatricesCache( );

            const Eigen::MatrixXd& cosineCoefficientsOfBody1 = effectiveMutualPotentialField->getCosineCoefficientsOfBody1( );
            const Eigen::MatrixXd& sineCoefficientsOfBody1 = effectiveMutualPotentialField->getSineCoefficientsOfBody1( );
            const Eigen::MatrixXd& cosineCoefficientsOfBody2 = effectiveMutualPotentialField->getCosineCoefficientsOfBody2( );
            const Eigen::MatrixXd& sineCoefficientsOfBody2 = effectiveMutualPotentialField->getSineCoefficientsOfBody2( );
            computeTransformedAngularMomentumCoefficients( cosineCoefficientsOfBody2,
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

            Eigen::Vector3d body2TorqueInBodyFixedFrameOfBody1 = Eigen::Vector3d::Zero( );
            for( unsigned int i = 0; i < body2TorqueCombinationsToUse_.size( ); i++ )
            {
                const int degreeOfBody1 = std::get< 0 >( body2TorqueCombinationsToUse_.at( i ) );
                const int orderOfBody1 = std::get< 1 >( body2TorqueCombinationsToUse_.at( i ) );
                const int degreeOfBody2 = std::get< 2 >( body2TorqueCombinationsToUse_.at( i ) );
                const int orderOfBody2 = std::get< 3 >( body2TorqueCombinationsToUse_.at( i ) );

                const double equatorialRadiusRatioPower = radius1Powers.at( degreeOfBody1 ) * radius2Powers.at( degreeOfBody2 );
                const int totalDegree = degreeOfBody1 + degreeOfBody2;

                for( int j = 0; j < 4; j++ )
                {
                    int signedOrderOfBody1 = 0;
                    int signedOrderOfBody2 = 0;
                    if( !getSignedOrdersForCombinationCase( j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
                    {
                        continue;
                    }

                    const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                    const double legendrePolynomial =
                            sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial( totalDegree, totalOrder );
                    const double cosineOfMultipleLongitude = sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder );
                    const double sineOfMultipleLongitude = sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder );

                    const double signOrderBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
                    const double signOrderBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
                    const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;

                    const double body1CosineCoefficient = cosineCoefficientsOfBody1( degreeOfBody1, std::abs( signedOrderOfBody1 ) );
                    const double body1SineCoefficient = sineCoefficientsOfBody1( degreeOfBody1, std::abs( signedOrderOfBody1 ) );
                    const double multiplier =
                            getMutualPotentialEffectiveCoefficientMultiplier( degreeOfBody1,
                                                                              signedOrderOfBody1,
                                                                              degreeOfBody2,
                                                                              signedOrderOfBody2,
                                                                              accelerationBetweenBodies_->getAreCoefficientsNormalized( ) );

                    Eigen::Vector3d angularMomentumTransformedCosineCoefficientsBody2;
                    Eigen::Vector3d angularMomentumTransformedSineCoefficientsBody2;
                    for( int k = 0; k < 3; k++ )
                    {
                        angularMomentumTransformedCosineCoefficientsBody2( k ) =
                                transformedCosineCoefficientsBody2AngularMomentum_.at( k )( degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                        angularMomentumTransformedSineCoefficientsBody2( k ) =
                                transformedSineCoefficientsBody2AngularMomentum_.at( k )( degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                    }

                    const Eigen::Vector3d effectiveAngularMomentumCosineCoefficients =
                            ( body1CosineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 -
                              signOrderBody1 * signOrderBody2 * body1SineCoefficient * angularMomentumTransformedSineCoefficientsBody2 ) *
                            multiplier;
                    const Eigen::Vector3d effectiveAngularMomentumSineCoefficients =
                            ( signOrderBody2 * body1CosineCoefficient * angularMomentumTransformedSineCoefficientsBody2 +
                              signOrderBody1 * body1SineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 ) *
                            signTotalOrder * multiplier;

                    const Eigen::Vector3d currentContribution = equatorialRadiusRatioPower * legendrePolynomial *
                            ( effectiveAngularMomentumCosineCoefficients * cosineOfMultipleLongitude +
                              effectiveAngularMomentumSineCoefficients * sineOfMultipleLongitude );
                    body2TorqueInBodyFixedFrameOfBody1 += currentContribution;
                }
            }
            body2TorqueInBodyFixedFrameOfBody1 *= -preMultiplier;

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
            currentTorque_ *= bodyUndergoingTorqueMassFunction_( );

            currentTime_ = currentTime;
        }
    }

    //! Return the current torque (body-fixed frame of body undergoing torque).
    Eigen::Vector3d getTorque( ) override
    {
        return currentTorque_;
    }

    //! Reset internal time/cache state.
    virtual void resetCurrentTime( ) override
    {
        accelerationBetweenBodies_->resetCurrentTime( );
        currentTime_ = TUDAT_NAN;
    }
    //! Return selected coefficient combinations used by this torque model.
    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& getCoefficientCombinationsToUse( ) const
    {
        return coefficientCombinationsToUse_;
    }

    //! Return associated full two-body acceleration model.
    std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > getAccelerationBetweenBodies( ) const
    {
        return accelerationBetweenBodies_;
    }

    //! Return whether this instance outputs torque for body 1.
    bool getAcceleratedBodyIsBody1( ) const
    {
        return acceleratedBodyIsBody1_;
    }

    double getBodyUndergoingTorqueMass( ) const
    {
        return bodyUndergoingTorqueMassFunction_( );
    }

    //! Compute \hat{J}\bar{C}_{l,m}^{2,F1} and \hat{J}\bar{S}_{l,m}^{2,F1} from body-2 coefficients and supplied Wigner D-matrices.
    /*!
     * Helper used by this model and analytical partials to evaluate the \hat{J}-mapped coefficient fields
     * entering Dirkx et al. (2019), Eq. (67) (through Eq. (60)).
     */
    void computeTransformedAngularMomentumCoefficients( const Eigen::MatrixXd& cosineCoefficientsBody2,
                                                        const Eigen::MatrixXd& sineCoefficientsBody2,
                                                        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
                                                        const bool areCoefficientsNormalized,
                                                        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
                                                        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum );

private:
    Eigen::Vector3d currentTorque_;

    std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > accelerationBetweenBodies_;

    std::function< double( ) > bodyUndergoingTorqueMassFunction_;

    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > body2TorqueCombinationsToUse_;

    std::array< Eigen::MatrixXd, 3 > transformedCosineCoefficientsBody2AngularMomentum_;
    std::array< Eigen::MatrixXd, 3 > transformedSineCoefficientsBody2AngularMomentum_;

    bool acceleratedBodyIsBody1_;
};

}  // namespace gravitation

}  // namespace tudat

#endif  // MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
