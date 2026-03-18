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
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
                coefficientCombinationsToUse );

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
class FullTwoBodySphericalHarmonicTorque: public basic_astrodynamics::TorqueModel
{

public:

    //! Constructor.
    FullTwoBodySphericalHarmonicTorque(
            const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > accelerationBetweenBodies,
            const bool acceleratedBodyIsBody1 ):
        accelerationBetweenBodies_( accelerationBetweenBodies ),
        acceleratedBodyIsBody1_( acceleratedBodyIsBody1 )
    {
        coefficientCombinationsToUse_ = accelerationBetweenBodies_->getEffectiveMutualPotentialField( )->getCoefficientCombinationsToUse( );
        body2TorqueCombinationsToUse_ = gravitation::getBody2TorqueCombinationsToUse( coefficientCombinationsToUse_ );
    }

    //! Update torque to current epoch.
    void updateMembers( const double currentTime = TUDAT_NAN );

    //! Return the current torque (body-fixed frame of body undergoing torque).
    Eigen::Vector3d getTorque( )
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
    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
    getCoefficientCombinationsToUse( ) const
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

    //! Compute \hat{J}\bar{C}_{l,m}^{2,F1} and \hat{J}\bar{S}_{l,m}^{2,F1} from body-2 coefficients and supplied Wigner D-matrices.
    /*!
     * Helper used by this model and analytical partials to evaluate the \hat{J}-mapped coefficient fields
     * entering Dirkx et al. (2019), Eq. (67) (through Eq. (60)).
     */
    void computeTransformedAngularMomentumCoefficients(
            const Eigen::MatrixXd& cosineCoefficientsBody2,
            const Eigen::MatrixXd& sineCoefficientsBody2,
            const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
            const bool areCoefficientsNormalized,
            std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
            std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum );

private:
    Eigen::Vector3d currentTorque_;

    std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > accelerationBetweenBodies_;

    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > body2TorqueCombinationsToUse_;

    std::array< Eigen::MatrixXd, 3 > transformedCosineCoefficientsBody2AngularMomentum_;
    std::array< Eigen::MatrixXd, 3 > transformedSineCoefficientsBody2AngularMomentum_;

    bool acceleratedBodyIsBody1_;

};

}

}

#endif // MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
