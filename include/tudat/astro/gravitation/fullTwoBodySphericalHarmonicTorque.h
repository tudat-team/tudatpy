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

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getBody2TorqueCombinationsToUse(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
                coefficientCombinationsToUse );

void computeTransformedAngularMomentumCoefficientsFromWignerMatrices(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum );

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

    FullTwoBodySphericalHarmonicTorque(
            const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > accelerationBetweenBodies,
            const bool acceleratedBodyIsBody1 ):
        accelerationBetweenBodies_( accelerationBetweenBodies ),
        acceleratedBodyIsBody1_( acceleratedBodyIsBody1 )
    {
        coefficientCombinationsToUse_ = accelerationBetweenBodies_->getEffectiveMutualPotentialField( )->getCoefficientCombinationsToUse( );
        body2TorqueCombinationsToUse_ = gravitation::getBody2TorqueCombinationsToUse( coefficientCombinationsToUse_ );
    }

    void updateMembers( const double currentTime = TUDAT_NAN );

    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    virtual void resetCurrentTime( ) override
    {
        accelerationBetweenBodies_->resetCurrentTime( );
        currentTime_ = TUDAT_NAN;
    }



    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
    getCoefficientCombinationsToUse( ) const
    {
        return coefficientCombinationsToUse_;
    }

    std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > getAccelerationBetweenBodies( ) const
    {
        return accelerationBetweenBodies_;
    }

    bool getAcceleratedBodyIsBody1( ) const
    {
        return acceleratedBodyIsBody1_;
    }

    //! Compute \hat{J}\bar{C}_{l,m}^{2,F1} and \hat{J}\bar{S}_{l,m}^{2,F1} from body-2 coefficients and supplied Wigner D-matrices.
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
