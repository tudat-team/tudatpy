#ifndef MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
#define MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP

#include <array>
#include <complex>
#include <tuple>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicAcceleration.h"
#include "tudat/math/basic/wignerDMatrices.h"

namespace tudat
{

namespace gravitation
{

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



    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getCoefficientCombinationsToUse( )
    {
        return coefficientCombinationsToUse_;
    }

private:
    //! Compute \hat{J} D_{m,k}^l in Cartesian basis of frame F1.
    Eigen::Vector3cd computeAngularMomentumOperatorOnWignerCoefficient(
            const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
            const int degree,
            const int orderM,
            const int orderK );

    //! Compute \hat{J}\bar{C}_{l,m}^{2,F1} and \hat{J}\bar{S}_{l,m}^{2,F1} from body-2 coefficients and current Wigner D-matrices.
    void computeTransformedAngularMomentumCoefficients(
            const Eigen::MatrixXd& cosineCoefficientsBody2,
            const Eigen::MatrixXd& sineCoefficientsBody2,
            const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
            const bool areCoefficientsNormalized,
            std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
            std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum );

    Eigen::Vector3d currentTorque_;

    std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > accelerationBetweenBodies_;

    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    bool acceleratedBodyIsBody1_;

};

}

}

#endif // MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
