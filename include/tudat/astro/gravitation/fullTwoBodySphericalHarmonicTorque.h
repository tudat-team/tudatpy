#ifndef MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
#define MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP

#include <tuple>
#include <memory>


#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/legendrePolynomials.h"

#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicAcceleration.h"

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



    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getCoefficientCombinationsToUse( )
    {
        return coefficientCombinationsToUse_;
    }

private:

    void calculateDirectTorque( )
    {

    }

    Eigen::Vector3d currentTorque_;

    std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > accelerationBetweenBodies_;

    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    bool acceleratedBodyIsBody1_;

};

}

}

#endif // MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
