#ifndef MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
#define MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP

#include <boost/tuple/tuple.hpp>
#include <memory>


#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/legendrePolynomials.h"

#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/gravitation/mutualExtendedBodySphericalHarmonicAcceleration.h"

namespace tudat
{

namespace gravitation
{

class MutualExtendedBodySphericalHarmonicTorque: public basic_astrodynamics::TorqueModel
{

public:

    MutualExtendedBodySphericalHarmonicTorque(
            const std::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > accelerationBetweenBodies,
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



    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getCoefficientCombinationsToUse( )
    {
        return coefficientCombinationsToUse_;
    }

private:

    void calculateDirectTorque( )
    {

    }

    Eigen::Vector3d currentTorque_;

    std::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > accelerationBetweenBodies_;

    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    bool acceleratedBodyIsBody1_;

};

}

}

#endif // MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
