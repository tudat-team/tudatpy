#ifndef MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
#define MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>


#include <Eigen/Core>
#include <Eigen/Geometry>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebra.h>
#include <Tudat/Mathematics/BasicMathematics/legendrePolynomials.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Tudat/Astrodynamics/Gravitation/mutualExtendedBodySphericalHarmonicAcceleration.h"

namespace tudat
{

namespace gravitation
{

class MutualExtendedBodySphericalHarmonicTorque: public basic_astrodynamics::TorqueModel
{

public:

    MutualExtendedBodySphericalHarmonicTorque(
            const boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > accelerationBetweenBodies,
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

    boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > accelerationBetweenBodies_;

    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    bool acceleratedBodyIsBody1_;

};

}

}

#endif // MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
