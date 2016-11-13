#include "Tudat/Astrodynamics/Gravitation/mutualSphericalHarmonicGravityModel.h"

namespace tudat
{

namespace gravitation
{

//! Function to manually remove the C(0,0) term from cosine coefficients,
Eigen::MatrixXd setDegreeAndOrderCoefficientToZero(
        const boost::function< Eigen::MatrixXd( ) > originalCosineCoefficientFunction )
{
    Eigen::MatrixXd newCoefficients = originalCosineCoefficientFunction( );
    newCoefficients( 0, 0 ) = 0.0;
    return newCoefficients;
}

}

}
