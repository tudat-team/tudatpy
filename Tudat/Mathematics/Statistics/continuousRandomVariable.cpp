#include <boost/math/special_functions/erf.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/Statistics/continuousRandomVariable.h"

namespace tudat
{

namespace statistics
{

//! Function to evaluate pdf of Gaussian distribution.
double evaluateGaussianPdf( const double independentVariable, const double mean, const double standardDeviation )
{
    double offsetFromMean = independentVariable - mean;
    return 1.0 / ( std::sqrt( 2.0 * mathematical_constants::PI ) * standardDeviation ) *
            std::exp( -( offsetFromMean * offsetFromMean ) / ( 2.0 * standardDeviation * standardDeviation ) );
}

//! Function to evaluate cdf of Gaussian distribution.
double calculateGaussianCdf( const double independentVariable, const double mean, const double standardDeviation )
{
    return 0.5 * ( 1.0 + boost::math::erf( ( independentVariable - mean ) / ( std::sqrt( 2.0 ) * ( standardDeviation ) ) ) );
}

}

}
