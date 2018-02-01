/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/math/special_functions/erf.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/Statistics/continuousProbabilityDistributions.h"

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
