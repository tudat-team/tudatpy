/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      090807    J. Melman         File created.
 *      100930    D. Dirkx          Modified to comply with Tudat standards
 *      100930    J. Melman         Implemented namespace, minor comment
 *                                  changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120202    K. Kumar          Moved functions from linearAlgebra.cpp.
 *
 *    References
 *
 *    Notes
 *
 */

#include <cmath>
#include <numeric>

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>

#include <Tudat/Mathematics/Statistics/multiVariateGaussianProbabilityDistributions.h>


namespace tudat
{
namespace statistics
{


//! Function to evaluate pdf of Gaussian cupola distribution
double GaussianCopulaDistributionXd::evaluatePdf(
        const Eigen::VectorXd& independentVariables )
{
    double probabilityDensity = 0.0 ;

    // Check if vector independentVariables is inside [0,1]
    int inBound = 0 ;
    for( int i = 0 ; i < dimension_ ; i++ )
    {
        if( independentVariables(i) > 0.0 && independentVariables(i) < 1.0 )
        {
            inBound++;
        }
    }

    // If data is in bounds
    if( inBound == dimension_ )
    {
        // Convert U[0,1] to N[0,1] using inverse CDF of standard normal distribution
        Eigen::VectorXd gaussianQuantiles( dimension_ ) ;
        boost::math::normal distribution( 0.0 , 1.0 );

        for( int i = 0 ; i < dimension_ ; i++ )
        {
            gaussianQuantiles( i ) = boost::math::quantile( distribution , independentVariables( i ) ); // Inverse cdf
        }

        // Calculate probability density
        Eigen::MatrixXd location = - 0.5 * ( gaussianQuantiles.transpose() *
                       ( inverseCorrelationMatrix_ - Eigen::MatrixXd::Identity( dimension_ , dimension_ ) ) *
                                             gaussianQuantiles );

        probabilityDensity = ( ( 1.0 / ( std::sqrt( determinant_ ) ) ) * std::exp( location( 0, 0 ) ) ) ;
    }

    return probabilityDensity ;
}

} // namespace statistics
} // namespace tudat
