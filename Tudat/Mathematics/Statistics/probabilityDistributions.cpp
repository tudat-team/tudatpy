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

#include <Tudat/Mathematics/Statistics/probabilityDistributions.h>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>

namespace tudat
{
namespace statistics
{

using tudat::mathematical_constants::PI;

//! Constructor.
GaussianDistributiond::GaussianDistributiond(double Mean, double StandardDeviation){ // constructor
    mean_ = Mean ;
    standardDeviation_ = StandardDeviation ;
    variance_ = std::pow( standardDeviation_ , 2.0 );
}

//! Get probability density of 1D Gaussian distribution
double GaussianDistributiond::getProbabilityDensity(const double& x){
    return (std::exp( ( - std::pow( x - mean_ , 2.0 ) ) / ( 2.0 * variance_ ) )
            /( std::sqrt( 2.0 * PI ) * standardDeviation_ ) ) ;
}

//! Get cumulative probability of 1D Gaussian distribution
double GaussianDistributiond::getCumulativeProbability(const double &x){
    boost::math::normal distribution(mean_, standardDeviation_);
    return boost::math::cdf( distribution , x );
}

//! Get Quantile (Inverse CDF).
double GaussianDistributiond::getQuantile(const double &x )
{
    boost::math::normal distribution(mean_, standardDeviation_);
    return boost::math::quantile( distribution , x );
}

//! Constructor.
UniformDistributiond::UniformDistributiond(double LowerBound, double UpperBound){ // constructor
    lowerBound_ = LowerBound ;
    upperBound_ = UpperBound ;
    probabilityDensity_ = 1.0/( upperBound_ - lowerBound_ ) ;
}


//! Get probability density of 1D Uniform distribution
double UniformDistributiond::getProbabilityDensity(const double& x){
    if ( (x >= lowerBound_ && x <= upperBound_ ) ){
        return probabilityDensity_ ;
    }
    else{
        return 0.0 ;
    }
}

//! Compute mean and standard deviation 1D Uniform distribution.
void UniformDistributiond::computeMeanAndStandardDeviation(){
    mean_ = (upperBound_ + lowerBound_)/2.0 ;
    variance_ = std::pow(upperBound_ - lowerBound_,2.0) / 12.0;
    standardDeviation_ = std::sqrt( variance_ );
}

//! Constructor of the lognormal distribution
LogNormalDistributiond::LogNormalDistributiond( double Mean, double StandardDeviation )
{
    mean_ = Mean;
    standardDeviation_ = StandardDeviation;
    variance_ = std::pow( standardDeviation_ , 2.0 );

    locationParameter_ = std::log( mean_ / ( std::sqrt( variance_ / std::pow( mean_ , 2.0 ) + 1.0 ) ) );
    scaleParameter_ = std::sqrt( std::log( variance_ / std::pow( mean_ , 2.0 ) + 1.0 ) ) ;
}

//! Get probability density.
double LogNormalDistributiond::getProbabilityDensity( const double& x )
{
    boost::math::lognormal distribution( locationParameter_ , scaleParameter_ );
    return boost::math::pdf( distribution , x );
}

//! Get cumulative probability.
double LogNormalDistributiond::getCumulativeProbability(const double &x )
{
    boost::math::lognormal distribution( locationParameter_ , scaleParameter_ );
    return boost::math::cdf( distribution , x );
}

//! Get Quantile (Inverse CDF).
double LogNormalDistributiond::getQuantile(const double &x )
{
    boost::math::lognormal distribution( locationParameter_ , scaleParameter_ );
    return boost::math::quantile( distribution , x );
}

//! Get probability density of gaussian copula.
double GaussianCopulaDistributionXd::getProbabilityDensity( const Eigen::VectorXd& x ){
    double probabilityDensity = 0.0 ;

    // Check if vector x is inside [0,1]
    int InBound = 0 ;
    for( int i = 0 ; i < dimension_ ; i++ ){ // check in bounds
        if( x(i) > 0.0 && x(i) < 1.0 ){
            InBound++ ;
        }
    }

    if( InBound == dimension_ ){
        // Convert U[0,1] to N[0,1] using inverse CDF of standard normal distribution
        Eigen::VectorXd y( dimension_ ) ;
        boost::math::normal distribution( 0.0 , 1.0 );

        for(int i = 0 ; i < dimension_ ; i++){
            y(i) = boost::math::quantile( distribution , x(i) ); // Inverse cdf
        }

        // Calculate probability density
        Eigen::MatrixXd location = - 0.5 * ( y.transpose() *
                       (inverseCorrelationMatrix_ - Eigen::MatrixXd::Identity( dimension_ , dimension_ ) ) * y ) ;

        probabilityDensity = ( ( 1.0/( std::sqrt( determinant_ ) ) )*std::exp( location(0,0) ) ) ;
    }

    return probabilityDensity ;
}

} // namespace statistics
} // namespace tudat
