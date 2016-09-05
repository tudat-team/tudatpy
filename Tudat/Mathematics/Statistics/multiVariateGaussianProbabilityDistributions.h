/*    Copyright (c) 2010-2016, Delft University of Technology
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
 *      160429    R. Hoogendoorn    File created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_PROBABILITY_DISTRIBUTIONS_H
#define TUDAT_PROBABILITY_DISTRIBUTIONS_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Householder>
#include <Eigen/QR>
#include <Eigen/Sparse>

#include <boost/shared_ptr.hpp>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Mathematics/Statistics/continuousProbabilityDistributions.h>

namespace tudat
{
namespace statistics
{

//! Multi-dimensional Gaussian Distribution class.
/*!
 * Multi-dimensional Gaussian Distribution class.
 * Source: Tong, Y. The Multivariate Normal Distribution Springer-Verslag, 1990
 */
class GaussianDistributionXd: public ContinuousProbabilityDistribution< Eigen::VectorXd >
{
public:

    //! Constructor
    GaussianDistributionXd(Eigen::VectorXd Mean , Eigen::MatrixXd CovarianceMatrix ){ // constructor
        mean_ = Mean ;
        covarianceMatrix_ = CovarianceMatrix ;
        dimension_ = double(mean_.rows()) ;
        determinant_ = CovarianceMatrix.determinant() ;
        inverseCovarianceMatrix_ = covarianceMatrix_.inverse() ;
    }

    //! Get probability density
    double evaluatePdf( const Eigen::VectorXd& x )
    {
        using tudat::mathematical_constants::PI;
        Eigen::VectorXd u = (x - mean_) ;
        Eigen::MatrixXd location = -0.5*( u.transpose()*inverseCovarianceMatrix_*u ) ;

        double probability = std::exp( location(0,0) ) /( std::pow(2.0*PI,dimension_/2.0) * std::sqrt(determinant_) ) ;
        return probability ;
    }

    double evaluateCdf( const Eigen::VectorXd& independentVariable )
    {
        throw std::runtime_error( "Cdf of GaussianDistributionXd not yet implemented" );

        return TUDAT_NAN;
    }

private:
    //! Dimension of the random variable X
    double              dimension_           ;

    //! Mean vector of random variable X
    Eigen::VectorXd     mean_                  ;

    //! Covariance matrix of random variable X
    Eigen::MatrixXd     covarianceMatrix_    ;

    //! Determinant of covariance matrix
    double              determinant_         ;

    //! Inverse of the covariance matrix
    Eigen::MatrixXd     inverseCovarianceMatrix_ ;
};


//! Gaussian Copula Distribution class.
/*!
 * A Gaussian copula can be used to link several marginal distributions to form a joint distribution.
 * Source: Song, P. X.-K. Multivariate Dispersion Models Generated from Gaussian Copula Scandinavian Journal of Statistics, 2000, 27, 305-320
 */
class GaussianCopulaDistributionXd: public ContinuousProbabilityDistribution< Eigen::VectorXd >
{
public:

    //! Constructor
    GaussianCopulaDistributionXd( int dimension , Eigen::MatrixXd correlationMatrix )
    {
        dimension_ = dimension;
        correlationMatrix_ = correlationMatrix;

        inverseCorrelationMatrix_ = correlationMatrix_.inverse() ;
        determinant_ = correlationMatrix_.determinant() ;
    }

    //! Get probability density.
    double evaluatePdf( const Eigen::VectorXd& x );

    double evaluateCdf( const Eigen::VectorXd& independentVariable )
    {
        throw std::runtime_error( "Cdf of GaussianCopulaDistributionXd not yet implemented" );
        return TUDAT_NAN;
    }

protected:

private:

    //! Dimension of the copula
    int dimension_;

    //! Correlation matrix
    Eigen::MatrixXd correlationMatrix_;

    //! Inverse of the correlation matrix
    Eigen::MatrixXd inverseCorrelationMatrix_;

    //! determinant of the correlation matrix
    double determinant_;
};

} // namespace statistics
} // namespace tudat

#endif // TUDAT_PROBABILITY_DISTRIBUTIONS_H
