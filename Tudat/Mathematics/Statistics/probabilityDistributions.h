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

namespace tudat
{
namespace statistics
{

//! Base class for probability distributions.
/*!
 * Base class for general probability distributions.
 * \tparam IndependentVariableType The type of the independent (random) variable.
 */
template< typename IndependentVariableType >
class ProbabilityDistribution
{
public:

    //! Destructor
    virtual ~ProbabilityDistribution(){}

    //! Get probability density.
    /*!
     * This function computes and returns the probability density at X = x.
     * The function uses the probability density function to compute this value.
     * This function is a pure virtual function, so it must be implemented by the derived class.
     * \param IndependentVariableType x sample of random variable.
     * \return probability density.
     */
    virtual double getProbabilityDensity( const IndependentVariableType& x ) = 0 ;

    //! Get Cumulative Probability
    /*!
     * The function uses the cumulative distribution function to compute the cumulative probability.
     * F(x) = P(X < x)
     * \param IndependentVariableType x sample of random variable.
     * \return cumulative probability.
     */
    virtual double getCumulativeProbability( const IndependentVariableType& x )
    {
        return TUDAT_NAN;
    }

    //! Get Quantile (REMOVE ONLY 1D)
    /*!
     * This function uses the cumulative probability.
     * Inverse CDF
     * \param IndependentVariableType x sample of random variable.
     * \return quantile.
     */
    virtual double getQuantile( const IndependentVariableType& x )
    {
        return TUDAT_NAN;
    }

protected:

private:
};

typedef boost::shared_ptr< ProbabilityDistribution<double> > ProbabilityDistributionDoublePointer;

typedef boost::shared_ptr< ProbabilityDistribution<Eigen::VectorXd> > ProbabilityDistributionXdPointer;

//! One-dimensional Gaussian distribution class.
/*!
 * One-dimensional Gaussian distribution class.
 * Source: Dekking, F.; Kraaikamp, C.; Lopuhaa, H. & Meester, L. A Modern Introduction to Probability and Statistics Springer, 2005
 */
class GaussianDistributiond: public ProbabilityDistribution<double>
{
public:

    //! Constructor.
    GaussianDistributiond( double Mean, double StandardDeviation );

    //! Get probability density.
    double getProbabilityDensity( const double& x );

    //! Get cumulative probability.
    double getCumulativeProbability( const double& x );

    //! Get Quantile (Inverse CDF)
    double getQuantile( const double& x );

    //! Set the properties (mean and standard deviation) of the Gaussian distribution
    void setProperties( const double& mean, double standardDeviation )
    {
        mean_ = mean;
        standardDeviation_ = standardDeviation;
    }

protected:
private:
    //! Mean value of the distribution.
    double mean_ ;

    //! Standard deviation value of the distribution.
    double standardDeviation_ ;

    //! Variance value of the distribution.
    double variance_ ;
};

//! One-dimensional Uniform distribution class.
/*!
 * One-dimensional Uniform distribution class.
 * Source: Dekking, F.; Kraaikamp, C.; Lopuhaa, H. & Meester, L. A Modern Introduction to Probability and Statistics Springer, 2005
 */
class UniformDistributiond: public ProbabilityDistribution<double>
{
public:

    //! Constructor.
    UniformDistributiond( double LowerBound, double UpperBound );

    //! Compute mean and standard deviation of distribution.
    void computeMeanAndStandardDeviation();

    //! Get mean value of distribution
    double getMean(){
        computeMeanAndStandardDeviation();
        return mean_;
    }

    //! Get standard deviation of distribution
    double getStandardDeviation(){
        computeMeanAndStandardDeviation();
        return standardDeviation_;
    }

    //! Get probability density.
    double getProbabilityDensity(const double& x);

protected:
private:
    //! Lower bound value of the distribution.
    double lowerBound_ ;

    //! Upper bound value of the distribution.
    double upperBound_ ;

    //! Mean value of the distribution.
    double mean_ ;

    //! Standard deviation value of the distribution.
    double standardDeviation_ ;

    //! Variance value of the distribution.
    double variance_ ;

    //! Probability density value of distribution (constant).
    double probabilityDensity_ ;
};

//! One-dimensional Lognormal distribution class.
/*!
 * One-dimensional Lognormal distribution class.
 * Source: Montgomery, D. C. & Runger, G. C. Applied Statistics and Probability for engineers Wiley, 2014
 */
class LogNormalDistributiond: public ProbabilityDistribution<double>
{
public:

    //! Constructor.
    LogNormalDistributiond( double Mean, double StandardDeviation );

    //! Get probability density.
    double getProbabilityDensity( const double& x );

    //! Get cumulative probability.
    double getCumulativeProbability( const double& x );

    //! Get Quantile (Inverse CDF)
    double getQuantile( const double& x );

    //! Get location parameter.
    double getLocationParameter()
    {
        return locationParameter_;
    }

    //! Get scale parameter.
    double getScaleParameter()
    {
        return scaleParameter_;
    }


protected:
private:
    //! Mean value of the distribution.
    double mean_ ;

    //! Standard deviation value of the distribution.
    double standardDeviation_ ;

    //! Variance value of the distribution.
    double variance_ ;

    //! Location parameter: Mean value of the log of the random variable.
    double locationParameter_ ;

    //! Scale parameter: Standard deviation of the log of the random variable.
    double scaleParameter_ ;
};


//! Multi-dimensional Gaussian Distribution class.
/*!
 * Multi-dimensional Gaussian Distribution class.
 * Source: Tong, Y. The Multivariate Normal Distribution Springer-Verslag, 1990
 */
class GaussianDistributionXd: public ProbabilityDistribution<Eigen::VectorXd>
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
    double getProbabilityDensity( const Eigen::VectorXd& x ){
        using tudat::mathematical_constants::PI;
        Eigen::VectorXd u = (x - mean_) ;
        Eigen::MatrixXd location = -0.5*( u.transpose()*inverseCovarianceMatrix_*u ) ;

        double probability = std::exp( location(0,0) ) /( std::pow(2.0*PI,dimension_/2.0) * std::sqrt(determinant_) ) ;
        return probability ;
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
class GaussianCopulaDistributionXd: public ProbabilityDistribution< Eigen::VectorXd >
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
    double getProbabilityDensity( const Eigen::VectorXd& x );

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








// //! Distribution of a random variable with uniform marginals and correlation using a Gaussian copula.
//template< int Dimension >
//class UniformCorrelatedDistributionXd: public ProbabilityDistribution<Eigen::VectorXd> {
//public:

//    UniformCorrelatedDistributionXd(Eigen::VectorXd Mean , Eigen::MatrixXd CovarianceMatrix ){
//        mean_ = Mean ;
//        covarianceMatrix_ = CovarianceMatrix ;

//        // correlationmatrix
//        correlationMatrix_ = Thesis::Statistics::Basics::Covariance2CorrelationMatrix(covarianceMatrix_) ;

//        dimension_ = mean_.rows() ;

//        // determinant correlation matrix
//        determinantCor = correlationMatrix_.determinant() ;

//        // inverse function needs to know what the dimension of the matrix is at compile time
//        typedef Eigen::Matrix<double,Dimension,Dimension> MatrixNd;
//        MatrixNd Correlation = correlationMatrix_ ;
//        MatrixNd Inversecorrelation = Correlation.inverse() ;
//        inverseCorrelationMatrix_ = Inversecorrelation ;

//        // Calculate bounds
//        lowerBound_.resize( Dimension ) ;
//        upperBound_.resize( Dimension ) ;
//        for(int i = 0 ; i < Dimension ; i++){
//            lowerBound_(i) = mean_(i) - std::sqrt(3.0 * covarianceMatrix_(i,i) ) ;
//            upperBound_(i) = mean_(i) + std::sqrt(3.0 * covarianceMatrix_(i,i) ) ;
//        }
//        width = upperBound_ - lowerBound_ ;
//        volume = width.prod() ; // product of marginal PDFs
//    }

//    double getProbabilityDensity(Eigen::VectorXd x){
//        double probability = 0.0 ;

//        int InBound = 0 ;
//        for(int i = 0 ; i < dimension_ ; i++){ // check in bounds
//            if( x(i) > lowerBound_(i) && x(i) < upperBound_(i) ){
//                InBound++ ;
//            }
//        }
//        if(InBound == dimension_){
//            // Convert uniform to U = [0,1] :  F(x) = x - a / (b-a) -> U = F_X(X)
//            Eigen::VectorXd Fx = (x - lowerBound_).cwiseQuotient( width ) ;

//            // Convert U[0,1] to N[0,1] using inverse CDF of standard normal distribution
//            Eigen::VectorXd y(Fx.rows()) ;
//            for(int i = 0 ; i < Fx.rows() ; i++){
//                y(i) = gsl_cdf_ugaussian_Pinv(Fx(i)) ; // Inverse standard normal CDF
//            }

//            // Calculate probability density
//            Eigen::MatrixXd location = -0.5*( y.transpose()*(inverseCorrelationMatrix_ - Eigen::MatrixXd::Identity(Dimension,Dimension) )*y ) ;

//            probability = ( ( 1.0/( volume * std::sqrt(determinantCor) ) )*std::exp( location(0,0) ) ) ;
//        }

//        return probability ;
//    }

//private:

//    //! Dimension of the distribution
//    int                 dimension_           ;

//    //! Vector of mean values of X
//    Eigen::VectorXd     mean_                  ;

//    //! Covariance matrix of X
//    Eigen::MatrixXd     covarianceMatrix_    ;

//    //! Correlation matrix of X
//    Eigen::MatrixXd     correlationMatrix_   ;

//    //! Determinant of correlation matrix
//    double              determinantCor      ;

//    //! Inverse of the correlation matrix
//    Eigen::MatrixXd     inverseCorrelationMatrix_  ;

//    //! Lowerbound of the uniform marginals
//    Eigen::VectorXd     lowerBound_           ;

//    //! Upperbound of the uniform marginals
//    Eigen::VectorXd     upperBound_          ;

//    //! Width of the uniform marginals
//    Eigen::VectorXd     width               ;

//    //! Volume of the distribution (product of the marginal PDFs)
//    double              volume              ;
//};


} // namespace statistics
} // namespace tudat

#endif // TUDAT_PROBABILITY_DISTRIBUTIONS_H
