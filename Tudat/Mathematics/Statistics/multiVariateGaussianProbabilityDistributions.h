/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_PROBABILITY_DISTRIBUTIONS_H
#define TUDAT_PROBABILITY_DISTRIBUTIONS_H

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/Statistics/continuousProbabilityDistributions.h"
namespace tudat
{
namespace statistics
{

//! Implementation of multi-dimensional Gaussian Probability Distribution.
/*!
 *  Implementation of multi-dimensional Gaussian Probability Distribution. Only pdf is presently implemented,
 *  a runtime error is thrown if the cdf function is called.
 *  Model taken from Tong, Y. The Multivariate Normal Distribution Springer-Verslag, 1990
 */
class GaussianDistributionXd: public ContinuousProbabilityDistribution< Eigen::VectorXd >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param mean Mean values of distribution
     * \param covarianceMatrix Covariance matrix of distribution
     */
    GaussianDistributionXd(
            const Eigen::VectorXd& mean , const Eigen::MatrixXd& covarianceMatrix ):
        mean_( mean ), covarianceMatrix_( covarianceMatrix )
    {
        if( covarianceMatrix.rows( ) != covarianceMatrix.cols( ) )
        {
            throw std::runtime_error( "Error, covarianceMatrix input to GaussianDistributionXd is not square" );
        }

        dimension_ = static_cast< double >( mean_.rows( ) );
        determinant_ = covarianceMatrix_.determinant( );
        inverseCovarianceMatrix_ = covarianceMatrix_.inverse( );
    }

    //! Function to evaluate pdf of multivariate Gaussian distribution
    /*!
     *  Function to evaluate probability distribution function of multivariate Gaussian distribution at given
     *  list of independent variable values.
     *  \param independentVariables List of independent variable values
     *  \return Evaluated multivariate Gaussian pdf
     */
    double evaluatePdf( const Eigen::VectorXd& independentVariables )
    {
        Eigen::VectorXd distanceFromMean = ( independentVariables - mean_ ) ;
        Eigen::MatrixXd location = -0.5 * (
                    distanceFromMean.transpose( ) * inverseCovarianceMatrix_ * distanceFromMean ) ;

        return std::exp( location( 0, 0 ) ) / ( std::pow( 2.0 * mathematical_constants::PI, dimension_ / 2.0 ) *
                                                std::sqrt( determinant_ ) );
    }

    //! Function to evaluate cdf of multivariate Gaussian distribution (not yet implemented)
    /*!
     *  Function to evaluate cumulative distribution function of multivariate Gaussian distribution at given
     *  list of independent variable values.
     *  NOTE: function not yet implemented
     *  \param independentVariables List of independent variable values
     *  \return Evaluated multivariate Gaussian cdf.
     */
    double evaluateCdf( const Eigen::VectorXd& independentVariables )
    {
        throw std::runtime_error( "Cdf of GaussianDistributionXd not yet implemented" );

        return TUDAT_NAN;
    }

private:
    //! Dimension of the random variable X
    double dimension_;

    //! Mean vector of random variable X
    Eigen::VectorXd mean_       ;

    //! Covariance matrix of random variable X
    Eigen::MatrixXd covarianceMatrix_    ;

    //! Determinant of covariance matrix
    double determinant_;

    //! Inverse of the covariance matrix
    Eigen::MatrixXd inverseCovarianceMatrix_ ;
};


//! Implementation of Gaussian Copula Distribution.
/*!
 *  Implementation of Gaussian Copula Distribution. A Gaussian copula can be used to link several marginal distributions
 *  to form a joint distribution. Only pdf is presently implemented,
 *  a runtime error is thrown if the cdf function is called.
 *  Source: Song, P. X.-K. Multivariate Dispersion Models Generated from Gaussian Copula Scandinavian Journal of
 *  Statistics, 2000, 27, 305-320
 */
class GaussianCopulaDistributionXd: public ContinuousProbabilityDistribution< Eigen::VectorXd >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param correlationMatrix Correlation matrix of distribution
     */
    GaussianCopulaDistributionXd( const Eigen::MatrixXd& correlationMatrix ):
        correlationMatrix_( correlationMatrix )
    {
        if( correlationMatrix.rows( ) != correlationMatrix.cols( ) )
        {
            throw std::runtime_error( "Error, correlationMatrix input to GaussianCopulaDistributionXd is not square" );
        }

        dimension_ = correlationMatrix.rows( );

        inverseCorrelationMatrix_ = correlationMatrix_.inverse( ) ;
        determinant_ = correlationMatrix_.determinant( ) ;
    }

    //! Function to evaluate pdf of Gaussian cupola distribution
    /*!
     *  Function to evaluate probability distribution function of Gaussian cupola distribution at given
     *  list of independent variable values.
     *  \param independentVariables List of independent variable values
     *  \return Evaluated multivariate Gaussian pdf
     */
    double evaluatePdf( const Eigen::VectorXd& independentVariables );

    //! Function to evaluate cdf of Gaussian cupola distribution (not yet implemented)
    /*!
     *  Function to evaluate cumulative distribution function of Gaussian cupola distribution at given
     *  list of independent variable values.
     *  NOTE: function not yet implemented
     *  \param independentVariables List of independent variable values
     *  \return Evaluated Gaussian cupola cdf.
     */
    double evaluateCdf( const Eigen::VectorXd& independentVariables )
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
