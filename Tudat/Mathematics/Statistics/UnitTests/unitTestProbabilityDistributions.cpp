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
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical( ) function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean( )
 *                                  and computeSampleVariance( ) functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120202    K. Kumar          Moved unit tests from unitTestBasicMathematicsFunctions.h/.cpp;
 *                                  rewrote unit tests using Boost unit test framework.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include <boost/math/distributions/lognormal.hpp>

#include "Tudat/Mathematics/Statistics/probabilityDistributions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_probability_distributions )

//! Test if 1D uniform distribution class works correctly.
BOOST_AUTO_TEST_CASE( testUniformDistributionDouble )
{
    using namespace tudat::statistics;
    UniformDistributiond uniformDistribution(0.0,2.0) ;

    ProbabilityDistributionDoublePointer distributionPointer = boost::make_shared< UniformDistributiond >(0.0,3.0) ;

    double probabilityDensity = uniformDistribution.getProbabilityDensity(0.8) ;

    // Check if probability density is correctly calculated
    BOOST_CHECK_SMALL( std::fabs( probabilityDensity - 0.5 ) ,
                       std::numeric_limits<double>::epsilon() ) ;

    BOOST_CHECK_SMALL( std::fabs( distributionPointer->getProbabilityDensity(1.0) - (1.0/3.0) ) ,
                       std::numeric_limits<double>::epsilon() ) ;

    BOOST_CHECK_SMALL( std::fabs( distributionPointer->getProbabilityDensity(-1.0) - 0.0 ) ,
                       std::numeric_limits<double>::epsilon() ) ;

    BOOST_CHECK_SMALL( std::fabs( distributionPointer->getProbabilityDensity(3.1) - 0.0 ) ,
                       std::numeric_limits<double>::epsilon() ) ;

    BOOST_CHECK_SMALL( std::fabs( uniformDistribution.getMean() - 1.0 ) ,
                       std::numeric_limits<double>::epsilon() ) ;

    // CHECK STANDARD DEVIATION AND VARIANCE
}

//! Test if 1D Gaussian distribution class works correctly.
BOOST_AUTO_TEST_CASE( testGaussianDistributionDouble )
{
    // Compare PDF values with normpdf / normcdf function of MATLAB.
    double solutionMatlab;
    double solutionMatlabMass;
    double mean ;
    double standardDeviation;
    double x ;

    // Case 1
    {
        mean = 0.0 ;
        standardDeviation  = 1.0 ;
        x = 0.0 ;
        solutionMatlab = 3.989422804014327e-01 ;
        solutionMatlabMass = 0.5 ;

        tudat::statistics::GaussianDistributiond distribution(mean,standardDeviation);

        BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity(x) - solutionMatlab ) ,
                           std::numeric_limits<double>::epsilon() ) ;
        BOOST_CHECK_SMALL( std::fabs( distribution.getCumulativeProbability(x) - solutionMatlabMass ) ,
                           std::numeric_limits<double>::epsilon() ) ;
    }

    // Case 2
    {
        mean = 1.0 ;
        standardDeviation  = 2.0 ;
        x = -1.0 ;
        solutionMatlab = 1.209853622595717e-01 ;
        solutionMatlabMass = 1.586552539314571e-01 ;

        tudat::statistics::GaussianDistributiond distribution(mean,standardDeviation);

        BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity(x) - solutionMatlab ) ,
                           std::numeric_limits<double>::epsilon() ) ;
        BOOST_CHECK_SMALL( std::fabs( distribution.getCumulativeProbability(x) - solutionMatlabMass ) ,
                           std::numeric_limits<double>::epsilon() ) ;
    }

    // Case 3
    {
        mean = 3.0 ;
        standardDeviation  = 6.0 ;
        x = 2.0 ;
        solutionMatlab = 6.557328601698999e-02 ;
        solutionMatlabMass = 4.338161673890963e-01 ;

        tudat::statistics::GaussianDistributiond distribution(mean,standardDeviation);

        BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity(x) - solutionMatlab ) ,
                           std::numeric_limits<double>::epsilon() ) ;
        BOOST_CHECK_SMALL( std::fabs( distribution.getCumulativeProbability(x) - solutionMatlabMass ) ,
                           std::numeric_limits<double>::epsilon() ) ;
    }



}

//! Test if ND Gaussian distribution class works correctly.
BOOST_AUTO_TEST_CASE( testGaussianDistributionMultiD )
{
    using namespace tudat::statistics;
    using tudat::mathematical_constants::PI;
    using std::pow;
    using std::exp;

    Eigen::VectorXd mean(2);
    mean << 0.0 , 1.0 ;
    Eigen::MatrixXd covariance(2,2);
    covariance <<   3.0 , -1.0 ,
                    -1.0 , 3.0 ;
    GaussianDistributionXd distribution(mean,covariance);

    Eigen::VectorXd x(2);
    x << 0.0 , 1.0 ;

    double sigmax = std::sqrt(covariance(0,0));
    double sigmay = std::sqrt(covariance(1,1));
    double mux = mean(0);
    double muy = mean(1);
    double rho = covariance(0,1) / (sigmax*sigmay) ;

    // Exact solution of 2D Gaussian distribution
    // Montgomery, D. C. & Runger, G. C. Applied Statistics and Probability for engineers Wiley, 2014
    double exactSolution = (1.0/( 2.0*PI*sigmax*sigmay*sqrt(1.0-pow(rho,2.0)) )) *
                              exp( (-1.0/(2.0*(1.0-pow(rho,2.0)) ))*(
                                       (pow(x(0)-mux,2.0))/(pow(sigmax,2.0))
                                       - ( 2.0*rho*(x(0)-mux)*(x(1)-muy) )/(sigmax*sigmay)
                                       + (pow(x(1)-muy,2.0))/(pow(sigmay,2.0))
                                       )  ) ;

    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity(x) - exactSolution ) ,
                       std::numeric_limits<double>::epsilon() ) ;

    // Second test
    mean << -1.0 , 2.0 ;
    covariance <<   4.0 , 1.5 ,
                    1.5 , 4.0 ;
    x << 2.5 , 0.0 ;

    sigmax = std::sqrt(covariance(0,0));
    sigmay = std::sqrt(covariance(1,1));
    mux = mean(0);
    muy = mean(1);
    rho = covariance(0,1) / (sigmax*sigmay) ;

    GaussianDistributionXd distribution2(mean,covariance);

    // Exact solution of 2D Gaussian distribution
    // Montgomery, D. C. & Runger, G. C. Applied Statistics and Probability for engineers Wiley, 2014
    exactSolution = (1.0/( 2.0*PI*sigmax*sigmay*sqrt(1.0-pow(rho,2.0)) )) *
                              exp( (-1.0/(2.0*(1.0-pow(rho,2.0)) ))*(
                                       (pow(x(0)-mux,2.0))/(pow(sigmax,2.0))
                                       - ( 2.0*rho*(x(0)-mux)*(x(1)-muy) )/(sigmax*sigmay)
                                       + (pow(x(1)-muy,2.0))/(pow(sigmay,2.0))
                                       )  ) ;

    BOOST_CHECK_SMALL( std::fabs( distribution2.getProbabilityDensity(x) - exactSolution ) ,
                       std::numeric_limits<double>::epsilon() ) ;
}

//! Test if Gaussian Copula distribution class works correctly.
BOOST_AUTO_TEST_CASE( testGaussianCopula )
{
    using namespace tudat::statistics;
    using tudat::mathematical_constants::PI;

    int dimension = 2 ;
    Eigen::MatrixXd correlationMatrix( dimension , dimension ) ;
    correlationMatrix << 1.0 , 0.3 ,
                         0.3 , 1.0 ;
    GaussianCopulaDistributionXd distribution( dimension , correlationMatrix );

    Eigen::VectorXd x(2);
    x << 0.5 , 0.3 ;

//    std::cout << distribution.getProbabilityDensity( x ) << std::endl;

    // Test out of bounds
    x << 1.1 , 0.5 ;
    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( x ) - 0.0 ) , 1E-15 );

    x << 0.1 , -0.5 ;
    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( x ) - 0.0 ) , 1E-15 );

    x << -0.1 , -0.5 ;
    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( x ) - 0.0 ) , 1E-15 );

    x << 0.1 , 1.5 ;
    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( x ) - 0.0 ) , 1E-15 );

    // ADD MORE TESTS!!!!!!!!!!!

}

//! Test if Log normal distribution class works correctly.
BOOST_AUTO_TEST_CASE( test_log_normal_distribution )
{

    tudat::statistics::LogNormalDistributiond distribution( 2.0 , 0.5 );

    double locationParameter = distribution.getLocationParameter();
    double scaleParameter = distribution.getScaleParameter();

    BOOST_CHECK_SMALL( std::fabs( locationParameter - 6.628348696517279e-01 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( scaleParameter - 2.462206770692398e-01 ) , 1E-15 );

    boost::math::lognormal distributionBoost( locationParameter , scaleParameter );

    double computedProbabilityDensity = distribution.getProbabilityDensity( 1.0 );

    // Computed using MATLAB
    double expectedProbabilityDensity = 4.324214246844349e-02 ;
//    double expectedProbabilityDensity = boost::math::pdf( distributionBoost , 1.0 ) ;

    BOOST_CHECK_SMALL( std::fabs( computedProbabilityDensity - expectedProbabilityDensity ) , 1.0E-15 );

    double computedCumulativeProbability = distribution.getCumulativeProbability( 1.0 );

    // Computed using MATLAB
    double expectedCumulativeProbability = 3.550866413284351e-03 ;

    BOOST_CHECK_SMALL( std::fabs( computedCumulativeProbability - expectedCumulativeProbability ) , 1.0E-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
