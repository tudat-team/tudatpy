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
 *      160902    R. Hoogendoorn    File created.
 *
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <iostream> // cout sometimes needs this
#include <vector>
#include <limits>

#include <Eigen/Core>


#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <boost/test/unit_test.hpp>
#include "tudat/Basics/testMacros.h"

#include "tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "tudat/Mathematics/Statistics/randomSampling.h"
#include <tudat/Mathematics/Statistics/basicStatistics.h>

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_Random_Sampler )

//using tudat::mathematical_constants::PI;

//BOOST_AUTO_TEST_CASE( test_Sobol_Sampler )
//{
//    int dimension = 3 ;
//    int numberOfSamples = 1E6 ;

//    Eigen::VectorXd lower(dimension);
//    Eigen::VectorXd upper(dimension);
//    lower << 0.0 , 1.0 , -2.0 ;
//    upper << 1.0 , 3.0 , 4.0 ;

//    Eigen::VectorXd width = upper - lower ;
//    Eigen::VectorXd average = (upper + lower) / 2.0 ;

//    std::vector<Eigen::VectorXd> sobolSamples = Thesis::Statistics::sobolSamplerXd(dimension,numberOfSamples,lower,upper);

//    // Compute sample average
//    Eigen::VectorXd sampleAverage(dimension);
//    sampleAverage.setZero(dimension,1);
//    for( unsigned int i = 0 ; i < sobolSamples.size() ; i++){
//        sampleAverage += sobolSamples[i] ;
//    }
//    sampleAverage = sampleAverage / double(numberOfSamples) ;

//    BOOST_CHECK_SMALL( std::fabs( average(0) - sampleAverage(0) ) , 2.0E-6 );
//    BOOST_CHECK_SMALL( std::fabs( average(1) - sampleAverage(1) ) , 2.0E-6 );
//    BOOST_CHECK_SMALL( std::fabs( average(2) - sampleAverage(2) ) , 2.0E-6 );

//    sobolSamples = Thesis::Statistics::sobolSamplerXd(dimension,numberOfSamples);

//    // Compute sample average
//    sampleAverage.setZero(dimension,1);
//    for( unsigned int i = 0 ; i < sobolSamples.size() ; i++){
//        sampleAverage += sobolSamples[i] ;
//    }
//    sampleAverage = sampleAverage / double(numberOfSamples) ;

//    BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleAverage(0) ) , 2.0E-6 );
//    BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleAverage(1) ) , 2.0E-6 );
//    BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleAverage(2) ) , 2.0E-6 );
//}

BOOST_AUTO_TEST_CASE( test_randomVectorUniform )
{
    int dimension = 3 ;
    int numberOfSamples = 1E6 ;
    int seed = 511;

    Eigen::VectorXd lower(dimension);
    Eigen::VectorXd upper(dimension);
    lower << 0.0 , 1.0 , -2.0 ;
    upper << 1.0 , 3.0 , 4.0 ;

    Eigen::VectorXd width = upper - lower ;
    Eigen::VectorXd average = (upper + lower) / 2.0 ;

    std::vector<Eigen::VectorXd> samples =
            tudat::statistics::generateRandomVectorUniform( seed , numberOfSamples, lower, upper );

    // Compute sample average
    Eigen::VectorXd sampleAverage(dimension);
    sampleAverage.setZero(dimension,1);
    for( unsigned int i = 0 ; i < samples.size() ; i++){
        sampleAverage += samples[i] ;
    }
    sampleAverage = sampleAverage / double(numberOfSamples) ;

    BOOST_CHECK_SMALL( std::fabs( average(0) - sampleAverage(0) ) , 2.0E-3 );
    BOOST_CHECK_SMALL( std::fabs( average(1) - sampleAverage(1) ) , 2.0E-3 );
    BOOST_CHECK_SMALL( std::fabs( average(2) - sampleAverage(2) ) , 2.0E-3 );

    samples = tudat::statistics::generateRandomVectorUniform( seed , numberOfSamples, dimension );

    // Compute sample average
    sampleAverage.setZero(dimension,1);
    for( unsigned int i = 0 ; i < samples.size() ; i++){
        sampleAverage += samples[i] ;
    }
    sampleAverage = sampleAverage / double(numberOfSamples) ;

    BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleAverage(0) ) , 2.0E-4 );
    BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleAverage(1) ) , 2.0E-4 );
    BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleAverage(2) ) , 2.0E-4 );
}

BOOST_AUTO_TEST_CASE( test_randomVectorGaussian )
{
    int dimension = 3 ;
    int numberOfSamples = 1E6 ;
    int seed = 511;

    Eigen::VectorXd mean(dimension);
    Eigen::VectorXd standardDev(dimension);
    mean << 0.0 , 1.0 , -2.0 ;
    standardDev << 1.0 , 3.0 , 4.0 ;

    std::vector<Eigen::VectorXd> samples =
            tudat::statistics::generateRandomVectorGaussian( seed , numberOfSamples, dimension, mean, standardDev );

    // Compute sample average
    Eigen::VectorXd sampleAverage(dimension);
    sampleAverage.setZero(dimension,1);
    for( unsigned int i = 0 ; i < samples.size() ; i++){
        sampleAverage += samples[i] ;
    }
    sampleAverage = sampleAverage / double(numberOfSamples) ;

    BOOST_CHECK_SMALL( std::fabs( mean(0) - sampleAverage(0) ) , 2.0E-3 );
    BOOST_CHECK_SMALL( std::fabs( mean(1) - sampleAverage(1) ) , 2.0E-3 );
    BOOST_CHECK_SMALL( std::fabs( mean(2) - sampleAverage(2) ) , 2.0E-3 );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
