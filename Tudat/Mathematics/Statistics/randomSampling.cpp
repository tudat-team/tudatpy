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
 *      160902    R. Hoogendoorn    File created.
 *
 *
 *    References
 *
 *    Notes
 *
 */



#include <boost/make_shared.hpp>
#include <cmath>

#if USE_GSL
#include <gsl/gsl_qrng.h>   // GSL Quasi Random generators -> Sobol
#endif

// Pseudo random sampler
#include <boost/random.hpp> // Random generator

#include <Tudat/Mathematics/Statistics/randomSampling.h>

namespace tudat
{

namespace statistics
{

//! Generator random vector using pseudo random generator
std::vector<Eigen::VectorXd> generateRandomVectorUniform(
        const int seed, const int numberOfSamples,
        const Eigen::VectorXd& lowerBound, const Eigen::VectorXd& upperBound)
{
    // Compute properties
    Eigen::VectorXd width = upperBound - lowerBound;
    Eigen::VectorXd average = (upperBound + lowerBound)/2.0 ;

    // Setup Random generator
    typedef boost::mt19937 RandomGeneratorType; // Mersenne Twister
    RandomGeneratorType randomGenerator(seed);              // Create random generator

    boost::uniform_real<> uniformDistribution(0.0 , 1.0); //
    boost::variate_generator< RandomGeneratorType, boost::uniform_real<> >
            Dice(randomGenerator, uniformDistribution); // define random generator

    std::vector< Eigen::VectorXd > randomSamples(numberOfSamples) ;
    Eigen::VectorXd randomSample( lowerBound.rows() ) ;

    // Sample
    for(int i = 0 ; i < numberOfSamples ; i++ ) // Generate N samples
    {
        for(int j = 0 ; j < randomSample.rows() ; j++) // Generate vector of samples
        {
            randomSample(j) = Dice() - 0.5 ;
        }
        randomSamples[i] = randomSample.cwiseProduct(width) + average ;
    }
    return randomSamples;
}


//! Generator random vector using pseudo random generator
std::vector<Eigen::VectorXd> generateRandomVectorUniform(
        const int seed, const int numberOfSamples, const int numberOfDimensions,
         const double lowerBound, const double upperBound )
{
    return generateRandomVectorUniform(
                seed, numberOfSamples,
                Eigen::VectorXd::Constant( numberOfDimensions, lowerBound ),
                Eigen::VectorXd::Constant( numberOfDimensions, upperBound ) );
}


//! Generator random vector using pseudo random generator with gaussian distribution (without correlation)
std::vector< Eigen::VectorXd > generateRandomVectorGaussian(
        const int seed, const int numberOfSamples,
        const Eigen::VectorXd& mean, const Eigen::VectorXd& standardDeviation )
{
    // Setup Random generator
    typedef boost::mt19937 RandomGeneratorType; // Mersenne Twister
    RandomGeneratorType randomGenerator(seed);              // Create random generator

    boost::normal_distribution<> standardNormalDistribution(0.0 , 1.0); // standard normal
    boost::variate_generator< RandomGeneratorType, boost::normal_distribution<> >
            Dice(randomGenerator, standardNormalDistribution); // define random generator

    std::vector< Eigen::VectorXd > randomSamples(numberOfSamples) ;
    Eigen::VectorXd randomSample( mean.rows( ) ) ;

    // Sample
    for(int i = 0 ; i < numberOfSamples ; i++ ) // Generate N samples
    {
        for(int j = 0 ; j < randomSample.rows() ; j++) // Generate vector of samples
        {
            randomSample(j) = Dice()*standardDeviation(j)  + mean(j)  ;
        }
        randomSamples[i] = randomSample ;
    }
    return randomSamples;
}


//! Generator random vector using pseudo random generator with gaussian distribution (without correlation)
std::vector< Eigen::VectorXd > generateRandomVectorGaussian(
        const int seed, const int numberOfSamples, const int numberOfDimensions,
        const double mean, const double standardDeviation )
{
    return generateRandomVectorUniform(
                seed, numberOfSamples,
                Eigen::VectorXd::Constant( numberOfDimensions, mean ),
                Eigen::VectorXd::Constant( numberOfDimensions, standardDeviation ) );
}

std::vector< Eigen::VectorXd > generateRandomVector(
        const int numberOfSamples,
        const std::vector< boost::shared_ptr< ContinuousRandomVariableGenerator > > randomVariableGenerators )
{
    std::vector< Eigen::VectorXd > randomSamples( numberOfSamples );
    Eigen::VectorXd randomSample( randomVariableGenerators.size( ) );

    for( int i = 0; i < numberOfSamples; i++ )
    {
        for( unsigned int j = 0; j < randomVariableGenerators.size( ); j++ )
        {
            randomSample( j ) = randomVariableGenerators.at( j )->getRandomVariableValue( );
        }
        randomSamples.push_back( randomSample );
    }

    return randomSamples;
}

std::vector< Eigen::VectorXd > generateRandomVector(
        const int numberOfSamples, const int numberOfDimensions,
        const boost::shared_ptr< ContinuousRandomVariableGenerator > randomVariableGenerator )
{
    std::vector< boost::shared_ptr< ContinuousRandomVariableGenerator > > randomVariableGenerators;
    for( int i = 0; i < numberOfDimensions; i++ )
    {
        randomVariableGenerators.push_back( randomVariableGenerator );
    }

    return generateRandomVector( numberOfSamples, randomVariableGenerators );
}

#if USE_GSL

//! Generator random vector using Sobol sampler
std::vector< Eigen::VectorXd > sobolSamplerXd(const int Dimension, int numberOfSamples,
                                              Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound){
    // Compute properties
    Eigen::VectorXd width = upperBound - lowerBound;
    Eigen::VectorXd average = (upperBound + lowerBound)/2.0 ;

    std::vector< Eigen::VectorXd > sobolSamples( numberOfSamples );

    Eigen::VectorXd randomSample( Dimension ) ;
    double randomSampleArray[ Dimension ];

    gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, Dimension );

    for( int j = 0 ; j < numberOfSamples ; j++ )
    { // Loop over samples
        gsl_qrng_get (q, randomSampleArray); // Generate sobol [0,1]

        for(int i = 0 ; i < Dimension ; i++){ // Fill vector
            randomSample(i) = randomSampleArray[i] - 0.5 ;
        }
        sobolSamples[j] = randomSample.cwiseProduct(width) + average ; // Save vector
    }

    gsl_qrng_free (q); // don't know why?

    return sobolSamples;
}

//! Generator random vector using Sobol sampler
std::vector< Eigen::VectorXd > sobolSamplerXd(const int Dimension, int numberOfSamples){
    std::vector< Eigen::VectorXd > sobolSamples( numberOfSamples );

    Eigen::VectorXd randomSample( Dimension ) ;
    double randomSampleArray[ Dimension ];

    gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, Dimension );

    for(int j = 0 ; j < numberOfSamples ; j++)
    { // Loop over samples
        gsl_qrng_get (q, randomSampleArray); // Generate sobol [0,1]

        // Fill vector
        for(int i = 0 ; i < Dimension ; i++)
        {
            randomSample(i) = randomSampleArray[i] ;
        }

        sobolSamples[ j ] = randomSample; // Save vector
    }

    gsl_qrng_free (q); // don't know why?

    return sobolSamples;
}
#endif

} // Close Namespace statistics

} // Close Namespace tudat

