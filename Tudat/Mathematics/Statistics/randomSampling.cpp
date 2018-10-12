/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cmath>

#include <boost/random.hpp>
#include <boost/make_shared.hpp>

#if USE_GSL
#include <gsl/gsl_qrng.h>
#endif

#include "Tudat/Mathematics/Statistics/randomSampling.h"
namespace tudat
{

namespace statistics
{

//! Generate sample of random vectors, with entries of each vector independently, but not identically, distributed.
std::vector< Eigen::VectorXd > generateRandomSampleFromGenerator(
        const int numberOfSamples,
        const std::vector< std::shared_ptr< RandomVariableGenerator< double > > > randomVariableGenerators )
{
    std::vector< Eigen::VectorXd > randomSamples;
    Eigen::VectorXd randomSample( randomVariableGenerators.size( ) );

    // Generate samples
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

//! Generate sample of random vectors, with entries of each vector independently and identically distributed.
std::vector< Eigen::VectorXd > generateRandomSampleFromGenerator(
        const int numberOfSamples, const int numberOfDimensions,
        const std::shared_ptr< RandomVariableGenerator< double > > randomVariableGenerator )
{
    std::vector< std::shared_ptr< RandomVariableGenerator< double > > > randomVariableGenerators;
    for( int i = 0; i < numberOfDimensions; i++ )
    {
        randomVariableGenerators.push_back( randomVariableGenerator );
    }

    return generateRandomSampleFromGenerator( numberOfSamples, randomVariableGenerators );
}



//! Generator random vector using pseudo random generator
std::vector< Eigen::VectorXd > generateUniformRandomSample(
        const int seed, const int numberOfSamples,
        const Eigen::VectorXd& lowerBound, const Eigen::VectorXd& upperBound )
{
    if( lowerBound.rows( ) != upperBound.rows( ) )
    {
        throw std::runtime_error( "Error when making uniformly distributed samples, input is inconsistent" );
    }

    // Create distributions
    std::vector< std::shared_ptr< RandomVariableGenerator< double > > > randomVariableGenerators;
    std::vector< double > currentParameters;
    for( int i = 0; i < lowerBound.rows( ); i++ )
    {
        currentParameters = { lowerBound( i ), upperBound( i ) };
        randomVariableGenerators.push_back(
                    createBoostContinuousRandomVariableGenerator(
                        uniform_boost_distribution, currentParameters, seed + i ) );
    }

    // Generate samples
    return generateRandomSampleFromGenerator( numberOfSamples, randomVariableGenerators );
}

//! Generator random vector using pseudo random generator
std::vector< Eigen::VectorXd > generateUniformRandomSample(
        const int seed, const int numberOfSamples, const int numberOfDimensions,
         const double lowerBound, const double upperBound )
{
    return generateUniformRandomSample(
                seed, numberOfSamples,
                Eigen::VectorXd::Constant( numberOfDimensions, lowerBound ),
                Eigen::VectorXd::Constant( numberOfDimensions, upperBound ) );
}



//! Generator random vector using pseudo random generator with gaussian distribution (without correlation)
std::vector< Eigen::VectorXd > generateGaussianRandomSample(
        const int seed, const int numberOfSamples,
        const Eigen::VectorXd& mean, const Eigen::VectorXd& standardDeviation )
{
    if( mean.rows( ) != standardDeviation.rows( ) )
    {
        throw std::runtime_error( "Error when making Gaussian distributed samples, input is inconsistent" );
    }

    // Create distributions
    std::vector< std::shared_ptr< RandomVariableGenerator< double > > > randomVariableGenerators;
    std::vector< double > currentParameters;
    for( int i = 0; i < mean.rows( ); i++ )
    {
        currentParameters = { mean( i ), standardDeviation( i ) };
        randomVariableGenerators.push_back(
                    createBoostContinuousRandomVariableGenerator(
                        normal_boost_distribution, currentParameters, seed + i ) );
    }

    // Generate samples
    return generateRandomSampleFromGenerator( numberOfSamples, randomVariableGenerators );
}


//! Generator random vector using pseudo random generator with gaussian distribution (without correlation)
std::vector< Eigen::VectorXd > generateGaussianRandomSample(
        const int seed, const int numberOfSamples, const int numberOfDimensions,
        const double mean, const double standardDeviation )
{
    return generateGaussianRandomSample(
                seed, numberOfSamples,
                Eigen::VectorXd::Constant( numberOfDimensions, mean ),
                Eigen::VectorXd::Constant( numberOfDimensions, standardDeviation ) );
}


#if USE_GSL

//! Generator random vector using Sobol sampler
std::vector< Eigen::VectorXd > generateVectorSobolSample(
        int numberOfSamples,
        const Eigen::VectorXd& lowerBound, const Eigen::VectorXd& upperBound )
{
    int numberOfDimensions = upperBound.rows( );

    // Compute propertie
    Eigen::VectorXd width = upperBound - lowerBound;
    Eigen::VectorXd average = (upperBound + lowerBound)/2.0 ;

    std::vector< Eigen::VectorXd > sobolSamples( numberOfSamples );

    Eigen::VectorXd randomSample( numberOfDimensions ) ;
    double randomSampleArray[ numberOfDimensions ];

    gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, numberOfDimensions );

    // Loop over samples
    for( int j = 0 ; j < numberOfSamples ; j++ )
    {
        gsl_qrng_get( q, randomSampleArray ); // Generate sobol [0,1]

        // Fill vector
        for( int i = 0 ; i < numberOfDimensions ; i++ )
        {
            randomSample( i ) = randomSampleArray[ i ] - 0.5 ;
        }
        sobolSamples[ j ] = randomSample.cwiseProduct( width ) + average ; // Save vector
    }

    gsl_qrng_free (q); // Deallocate GSL variables

    return sobolSamples;
}

//! Generator random vector using Sobol sampler
std::vector< Eigen::VectorXd > generateVectorSobolSample(
        const int numberOfDimensions, int numberOfSamples,
        const double lowerBound, const double upperBound )
{
    return generateVectorSobolSample( numberOfSamples,
                                      Eigen::VectorXd::Constant( numberOfDimensions, lowerBound ),
                                      Eigen::VectorXd::Constant( numberOfDimensions, upperBound ) );
}
#endif

} // Close Namespace statistics

} // Close Namespace tudat

