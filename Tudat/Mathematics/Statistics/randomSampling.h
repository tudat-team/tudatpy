/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RANDOM_SAMPLING_H
#define TUDAT_RANDOM_SAMPLING_H

#include <map>

#include <Eigen/Core>

#include <memory>

#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"

namespace tudat
{

namespace statistics
{

//! Generate sample of random vectors, with entries of each vector independently, but not identically, distributed.
/*!
 *  Function to generate sample of random vectors, with entries of each vector independently, but not identically,
 *  distributed according to probability distributions provided as input. Size of each sample is define by size of
 *  vector randomVariableGenerators.
 *  \param numberOfSamples Number of samples that are to be generated.
 *  \param randomVariableGenerators Probability distributions for the entries of the random vectors (i.e. entry i of
 *  this vector is distribution of entry i of each sample).
 *  \return Set of samples according to given distribution settings.
 */
std::vector< Eigen::VectorXd > generateRandomSampleFromGenerator(
        const int numberOfSamples,
        const std::vector< std::shared_ptr< RandomVariableGenerator< double > > > randomVariableGenerators );

//! Generate sample of random vectors, with entries of each vector independently and identically distributed.
/*!
 *  Function to generate sample of random vectors, with entries of each vector independently and identically,
 *  distributed according to probability distribution provided as input.
 *  \param numberOfSamples Number of samples that are to be generated.
 *  \param numberOfDimensions Size of each sample.
 *  \param randomVariableGenerator Probability distribution used for each of the entries of the random vectors.
 *  \return Set of samples according to given distribution settings.
 */
std::vector< Eigen::VectorXd > generateRandomSampleFromGenerator(
        const int numberOfSamples, const int numberOfDimensions,
        const std::shared_ptr< RandomVariableGenerator< double > > randomVariableGenerator );



//! Generate sample of random vectors, with entries of each vector independently, but not identically, uniformly distributed.
/*!
 *  Function to generate sample of random vectors, with entries of each vector independently, but not identically,
 *  uniformly distributed. The size of each sample is defined by the size of the lowerBound and upperBound vectors
 *  (which must have identical size).
 *  \param seed Seed of random number generator.
 *  \param numberOfSamples Number of samples that are to be generated.
 *  \param lowerBound Vector of lower bounds for the distributions for the entries of the random vectors (i.e. entry i of
 *  this vector is lower bound for distribution of entry i of each sample).
 *  \param upperBound Vector of upper bounds for the distributions for the entries of the random vectors (i.e. entry i of
 *  this vector is upper bound for distribution of entry i of each sample).
 *  \return Set of samples according to given distribution settings.
 */
std::vector< Eigen::VectorXd > generateUniformRandomSample(
        const int seed, const int numberOfSamples,
        const Eigen::VectorXd& lowerBound, const Eigen::VectorXd& upperBound);

//! Generate sample of random vectors, with entries of each vector independently and identically uniformly distributed.
/*!
 *  Function to generate sample of random vectors, with entries of each vector independently and identically
 *  uniformly distributed.
 *  \param seed Seed of random number generator.
 *  \param numberOfSamples Number of samples that are to be generated.
 *  \param numberOfDimensions Size of each sample.
 *  \param lowerBound Lower bound for the distributions for the entries of the random vectors
 *  \param upperBound Upper bound for the distributions for the entries of the random vectors
 *  \return Set of samples according to given distribution settings.
 */
std::vector< Eigen::VectorXd > generateUniformRandomSample(
        const int seed, const int numberOfSamples, const int numberOfDimensions,
        const double lowerBound = 0.0, const double upperBound = 1.0 );



//! Generate sample of random vectors, with entries of each vector independently, but not identically, Gaussian distributed.
/*!
 *  Function to generate sample of random vectors, with entries of each vector independently, but not identically,
 *  Gaussian distributed. The size of each sample is defined by the size of the lowerBound and upperBound vectors
 *  (which must have identical size).
 *  \param seed Seed of random number generator.
 *  \param numberOfSamples Number of samples that are to be generated.
 *  \param mean Vector of mean values for the distributions for the entries of the random vectors (i.e. entry i of
 *  this vector is mean for distribution of entry i of each sample).
 *  \param standardDeviation Vector of standard deviations for the distributions for the entries of the random vectors
 *  (i.e. entry i of this vector is standard deviation for distribution of entry i of each sample). *
 *  \return Set of samples according to given distribution settings.
 */
std::vector< Eigen::VectorXd > generateGaussianRandomSample(
        const int seed, const int numberOfSamples,
        const Eigen::VectorXd& mean , const Eigen::VectorXd& standardDeviation );

//! Generate sample of random vectors, with entries of each vector independently and identically Gaussian distributed.
/*!
  *  Function to generate sample of random vectors, with entries of each vector independently and identically
  *  Gaussian distributed.
  *  \param seed Seed of random number generator.
  *  \param numberOfSamples Number of samples that are to be generated.
  *  \param numberOfDimensions Size of each samples
  *  \param mean Mean value for the distributions for the entries of the random vectors
  *  \param standardDeviation Standard deviation for the distributions for the entries of the random vectors
  *  \return Set of samples according to given distribution settings.
  */
std::vector< Eigen::VectorXd > generateGaussianRandomSample(
        const int seed, const int numberOfSamples, const int numberOfDimensions,
        const double mean = 0.0, const double standardDeviation = 1.0 );



#if USE_GSL

//! Generate sample of random vectors, using a Sobol sampling algorithm.
/*!
 *  Generate sample of random vectors, using a Sobol sampling algorithm. The size of each sample is defined by the size of
 *  the lowerBound and upperBound vectors (which must have identical size).
 *  \param numberOfSamples Number of samples that are to be generated.
 *  \param lowerBound Vector of lower bounds for the entries of the samples (i.e. entry i of
 *  this vector is lower bound for distribution of entry i of each sample).
 *  \param upperBound Vector of upper bounds for the entries of the samples(i.e. entry i of
 *  this vector is upper bound for distribution of entry i of each sample).
 *  \return Set of samples generated with Sobol algorithm
 */
std::vector< Eigen::VectorXd > generateVectorSobolSample(
        int numberOfSamples,
        const Eigen::VectorXd& lowerBound, const Eigen::VectorXd& upperBound );

//! Generate sample of random vectors, using a Sobol sampling algorithm.
/*!
 *  Generate sample of random vectors, using a Sobol sampling algorithm.
 *  \param numberOfSamples Number of samples that are to be generated.
 *  \param numberOfDimensions Size of each sample.
 *  \param lowerBound Lower bound for the distributions for the entries of the random vectors
 *  \param upperBound Upper bound for the distributions for the entries of the random vectors
 *  \return Set of samples generated with Sobol algorithm
 */
std::vector< Eigen::VectorXd > generateVectorSobolSample(
        const int numberOfDimensions, int numberOfSamples,
        const double lowerBound = 0.0, const double upperBound = 1.0 );
#endif

} // Close Namespace statistics

} // Close Namespace tudat

#endif // TUDAT_RANDOM_SAMPLING_H
