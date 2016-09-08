#ifndef TUDAT_RANDOM_SAMPLING_H
#define TUDAT_RANDOM_SAMPLING_H

#include <iostream>
#include <map>

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"

namespace tudat
{

namespace statistics
{

//! Generator random vector using pseudo random generator
std::vector< Eigen::VectorXd > generateRandomVectorUniform(
        const int seed, const int numberOfSamples,
        const Eigen::VectorXd& lowerBound, const Eigen::VectorXd& upperBound);

//! Generator random vector using pseudo random generator
std::vector< Eigen::VectorXd > generateRandomVectorUniform(
        const int seed, const int numberOfSamples, const int numberOfDimensions,
        const double lowerBound = 0.0, const double upperBound = 1.0 );

//! Generator random vector using pseudo random generator with gaussian distribution (without correlation)
std::vector< Eigen::VectorXd > generateRandomVectorGaussian(
        const int seed, const int numberOfSamples, const int numberOfDimensions,
        const double mean, const double standardDeviation );

//! Generator random vector using pseudo random generator with gaussian distribution (without correlation)
std::vector< Eigen::VectorXd > generateRandomVectorGaussian(
        const int seed, const int numberOfSamples,
        const Eigen::VectorXd& mean, const Eigen::VectorXd& standardDeviation );

std::vector< Eigen::VectorXd > generateRandomVector(
        const int numberOfSamples,
        const std::vector< boost::shared_ptr< ContinuousRandomVariableGenerator > > randomVariableGenerators );

std::vector< Eigen::VectorXd > generateRandomVector(
        const int numberOfSamples, const int numberOfDimensions,
        const boost::shared_ptr< ContinuousRandomVariableGenerator > randomVariableGenerator );

#if USE_GSL
//! Generator random vector using Sobol sampler
std::vector< Eigen::VectorXd > sobolSamplerXd(const int Dimension, int numberOfSamples,
                                              Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound);

//! Generator random vector using Sobol sampler
std::vector< Eigen::VectorXd > sobolSamplerXd(const int Dimension, int numberOfSamples);
#endif

} // Close Namespace statistics

} // Close Namespace tudat

#endif // TUDAT_RANDOM_SAMPLING_H
