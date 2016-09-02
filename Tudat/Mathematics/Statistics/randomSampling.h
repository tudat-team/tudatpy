#ifndef TUDAT_RANDOM_SAMPLING_H
#define TUDAT_RANDOM_SAMPLING_H

#include <iostream> // cout sometimes needs this
#include <map>

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

namespace tudat
{

namespace statistics
{

//! Generator random vector using pseudo random generator
std::vector< Eigen::VectorXd > generateRandomVectorUniform(int seed, int numberOfSamples, int Dimension);

//! Generator random vector using pseudo random generator
std::vector< Eigen::VectorXd > generateRandomVectorUniform(int seed, int numberOfSamples,
                                             Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound);

////! Generator random vector using Sobol sampler
//std::vector< Eigen::VectorXd > sobolSamplerXd(const int Dimension, int numberOfSamples,
//                                              Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound);

////! Generator random vector using Sobol sampler
//std::vector< Eigen::VectorXd > sobolSamplerXd(const int Dimension, int numberOfSamples);

//! Generator random vector using pseudo random generator with gaussian distribution (without correlation)
std::vector< Eigen::VectorXd > generateRandomVectorGaussian(int seed, int numberOfSamples, int Dimension,
                                    Eigen::VectorXd mean , Eigen::VectorXd standardDeviation );

} // Close Namespace statistics

} // Close Namespace tudat

#endif // TUDAT_RANDOM_SAMPLING_H
