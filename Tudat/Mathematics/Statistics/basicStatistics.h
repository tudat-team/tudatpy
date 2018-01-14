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

#ifndef TUDAT_BASIC_STATISTICS_H
#define TUDAT_BASIC_STATISTICS_H

#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace statistics
{

//! Compute average of the components of a vector.
/*!
 * Computes the average (arithmetic mean) of the components of a vector.
 * \param vectorOfData Vector containing data to be averaged.
 * \return Average of data in vector.
 */
double computeAverageOfVectorComponents( const Eigen::VectorXd& vectorOfData );

//! Compute standard deviation of the components of a vector.
/*!
 * Computes the standard deviation of the components of a vector.
 * \param vectorOfData Vector containing data to be averaged.
 * \return Standard deviation of data in vector.
 */
double computeStandardDeviationOfVectorComponents( const Eigen::VectorXd& vectorOfData );

//! Compute sample mean.
/*!
 * Computes sample mean based on the following unbiased estimator
 * (Spiegel and Stephens, 2008):
 * \f[
 *      \mu_{s} = \frac{ \sum_{i=1}^{N} X_{i} } { N }
 * \f]
 * where \f$\mu_{s}\f$ is the unbiased estimate of the sample mean,
 * \f$ N \f$ is the number of samples, and \f$ X \f$ is the sample value.
 * \param sampleData Sample data.
 * \return Sample mean.
 */
double computeSampleMean( const std::vector< double >& sampleData );

//! Compute sample mean for a sample of VectorXd
/*!
 * Computes sample mean for a sample of VectorXd based on the following unbiased estimator
 * (Spiegel and Stephens, 2008):
 * \f[
 *      \mu_{s} = \frac{ \sum_{i=1}^{N} X_{i} } { N }
 * \f]
 * where \f$\mu_{s}\f$ is the unbiased estimate of the sample mean,
 * \f$ N \f$ is the number of samples, and \f$ X \f$ is the sample value.
 * \param sampleData Sample data.
 * \return Sample mean.
 */
Eigen::VectorXd computeSampleMean( const std::vector< Eigen::VectorXd >& sampleData );

//! Compute sample variance .
/*!
 * Computes sample variance based on the following unbiased estimator
 * (Spiegel and Stephens, 2008):
 * \f[
 *      s^{2}_{s} = \frac{ 1 }{ N - 1 } * \sum_{i=1}^{N} X_{i}
 *                  ( X_{i} - \bar{ X } )^{ 2 } )
 * \f]
 * where \f$ s^{2}_{s} \f$ is the unbiased estimate of the sample variance,
 * \f$ N \f$ is the number of samples, \f$ X \f$ is the sample value, and
 * \f$ \bar{ X } \f$ is the sample mean.
 * \param sampleData Map containing sample data.
 * \return Sample variance.
 */
double computeSampleVariance( const std::vector< double >& sampleData );

//! Compute Sample median
/*!
 * Computes the median of a set of data samples.
 * (Montgomery, D. C. & Runger, G. C. Applied Statistics and Probability for engineers Wiley, 2014)
 *
 * \param sampleData Map containing sample data.
 * \return Sample variance.
 */
double computeSampleMedian( std::vector< double > sampleData );

//! Compute sample variance for a sample of VectorXd.
/*!
 * Computes sample variance for a sample of VectorXd based on the following unbiased estimator
 * (Spiegel and Stephens, 2008):
 * \f[
 *      s^{2}_{s} = \frac{ 1 }{ N - 1 } * \sum_{i=1}^{N} X_{i}
 *                  ( X_{i} - \bar{ X } )^{ 2 } )
 * \f]
 * where \f$ s^{2}_{s} \f$ is the unbiased estimate of the sample variance,
 * \f$ N \f$ is the number of samples, \f$ X \f$ is the sample value, and
 * \f$ \bar{ X } \f$ is the sample mean.
 * \param sampleData Map containing sample data.
 * \return Sample variance.
 */
Eigen::VectorXd computeSampleVariance( const std::vector< Eigen::VectorXd >& sampleData );

} // namespace statistics
} // namespace tudat

#endif // TUDAT_BASIC_STATISTICS_H
