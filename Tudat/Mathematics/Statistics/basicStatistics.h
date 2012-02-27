/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      090807    J. Melman         First creation of code.
 *      100930    D. Dirkx          Modified to comply with Tudat standards.
 *      100930    J. Melman         Implemented namespace, minor comment changes.
 *      120202    K. Kumar          Moved functions from linearAlgebra.h.
 *
 *    References
 *
 */

#ifndef TUDAT_BASIC_STATISTICS_H
#define TUDAT_BASIC_STATISTICS_H

#include <Eigen/Core>
#include <vector>

namespace tudat
{
namespace mathematics
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

//! Compute sample variance.
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

} // namespace statistics
} // namespace mathematics
} // namespace tudat

#endif // TUDAT_BASIC_STATISTICS_H
