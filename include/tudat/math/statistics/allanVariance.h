#ifndef ALLANVARIANCE_H
#define ALLANVARIANCE_H

#include <map>
#include <vector>

namespace tudat
{

namespace statistics
{

//! Function to calculate the Allan variance at a given time interval from equi-spaced vector of (timing) errors in time domain.
/*!
 *  Function to calculate the Allan variance at a given time interval from equi-space vector of (timing) errors in time domain. The Allan variance
 *  is calculated at a time interval of an integer multiple of the time interval between two entries of the input vector.
 *  \param binSize Number of entries over input vector timeErrors that define interval for which Allan variance is to be calculated.
 *  \param singleTimeStepSize Time interval between two entries in timeErrors
 *  \param timeErrors Timing errors, equispaced at singleTimeStepSize between consecutive entries.
 *  \return Allan variance at specified bin size
 */
double calculateAllanVarianceOfTimeDataSetAtInterval(
        const int binSize, const double singleTimeStepSize, const std::vector< double >& timeErrors );

//! Function to calculate the Allan variance at all time intervals from equi-spaced vector of (timing) errors in time domain.
/*!
 *  Function to calculate the Allan variances at all time intervals from equi-space vector of (timing) errors in time domain. The Allan variance
 *  can be calculated at a time interval of an integer multiple of the time interval between two entries of the input vector.
 *  \param timingErrors Timing errors, equispaced at singleTimeStepSize between consecutive entries.
 *  \param timeStepSize Time interval between two entries in timeErrors
 *  \return Map of key: time interval and value: Allan variance at given time interval.
 */
std::map< double, double > calculateAllanVarianceOfTimeDataSet(
        const std::vector< double >& timingErrors, const double timeStepSize );

}

}

#endif // ALLANVARIANCE_H
