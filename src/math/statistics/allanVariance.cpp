#include "tudat/math/statistics/allanVariance.h"

namespace tudat
{

namespace statistics
{

//! Function to calculate the allan variance of a set of time offsets for a single bin size.
double calculateAllanVarianceOfTimeDataSetAtInterval(
        const int binSize, const double singleTimeStepSize, const std::vector< double >& timeErrors )
{
    // Determine time duration of a single bin.
    double binTimeStep = singleTimeStepSize * static_cast< double >( binSize );

    // Define start index.
    int startIndex = 0;

    // Declare variables for iteration
    double allanVariance = 0.0;
    double frequencyChange = 0.0;
    int numberOfSteps = 0;

    // Iterate over all bins
    while( startIndex + 2 * binSize < static_cast< int >( timeErrors.size( ) ) )
    {
        // Calculate frequency change (multiplied by time interval squared; division by time interval squared is one just before return)
        frequencyChange = ( timeErrors[ startIndex + 2 * binSize ] - 2.0 * timeErrors[ startIndex + binSize ] + timeErrors[ startIndex ] );

        // Add contribution from current bin
        allanVariance += frequencyChange * frequencyChange;

        // Update variables for next iteration.
        startIndex += binSize;
        numberOfSteps++;
    }

    // Normalize Allan variance and return
    return allanVariance / ( ( numberOfSteps - 1 ) * ( 2.0 * binTimeStep * binTimeStep ) );
}

//! Function to calculate the Allan variance at all time intervals from equi-spaced vector of (timing) errors in time domain.
std::map< double, double > calculateAllanVarianceOfTimeDataSet(
        const std::vector< double >& timingErrors, const double timeStepSize )
{
    // Define initial bin size
    int currentBinSize = 1;

    // Declare map of allan variances.
    std::map< double, double > allanVariances;

    // Increase bin size by factor of 2 until maximum size reached
    while( currentBinSize < static_cast< int >( timingErrors.size( ) ) / 2 )
    {
        // Calculate allan variance at current bin size
        allanVariances[ timeStepSize * static_cast< double >( currentBinSize ) ] =
                calculateAllanVarianceOfTimeDataSetAtInterval( currentBinSize, timeStepSize, timingErrors );

        // Increase bin size
        currentBinSize *= 2;
    }

    return allanVariances;
}

}

}
