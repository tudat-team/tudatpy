#ifndef READAMPLITUDEANDDOODSONNUMBER_H
#define READAMPLITUDEANDDOODSONNUMBER_H

#include <cmath>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace earth_orientation
{

//! Function to read the amplitudes and doodson multiplers for tidal corrections
/*!
 *  Function to read the amplitudes and doodson multiplers for tidal corrections to
 *  earth orientation etc. Two file names are required, one containing the Doodson numbers
 *  (one set per row) and one file containing the amplitudes (number of entries dependent on
 *  correction that is calculated). Function can filter amplitudes, so that only those entries
 *  with an RSS higher than a certain value are included in the return type.
 *  \param amplitudesFile File name of file containing correction amplitudes.
 *  \param doodsonMultipliersFile File name of file containing doodson numbers for corrections.
 *  \param numberOfAmplitudes Number of amplitudes that are to be read per entry.
 *  \param minimumAmplitude Minimum value of RSS of amplitudes for single entry for which entry is
 *  accepted.
 *  \return Pair of Matrices with 1: Doodson multipliers and 2: Amplitudes. A single entry is
 *  stored on a single row with same index for  Doodson multipliers and amplitude.
 */
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > readAmplitudesAndDoodsonMultipliers(
        const std::string amplitudesFile,
        const std::string DoodsonMultipliersFile,
        const double minimumAmplitude = 0.0 );

//! Function to filter tidal corrections based on RSS amplitude criteria.
/*!
 *  Function to filter tidal corrections based on RSS amplitude criteria.
 *  Only those entries with an RSS higher than a provided value are returned.
 *  \param rawAmplitudes Matrix of amplitudes, one correction set per row.
 *  \param rawDoodsonMultipliers Matrix of Doodson multiplers, one correction set per row.
 *  \param numberOfAmplitudes Number of amplitudes that are to be read per entry.
 *  \param minimumAmplitude Minimum value of RSS of amplitudes for single entry for which entry is
 *  accepted.
 *  \return Pair of Matrices with 1: Doodson multipliers and 2: Amplitudes for which RSS amplitude
 *   is sufficient. A single entry is stored on a single row with same index for Doodson
 *   multipliers and amplitude.
 */
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > filterRawDataForAmplitudes(
        const Eigen::MatrixXd rawAmplitudes,
        const Eigen::MatrixXd rawDoodsonMultipliers,
        const double minimumAmplitude );

}

}

#endif // READAMPLITUDEANDDOODSONNUMBER_H
