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

#ifndef TUDAT_READAMPLITUDEANDARGUMENTMULTIPLIERS_H
#define TUDAT_READAMPLITUDEANDARGUMENTMULTIPLIERS_H

#include <cmath>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace earth_orientation
{

//! Function to read the amplitudes and fundamental argument multipliers for tidal corrections
/*!
 *  Function to read the amplitudes and fundamental argument multipliers for tidal corrections to
 *  earth orientation etc. Two file names are required, one containing the fundamental argument multipliers
 *  (one set per row) and one file containing the amplitudes (number of entries dependent on
 *  correction that is calculated). Function can filter amplitudes, so that only those entries
 *  with an RSS higher than a certain value are included in the return type.
 *  \param amplitudesFile File name of file containing correction amplitudes.
 *  \param fundamentalArgumentMultipliersFile File name of file containing fundamental argument multipliers for corrections.
 *  \param minimumAmplitude Minimum value of RSS of amplitudes for single entry for which entry is
 *  accepted.
 *  \return Pair of Matrices with 1: fundamental argument multipliers and 2: Amplitudes. A single entry is
 *  stored on a single row with same index for  fundamental argument multipliers and amplitude.
 */
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > readAmplitudesAndFundamentalArgumentMultipliers(
        const std::string amplitudesFile,
        const std::string fundamentalArgumentMultipliersFile,
        const double minimumAmplitude = 0.0 );

//! Function to filter tidal corrections based on RSS amplitude criteria.
/*!
 *  Function to filter tidal corrections based on RSS amplitude criteria.
 *  Only those entries with an RSS higher than a provided value are returned.
 *  \param rawAmplitudes Matrix of amplitudes, one correction set per row.
 *  \param rawFundamentalArgumentMultipliers Matrix of fundamental argument multipliers, one correction set per row.
 *  \param minimumAmplitude Minimum value of RSS of amplitudes for single entry for which entry is
 *  accepted.
 *  \return Pair of Matrices with 1: fundamental argument multipliers and 2: Amplitudes for which RSS amplitude
 *   is sufficient. A single entry is stored on a single row with same index for fundamental argument
 *   multipliers and amplitude.
 */
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > filterRawDataForAmplitudes(
        const Eigen::MatrixXd rawAmplitudes,
        const Eigen::MatrixXd rawFundamentalArgumentMultipliers,
        const double minimumAmplitude );

}

}

#endif // TUDAT_READAMPLITUDEANDARGUMENTMULTIPLIERS_H
