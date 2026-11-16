/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_READDSNAMCFILE_H
#define TUDAT_READDSNAMCFILE_H

#include <map>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace input_output
{

//! Single record of a DSN Advanced Media Calibration (AMC) path delay file (see DSN document DORS_002_020501)
struct DsnAmcRecord {
    //! UTC time, in seconds since J2000
    double utcTime;

    //! Azimuth of the Advanced Water Vapor Radiometer (AWVR) [deg]
    double azimuth;

    //! Elevation of the AWVR [deg]
    double elevation;

    //! Zenith wet path delay [m] (9.99 m denotes a missing retrieval)
    double zenithWetDelay;

    //! Zenith dry (static) path delay [m] (9.99 m denotes a missing retrieval)
    double zenithDryDelay;

    //! Cloud liquid [micrometer]
    double cloudLiquid;

    //! Identifier of retrieval algorithm (e.g. 005/006: statistical, 200: Bayesian), with leading zeros retained
    std::string algorithmId;

    //! Cloud flag: 0: liquid water < 20 micrometer, 1: 20-100 micrometer, 2: 100-200 micrometer, 3: > 200 micrometer
    int cloudFlag;

    //! AWVR flag: 0: all AWVR observables within acceptable range, 1: AWVR observable(s) outside acceptable range
    int wvrFlag;

    //! Microwave Temperature Profiler (MTP) flag: 0: all MTP observables within acceptable range, 1: MTP observable(s) outside
    //! acceptable range
    int mtpFlag;

    //! Surface meteorology flag: 0: all surface meteorology observables within acceptable range, 1: observable(s) outside
    //! acceptable range
    int surfaceFlag;

    //! Wind flag: 0: wind observables within acceptable range, 1: wind observable(s) outside acceptable range
    int windFlag;

    //! Path delay flag: 0: retrieved zenith path delay within acceptable range, 1: retrieved zenith path delay outside acceptable
    //! range, 2: elevation outside prescribed range (no retrieval for Bayesian algorithm; retrieval possibly erroneous due to ground
    //! pickup for statistical algorithm)
    int pathDelayFlag;

    //! Wind speed [mph] (NaN if not in file)
    double windSpeed;

    //! Wind direction [deg] (NaN if not in file)
    double windDirection;
};

//! Function to read all records of a DSN AMC path delay file (without any quality screening)
/*!
 * Function to read all records of a DSN AMC path delay file (without any quality screening). The time is given in UTC seconds since
 * J2000, and path delays in m. Other quantities are in the units of the file. Both 16-column files (e.g. Juno) and 14-column files
 * without wind speed and direction (e.g. Cassini) are accepted, where the columns include the leading 'P'. For 14-column files, the
 * wind speed and direction are set to NaN.
 * \param fileName Name of AMC file
 * \return List of records in file
 */
std::vector< DsnAmcRecord > readDsnAmcFile( const std::string& fileName );

//! Function to retrieve valid zenith path delays from DSN AMC records
/*!
 * Function to retrieve valid zenith path delays from DSN AMC records. Records with an invalid wet path delay or bad quality flags are
 * rejected, where the input arguments define part of the rejection criteria. An invalid dry path delay is stored as NaN.
 * \param records AMC records (from readDsnAmcFile)
 * \param minimumRejectedCloudFlag Cloud flag value at and above which records are rejected (default 1: only records with cloud
 * flag 0, i.e. liquid water below 20 micrometer, are retained)
 * \param rejectOnElevationFlag Boolean denoting whether records with elevation outside the prescribed range are rejected
 * \param rejectOnMtpFlag Boolean denoting whether records with non-zero MTP flag are rejected
 * \param rejectOnSurfaceFlag Boolean denoting whether records with non-zero surface meteorology flag are rejected
 * \param rejectOnWindFlag Boolean denoting whether records with non-zero wind flag are rejected
 * \return Zenith [dry, wet] path delays [m], with UTC time (seconds since J2000) as key.
 */
std::map< double, Eigen::Vector2d > getValidDsnAmcZenithDelays( const std::vector< DsnAmcRecord >& records,
                                                                const int minimumRejectedCloudFlag = 1,
                                                                const bool rejectOnElevationFlag = true,
                                                                const bool rejectOnMtpFlag = true,
                                                                const bool rejectOnSurfaceFlag = true,
                                                                const bool rejectOnWindFlag = true );

}  // namespace input_output
}  // namespace tudat

#endif  // TUDAT_READDSNAMCFILE_H
