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
#include <utility>
#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace input_output
{

//! Single record of a DSN Advanced Media Calibration (AMC) path delay file (see DSN document DORS_002_020501)
struct DsnAmcRecord {
    //! UTC time, in seconds since J2000 (converted from UTC seconds since 1970 in the file)
    double utcTime;

    //! Azimuth of the Advanced Water Vapor Radiometer (AWVR) [deg]
    double azimuth;

    //! Elevation of the AWVR [deg]
    double elevation;

    //! Zenith wet path delay [m] (converted from cm in the file; 999 cm denotes a missing retrieval)
    double zenithWetDelay;

    //! Zenith dry (static) path delay [m] (converted from cm in the file; 999 cm denotes a missing retrieval)
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

    //! Wind speed [mph]
    double windSpeed;

    //! Wind direction [deg]
    double windDirection;
};

//! Function to read the records of a DSN AMC path delay file
/*!
 * Function to read the records of a DSN AMC path delay file. The time is converted from UTC seconds since 1970 to UTC seconds
 * since J2000, and path delays from cm to m. Other quantities are in the units of the file.
 * \param fileName Name of AMC file
 * \param filterInvalidRecords Boolean denoting whether only valid records are retained: those for which all quality flags are 0,
 * and the path delays are valid (non-negative, and not equal to the 999 cm placeholder for a missing retrieval)
 * \return Pair with list of records in file, and the zenith [dry, wet] path delays [m] of these records, with UTC time (seconds
 * since J2000) as key
 */
std::pair< std::vector< DsnAmcRecord >, std::map< double, Eigen::Vector2d > > readDsnAmcFile( const std::string& fileName,
                                                                                              const bool filterInvalidRecords = true );

}  // namespace input_output
}  // namespace tudat

#endif  // TUDAT_READDSNAMCFILE_H
