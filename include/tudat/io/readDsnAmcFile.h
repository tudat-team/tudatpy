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
    //! Azimuth and elevation of radiometer [deg]
    double azimuth;
    double elevation;
    //! Zenith wet and dry path delays [m]
    double zenithWetDelay;
    double zenithDryDelay;
    //! Cloud liquid [micrometer]
    double cloudLiquid;
    //! Identifier of retrieval algorithm
    std::string algorithmId;
    //! Quality flags (0 is nominal)
    int cloudFlag;
    int wvrFlag;
    int mtpFlag;
    int surfaceFlag;
    int windFlag;
    int pathDelayFlag;
    //! Wind speed [mph] and direction [deg]
    double windSpeed;
    double windDirection;
};

//! Function to read all records of a DSN AMC path delay file (without any quality screening)
/*!
 * Function to read all records of a DSN AMC path delay file (without any quality screening). The time is converted from UTC
 * seconds since 1970 to UTC seconds since J2000, and path delays from cm to m. Other quantities are in the units of the file.
 * \param fileName Name of AMC file
 * \return List of records in file
 */
std::vector< DsnAmcRecord > readDsnAmcFile( const std::string& fileName );

//! Function to retrieve valid zenith path delays from DSN AMC records
/*!
 * Function to retrieve valid zenith path delays from DSN AMC records, retaining only records for which all quality flags are 0, and
 * the path delays are valid (non-negative, and not equal to the 999 cm placeholder for a missing retrieval).
 * \param records AMC records (from readDsnAmcFile)
 * \return Zenith [dry, wet] path delays [m], with UTC time (seconds since J2000) as key.
 */
std::map< double, Eigen::Vector2d > getValidDsnAmcZenithDelays( const std::vector< DsnAmcRecord >& records );

}  // namespace input_output
}  // namespace tudat

#endif  // TUDAT_READDSNAMCFILE_H
