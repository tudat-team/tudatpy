/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved.
 *    This file is part of Tudat and is distributed under the Modified BSD license.
 */

#ifndef TUDAT_READDSNAMCFILE_H
#define TUDAT_READDSNAMCFILE_H

#include <memory>
#include <string>
#include <vector>

namespace tudat
{
namespace input_output
{

//! A DSN Advanced Media Calibration path-delay record, without quality screening.
struct DsnAmcRecord {
    double utcTime;           //!< UTC seconds since J2000 (2000-01-01 at noon).
    double azimuth;           //!< AWVR azimuth [deg], as recorded.
    double elevation;         //!< AWVR elevation [deg], as recorded.
    double zenithWetDelay;    //!< Zenith wet path length [m].
    double zenithDryDelay;    //!< Zenith dry path length [m].
    double cloudLiquid;       //!< Cloud liquid [micrometres].
    std::string algorithmId;  //!< Retrieval identifier, including leading zeroes.
    int cloudFlag;
    int wvrFlag;
    int mtpFlag;
    int surfaceFlag;
    int windFlag;
    int pathDelayFlag;
    double windSpeed;      //!< Wind speed [mph], as recorded.
    double windDirection;  //!< Wind direction [deg], as recorded.
};

//! Records from one AMC file. Unit numbers in filenames do not identify DSN stations.
class DsnAmcData
{
public:
    DsnAmcData( const std::vector< DsnAmcRecord >& records, const std::string& fileName = "" ): records_( records ), fileName_( fileName )
    {}

    const std::vector< DsnAmcRecord >& getRecords( ) const
    {
        return records_;
    }

    const std::string& getFileName( ) const
    {
        return fileName_;
    }

private:
    std::vector< DsnAmcRecord > records_;
    std::string fileName_;
};

/*!
 * Read a DSN AMC path-delay file (P records, with # header lines).
 * The 1970 UT epochs are converted to UTC seconds since J2000, and centimetres
 * of zenith wet/dry path length to metres. Other fields retain their file units.
 * No quality filtering, interpolation, station association or gap filling is done.
 * Malformed records, non-finite numbers, duplicate/decreasing epochs and files
 * without records throw std::runtime_error. Unknown integer flags and algorithm
 * identifiers, including the Juno identifier 010, are preserved.
 */
std::shared_ptr< DsnAmcData > readDsnAmcFile( const std::string& fileName );

}  // namespace input_output
}  // namespace tudat

#endif  // TUDAT_READDSNAMCFILE_H
