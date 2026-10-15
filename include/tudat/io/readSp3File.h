/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_READ_SP3_FILE_H
#define TUDAT_READ_SP3_FILE_H

#include <cstddef>
#include <istream>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/basics/timeType.h"

namespace tudat
{
namespace input_output
{

//! One SP3 satellite sample, with SI values and an epoch in the file's declared time coordinate.
struct Sp3Record {

    //! Original epoch and source line of the position record; no time conversion is applied.
    Time epoch;

    std::size_t sourceLine = 0;

    //! Position in metres and velocity in metres/second; zero-valued missing triples remain unavailable.
    Eigen::Vector3d position = Eigen::Vector3d::Constant( TUDAT_NAN );

    Eigen::Vector3d velocity = Eigen::Vector3d::Constant( TUDAT_NAN );

    bool positionAvailable = false;

    bool velocityAvailable = false;

    //! Clock correction in seconds and rate in seconds/second; missing sentinels become NaN.
    double clockCorrection = TUDAT_NAN;

    double clockRate = TUDAT_NAN;

    //! Published event and prediction flags; false also covers blank/unspecified flags.
    bool clockEvent = false;

    bool clockPredicted = false;

    bool maneuver = false;

    bool orbitPredicted = false;
};

//! Raw SP3c/d ephemeris data; the frame label and comments are preserved without assigning an inertial frame.
struct Sp3Data {

    //! Source and format/header metadata, including FIT versus prediction orbit type.
    std::string source;

    char version = '\0';

    bool hasVelocityRecords = false;

    int declaredEpochCount = 0;

    double epochInterval = TUDAT_NAN;

    std::string dataUsed;

    std::string coordinateSystem;

    std::string orbitType;

    std::string agency;

    std::string timeSystem;

    //! Provider comments and declared satellite identifiers in their original order.
    std::vector< std::string > comments;

    std::vector< std::string > satelliteIds;

    //! Ordered samples per satellite, using the printed satellite IDs (e.g. L51 for LAGEOS-1).
    std::map< std::string, std::vector< Sp3Record > > records;
};

//! Read SP3c/d records, converting km and dm/s to SI without converting epochs or interpolating an orbit.
//! Missing/predicted records remain distinguishable; optional EP/EV covariances are not consumed.
Sp3Data readSp3File( std::istream& stream, const std::string& source = "<stream>" );

//! Open an uncompressed SP3 file and retain its path as provenance.
Sp3Data readSp3File( const std::string& fileName );

}  // namespace input_output
}  // namespace tudat

#endif
