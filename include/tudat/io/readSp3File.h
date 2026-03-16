/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_READSP3FILE_H
#define TUDAT_READSP3FILE_H

#include <map>
#include <memory>
#include <string>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{

namespace input_output
{

//! Parsed contents of an SP3 orbit file.
struct Sp3FileContents {
    //! Start epoch from SP3 header [s since referenceJulianDay].
    double startEpoch = TUDAT_NAN;

    //! Declared number of epochs from SP3 header line '#'.
    int declaredNumberOfEpochs = -1;

    //! Declared epoch interval from SP3 header line '##' [s].
    double declaredEpochInterval = TUDAT_NAN;

    //! Reference frame name from the SP3 header.
    std::string frameName;

    //! Producer/analysis center tag from the SP3 header.
    std::string analysisCenter;

    //! Time scale name from the SP3 header.
    std::string timeScale;

    //! Satellite state history [m, m/s], keyed by satellite id and epoch [s since referenceJulianDay].
    std::map< std::string, std::map< double, Eigen::VectorXd > > satelliteStates;
};

//! Read an SP3 file and convert states to SI units (position in meters, velocity in m/s).
std::shared_ptr< Sp3FileContents > readSp3File( const std::string& fileName,
                                                const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

//! Legacy aliases kept for backwards compatibility with older code.
using SP3cFileContents = Sp3FileContents;

inline std::shared_ptr< SP3cFileContents > readSp3cFile( const std::string& fileName,
                                                         const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 )
{
    return readSp3File( fileName, referenceJulianDay );
}

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READSP3FILE_H
