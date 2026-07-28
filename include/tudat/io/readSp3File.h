/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
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
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{

namespace input_output
{

//! Parsed contents of an SP3 orbit file.
struct Sp3FileContents {
    //! SP3 format version ('a', 'b', 'c', or 'd').
    char formatVersion = '\0';

    //! True when the file declares measured/published velocity records ('V' mode).
    bool hasVelocityRecords = false;

    //! True when velocities were reconstructed from positions in a position-only ('P' mode) file.
    bool velocitiesWereDerived = false;

    //! Start epoch from SP3 header [s since referenceJulianDay].
    double startEpoch = TUDAT_NAN;

    //! Julian day relative to which all epochs in satelliteStates are expressed.
    double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000;

    //! Declared number of epochs from SP3 header line '#'.
    int declaredNumberOfEpochs = -1;

    //! Declared number of satellites from the first SP3 header line '+'.
    int declaredNumberOfSatellites = -1;

    //! Declared epoch interval from SP3 header line '##' [s].
    double declaredEpochInterval = TUDAT_NAN;

    //! Reference frame name from the SP3 header.
    std::string frameName;

    //! Producer/analysis center tag from the SP3 header.
    std::string analysisCenter;

    //! Time-system code from the SP3 header. The legacy "ccc" default is normalized to "GPS".
    std::string timeScale;

    //! Satellite state history [m, m/s], keyed by satellite id and epoch [s since referenceJulianDay in timeScale].
    /*!
     * Position-only files receive second-order finite-difference velocities. Three epochs are therefore required
     * for a position-only file. Missing position values propagate to the affected reconstructed velocities.
     */
    std::map< std::string, std::map< double, Eigen::Vector6d > > satelliteStates;

    //! Retrieve the state history for one satellite.
    /*!
     * \param satelliteId Identifier of the satellite whose states are to be returned.
     * \return Satellite state history in the frame declared by the SP3 file.
     */
    std::map< double, Eigen::Vector6d > getSatelliteStateHistory( const std::string& satelliteId ) const;
};

//! Read an SP3-a, SP3-b, SP3-c, or SP3-d file and convert states to SI units.
/*!
 * The official GPS, GLO, GAL, BDT, TAI, UTC, IRN, and QZS time-system tags are accepted and retained as metadata.
 * Epoch values are not converted between time systems. Coordinate-system tags are also retained verbatim; no
 * reference-frame transformation is performed by the reader. EP/EV covariance records and clock fields are currently ignored.
 */
std::shared_ptr< Sp3FileContents > readSp3File( const std::string& fileName,
                                                const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READSP3FILE_H
