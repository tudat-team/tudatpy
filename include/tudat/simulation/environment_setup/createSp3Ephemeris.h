/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATESP3EPHEMERIS_H
#define TUDAT_CREATESP3EPHEMERIS_H

#include <memory>
#include <string>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

namespace tudat
{

namespace input_output
{
struct Sp3FileContents;
}

namespace simulation_setup
{

//! Create tabulated ephemeris settings for one satellite from parsed SP3 file contents.
/*!
 * Position-only SP3 files are supported through the second-order velocities reconstructed by readSp3File. The
 * coordinate-system name is retained as frame metadata when frameOrientation is empty. If frameOrientation differs,
 * supported terrestrial-realization transformations or Earth-fixed/inertial rotations are applied to both position
 * and velocity. An exception is thrown when the requested transformation is unavailable. Epoch keys remain expressed
 * in the time system declared by the SP3 file.
 */
std::shared_ptr< EphemerisSettings > sp3EphemerisSettings( const std::shared_ptr< input_output::Sp3FileContents >& fileContents,
                                                           const std::string& satelliteId,
                                                           const std::string& frameOrigin = "Earth",
                                                           const std::string& frameOrientation = "" );

//! Read an SP3 file and create tabulated ephemeris settings for one satellite.
std::shared_ptr< EphemerisSettings > sp3EphemerisSettings( const std::string& fileName,
                                                           const std::string& satelliteId,
                                                           const std::string& frameOrigin = "Earth",
                                                           const std::string& frameOrientation = "",
                                                           const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATESP3EPHEMERIS_H
