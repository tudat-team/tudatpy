/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    CALCEPH-based ephemeris for direct binary SPK reading in WASM.
 *    CALCEPH is a pure C library that can read SPK files without CSPICE's
 *    furnsh_c() which crashes in WASM due to f2c FORTRAN I/O issues.
 */

#ifndef TUDAT_CALCEPH_EPHEMERIS_H
#define TUDAT_CALCEPH_EPHEMERIS_H

#ifdef TUDAT_BUILD_WITH_CALCEPH

#include <string>
#include <map>
#include <memory>
#include <mutex>

#include <Eigen/Core>

#include "tudat/astro/ephemerides/ephemeris.h"

// CALCEPH C API
extern "C" {
#include "calceph.h"
}

namespace tudat
{
namespace ephemerides
{

//! Class to provide ephemeris data from binary SPK files using CALCEPH.
/*!
 * This class uses the CALCEPH library (from IMCCE) to read binary SPK files
 * directly, without going through CSPICE's furnsh_c() function. This is
 * essential for WASM builds where CSPICE's f2c-generated FORTRAN I/O code
 * crashes due to incompatibilities with WebAssembly.
 *
 * CALCEPH supports SPK segment types: 1, 2, 3, 5, 8, 9, 12, 13, 14, 17, 18, 20, 21
 * (plus INPOP types 102, 103, 120).
 *
 * Usage:
 *   1. Create a CalcephEphemeris for each SPK file
 *   2. Query states using getCartesianState()
 *   3. The file handle is automatically closed on destruction
 */
class CalcephEphemeris : public Ephemeris
{
public:
    //! Constructor - opens an SPK file for reading.
    /*!
     * \param spkFilePath Path to the binary SPK file (or INPOP file)
     * \param targetNaifId NAIF ID of the target body (e.g., 399 for Earth)
     * \param observerNaifId NAIF ID of the observer/center body (e.g., 10 for Sun)
     * \param referenceFrameOrigin Name of the reference frame origin (for Tudat compatibility)
     * \param referenceFrameOrientation Name of the reference frame (for Tudat compatibility)
     */
    CalcephEphemeris(
        const std::string& spkFilePath,
        int targetNaifId,
        int observerNaifId,
        const std::string& referenceFrameOrigin = "SSB",
        const std::string& referenceFrameOrientation = "J2000" );

    //! Destructor - closes the ephemeris file handle.
    ~CalcephEphemeris( );

    //! Get Cartesian state at epoch.
    /*!
     * Returns the Cartesian state [x, y, z, vx, vy, vz] in meters and m/s.
     * \param secondsSinceEpoch Seconds since J2000 epoch (TDB)
     * \return State vector (6 elements: position in m, velocity in m/s)
     */
    Eigen::Vector6d getCartesianState( double secondsSinceEpoch ) override;

    //! Get the time bounds covered by this ephemeris.
    /*!
     * \return Pair of (startEpoch, endEpoch) in seconds since J2000
     */
    std::pair<double, double> getTimeBounds( ) const;

    //! Check if the ephemeris file was loaded successfully.
    bool isLoaded( ) const { return ephemerisHandle_ != nullptr; }

    //! Get the target NAIF ID.
    int getTargetNaifId( ) const { return targetNaifId_; }

    //! Get the observer NAIF ID.
    int getObserverNaifId( ) const { return observerNaifId_; }

private:
    //! CALCEPH ephemeris handle
    t_calcephbin* ephemerisHandle_;

    //! Target body NAIF ID
    int targetNaifId_;

    //! Observer body NAIF ID
    int observerNaifId_;

    //! Time bounds (seconds since J2000)
    double startEpoch_;
    double endEpoch_;

    //! Mutex for thread safety (CALCEPH file access)
    mutable std::mutex mutex_;
};

//! Global manager for CALCEPH-based ephemeris files.
/*!
 * This singleton class manages multiple SPK files loaded via CALCEPH,
 * providing a unified interface for querying body states.
 */
class CalcephEphemerisManager
{
public:
    //! Get the singleton instance
    static CalcephEphemerisManager& getInstance( );

    //! Load an SPK file for the specified target/observer pair.
    /*!
     * \param spkFilePath Path to the binary SPK file
     * \param targetName Name of the target body (converted to NAIF ID)
     * \param observerName Name of the observer body (converted to NAIF ID)
     * \param frame Reference frame name
     * \return True if loaded successfully
     */
    bool loadSpkFile(
        const std::string& spkFilePath,
        const std::string& targetName,
        const std::string& observerName,
        const std::string& frame = "J2000" );

    //! Load an SPK file using NAIF IDs directly.
    /*!
     * \param spkFilePath Path to the binary SPK file
     * \param targetNaifId NAIF ID of the target body
     * \param observerNaifId NAIF ID of the observer body
     * \param frame Reference frame name
     * \return True if loaded successfully
     */
    bool loadSpkFileByNaifId(
        const std::string& spkFilePath,
        int targetNaifId,
        int observerNaifId,
        const std::string& frame = "J2000" );

    //! Check if ephemeris is available for a target/observer pair.
    bool isAvailable(
        const std::string& targetName,
        const std::string& observerName,
        const std::string& frame = "J2000" ) const;

    //! Get state of target relative to observer.
    /*!
     * \param targetName Name of the target body
     * \param observerName Name of the observer body
     * \param frame Reference frame name
     * \param secondsSinceJ2000 Epoch in seconds since J2000
     * \return Cartesian state [x, y, z, vx, vy, vz] in m and m/s
     */
    Eigen::Vector6d getState(
        const std::string& targetName,
        const std::string& observerName,
        const std::string& frame,
        double secondsSinceJ2000 ) const;

    //! Get time bounds for a target/observer pair.
    std::pair<double, double> getTimeBounds(
        const std::string& targetName,
        const std::string& observerName,
        const std::string& frame = "J2000" ) const;

    //! Clear all loaded ephemeris files.
    void clearAll( );

    //! List all loaded ephemeris keys.
    std::vector<std::string> listLoaded( ) const;

    //! Convert body name to NAIF ID.
    static int bodyNameToNaifId( const std::string& name );

    //! Convert NAIF ID to body name.
    static std::string naifIdToBodyName( int naifId );

    //! Normalize body name for key lookup.
    /*!
     * Normalizes different variations of body names to a standard form.
     * E.g., "Earth", "EMB", "Earth-Moon Barycenter" all -> "earth"
     */
    static std::string normalizeBodyName( const std::string& name );

private:
    CalcephEphemerisManager( ) = default;
    ~CalcephEphemerisManager( ) = default;
    CalcephEphemerisManager( const CalcephEphemerisManager& ) = delete;
    CalcephEphemerisManager& operator=( const CalcephEphemerisManager& ) = delete;

    //! Create key for ephemeris lookup
    std::string makeKey(
        const std::string& target,
        const std::string& observer,
        const std::string& frame ) const;

    //! Loaded ephemeris objects
    std::map<std::string, std::shared_ptr<CalcephEphemeris>> ephemerides_;

    //! Mutex for thread safety
    mutable std::mutex mutex_;
};

}  // namespace ephemerides
}  // namespace tudat

#endif  // TUDAT_BUILD_WITH_CALCEPH

#endif  // TUDAT_CALCEPH_EPHEMERIS_H
