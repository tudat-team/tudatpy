/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    CALCEPH-based ephemeris implementation for direct binary SPK reading.
 */

#ifdef TUDAT_BUILD_WITH_CALCEPH

#include "tudat/astro/ephemerides/calcephEphemeris.h"

#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cctype>

namespace tudat
{
namespace ephemerides
{

// ============================================================================
// CalcephEphemeris Implementation
// ============================================================================

CalcephEphemeris::CalcephEphemeris(
    const std::string& spkFilePath,
    int targetNaifId,
    int observerNaifId,
    const std::string& referenceFrameOrigin,
    const std::string& referenceFrameOrientation )
    : Ephemeris( referenceFrameOrigin, referenceFrameOrientation ),
      ephemerisHandle_( nullptr ),
      targetNaifId_( targetNaifId ),
      observerNaifId_( observerNaifId ),
      startEpoch_( 0.0 ),
      endEpoch_( 0.0 )
{
    // Open the ephemeris file
    ephemerisHandle_ = calceph_open( spkFilePath.c_str( ) );

    if ( ephemerisHandle_ == nullptr )
    {
        throw std::runtime_error( "CalcephEphemeris: Failed to open SPK file: " + spkFilePath );
    }

    // Get time span covered by the ephemeris
    double firstJD, lastJD;
    int continuous;
    if ( calceph_gettimespan( ephemerisHandle_, &firstJD, &lastJD, &continuous ) )
    {
        // Convert Julian Date to seconds since J2000
        // J2000 epoch = JD 2451545.0
        const double J2000_JD = 2451545.0;
        const double SECONDS_PER_DAY = 86400.0;

        startEpoch_ = ( firstJD - J2000_JD ) * SECONDS_PER_DAY;
        endEpoch_ = ( lastJD - J2000_JD ) * SECONDS_PER_DAY;

        std::cout << "[CALCEPH] Loaded " << spkFilePath
                  << " (target=" << targetNaifId << ", observer=" << observerNaifId
                  << ", JD " << firstJD << " to " << lastJD << ")" << std::endl;
    }
    else
    {
        std::cerr << "[CALCEPH] Warning: Could not get time span for " << spkFilePath << std::endl;
    }

    // Prefetch data to memory for better performance
    if ( calceph_prefetch( ephemerisHandle_ ) == 0 )
    {
        std::cerr << "[CALCEPH] Warning: Prefetch failed for " << spkFilePath << std::endl;
    }
}

CalcephEphemeris::~CalcephEphemeris( )
{
    std::lock_guard<std::mutex> lock( mutex_ );
    if ( ephemerisHandle_ != nullptr )
    {
        calceph_close( ephemerisHandle_ );
        ephemerisHandle_ = nullptr;
    }
}

Eigen::Vector6d CalcephEphemeris::getCartesianState( double secondsSinceEpoch )
{
    std::lock_guard<std::mutex> lock( mutex_ );

    if ( ephemerisHandle_ == nullptr )
    {
        throw std::runtime_error( "CalcephEphemeris: Ephemeris file not loaded" );
    }

    // Convert seconds since J2000 to Julian Date
    const double J2000_JD = 2451545.0;
    const double SECONDS_PER_DAY = 86400.0;

    double julianDate = J2000_JD + secondsSinceEpoch / SECONDS_PER_DAY;

    // Split into integer and fractional parts for better precision
    // CALCEPH uses JD0 + time format where JD0 is the integer part
    double JD0 = std::floor( julianDate );
    double time = julianDate - JD0;

    // Output array: position and velocity
    double PV[6];

    // Use CALCEPH_UNIT_KM | CALCEPH_UNIT_SEC for km and km/s output
    // Use CALCEPH_USE_NAIFID to use NAIF IDs instead of CALCEPH's internal numbering
    int unit = CALCEPH_UNIT_KM | CALCEPH_UNIT_SEC | CALCEPH_USE_NAIFID;

    int result = calceph_compute_unit(
        ephemerisHandle_,
        JD0,
        time,
        targetNaifId_,
        observerNaifId_,
        unit,
        PV );

    if ( result == 0 )
    {
        throw std::runtime_error( "CalcephEphemeris: Failed to compute state at epoch " +
                                  std::to_string( secondsSinceEpoch ) + " s" );
    }

    // Convert km to m, km/s to m/s
    Eigen::Vector6d state;
    state( 0 ) = PV[0] * 1000.0;  // x
    state( 1 ) = PV[1] * 1000.0;  // y
    state( 2 ) = PV[2] * 1000.0;  // z
    state( 3 ) = PV[3] * 1000.0;  // vx
    state( 4 ) = PV[4] * 1000.0;  // vy
    state( 5 ) = PV[5] * 1000.0;  // vz

    return state;
}

std::pair<double, double> CalcephEphemeris::getTimeBounds( ) const
{
    return std::make_pair( startEpoch_, endEpoch_ );
}

// ============================================================================
// CalcephEphemerisManager Implementation
// ============================================================================

CalcephEphemerisManager& CalcephEphemerisManager::getInstance( )
{
    static CalcephEphemerisManager instance;
    return instance;
}

bool CalcephEphemerisManager::loadSpkFile(
    const std::string& spkFilePath,
    const std::string& targetName,
    const std::string& observerName,
    const std::string& frame )
{
    int targetId = bodyNameToNaifId( targetName );
    int observerId = bodyNameToNaifId( observerName );

    if ( targetId == -1 || observerId == -1 )
    {
        std::cerr << "[CALCEPH] Unknown body name: " << targetName << " or " << observerName << std::endl;
        return false;
    }

    return loadSpkFileByNaifId( spkFilePath, targetId, observerId, frame );
}

bool CalcephEphemerisManager::loadSpkFileByNaifId(
    const std::string& spkFilePath,
    int targetNaifId,
    int observerNaifId,
    const std::string& frame )
{
    std::lock_guard<std::mutex> lock( mutex_ );

    try
    {
        auto ephemeris = std::make_shared<CalcephEphemeris>(
            spkFilePath, targetNaifId, observerNaifId, naifIdToBodyName( observerNaifId ), frame );

        if ( !ephemeris->isLoaded( ) )
        {
            return false;
        }

        std::string key = makeKey(
            naifIdToBodyName( targetNaifId ),
            naifIdToBodyName( observerNaifId ),
            frame );

        ephemerides_[key] = ephemeris;
        return true;
    }
    catch ( const std::exception& e )
    {
        std::cerr << "[CALCEPH] Error loading SPK file: " << e.what( ) << std::endl;
        return false;
    }
}

bool CalcephEphemerisManager::isAvailable(
    const std::string& targetName,
    const std::string& observerName,
    const std::string& frame ) const
{
    std::lock_guard<std::mutex> lock( mutex_ );
    std::string key = makeKey( targetName, observerName, frame );
    return ephemerides_.find( key ) != ephemerides_.end( );
}

Eigen::Vector6d CalcephEphemerisManager::getState(
    const std::string& targetName,
    const std::string& observerName,
    const std::string& frame,
    double secondsSinceJ2000 ) const
{
    std::lock_guard<std::mutex> lock( mutex_ );
    std::string key = makeKey( targetName, observerName, frame );

    auto it = ephemerides_.find( key );
    if ( it == ephemerides_.end( ) )
    {
        throw std::runtime_error( "CalcephEphemerisManager: Ephemeris not loaded for " + key );
    }

    return it->second->getCartesianState( secondsSinceJ2000 );
}

std::pair<double, double> CalcephEphemerisManager::getTimeBounds(
    const std::string& targetName,
    const std::string& observerName,
    const std::string& frame ) const
{
    std::lock_guard<std::mutex> lock( mutex_ );
    std::string key = makeKey( targetName, observerName, frame );

    auto it = ephemerides_.find( key );
    if ( it == ephemerides_.end( ) )
    {
        return std::make_pair( 0.0, 0.0 );
    }

    return it->second->getTimeBounds( );
}

void CalcephEphemerisManager::clearAll( )
{
    std::lock_guard<std::mutex> lock( mutex_ );
    ephemerides_.clear( );
    std::cout << "[CALCEPH] Cleared all loaded ephemeris files" << std::endl;
}

std::vector<std::string> CalcephEphemerisManager::listLoaded( ) const
{
    std::lock_guard<std::mutex> lock( mutex_ );
    std::vector<std::string> keys;
    for ( const auto& pair : ephemerides_ )
    {
        keys.push_back( pair.first );
    }
    return keys;
}

std::string CalcephEphemerisManager::normalizeBodyName( const std::string& name )
{
    // Normalize body names to simple planet names
    // This ensures "Earth", "Earth-Moon Barycenter", "EMB" all map to "earth"
    std::string lowerName = name;
    std::transform( lowerName.begin( ), lowerName.end( ), lowerName.begin( ), ::tolower );

    // Map barycenter names to simple planet names
    if ( lowerName == "mercury barycenter" || lowerName == "mercury_barycenter" ) return "mercury";
    if ( lowerName == "venus barycenter" || lowerName == "venus_barycenter" ) return "venus";
    if ( lowerName == "earth-moon barycenter" || lowerName == "earth_moon_barycenter" ||
         lowerName == "emb" || lowerName == "earth barycenter" ) return "earth";
    if ( lowerName == "mars barycenter" || lowerName == "mars_barycenter" ) return "mars";
    if ( lowerName == "jupiter barycenter" || lowerName == "jupiter_barycenter" ) return "jupiter";
    if ( lowerName == "saturn barycenter" || lowerName == "saturn_barycenter" ) return "saturn";
    if ( lowerName == "uranus barycenter" || lowerName == "uranus_barycenter" ) return "uranus";
    if ( lowerName == "neptune barycenter" || lowerName == "neptune_barycenter" ) return "neptune";
    if ( lowerName == "pluto barycenter" || lowerName == "pluto_barycenter" ) return "pluto";

    // Map "xxx center" names to simple names
    if ( lowerName == "mercury center" ) return "mercury";
    if ( lowerName == "venus center" ) return "venus";
    if ( lowerName == "earth center" ) return "earth";
    if ( lowerName == "mars center" ) return "mars";
    if ( lowerName == "jupiter center" ) return "jupiter";
    if ( lowerName == "saturn center" ) return "saturn";
    if ( lowerName == "uranus center" ) return "uranus";
    if ( lowerName == "neptune center" ) return "neptune";
    if ( lowerName == "pluto center" ) return "pluto";

    // Map SSB variations
    if ( lowerName == "solar system barycenter" ) return "ssb";

    // Return lowercase version for other names
    return lowerName;
}

std::string CalcephEphemerisManager::makeKey(
    const std::string& target,
    const std::string& observer,
    const std::string& frame ) const
{
    // Normalize body names so "Earth", "EMB", "Earth-Moon Barycenter" all create the same key
    std::string normTarget = normalizeBodyName( target );
    std::string normObserver = normalizeBodyName( observer );
    std::string normFrame = frame;
    std::transform( normFrame.begin( ), normFrame.end( ), normFrame.begin( ), ::tolower );

    return normTarget + "_" + normObserver + "_" + normFrame;
}

int CalcephEphemerisManager::bodyNameToNaifId( const std::string& name )
{
    // Convert to lowercase for comparison
    std::string lowerName = name;
    std::transform( lowerName.begin( ), lowerName.end( ), lowerName.begin( ), ::tolower );

    // Solar System Barycenter and Sun
    if ( lowerName == "ssb" || lowerName == "solar system barycenter" ) return 0;
    if ( lowerName == "sun" ) return 10;

    // Planetary barycenters (explicit)
    if ( lowerName == "mercury barycenter" || lowerName == "mercury_barycenter" ) return 1;
    if ( lowerName == "venus barycenter" || lowerName == "venus_barycenter" ) return 2;
    if ( lowerName == "earth-moon barycenter" || lowerName == "earth_moon_barycenter" ||
         lowerName == "emb" || lowerName == "earth barycenter" ) return 3;
    if ( lowerName == "mars barycenter" || lowerName == "mars_barycenter" ) return 4;
    if ( lowerName == "jupiter barycenter" || lowerName == "jupiter_barycenter" ) return 5;
    if ( lowerName == "saturn barycenter" || lowerName == "saturn_barycenter" ) return 6;
    if ( lowerName == "uranus barycenter" || lowerName == "uranus_barycenter" ) return 7;
    if ( lowerName == "neptune barycenter" || lowerName == "neptune_barycenter" ) return 8;
    if ( lowerName == "pluto barycenter" || lowerName == "pluto_barycenter" ) return 9;

    // Planet names -> Barycenters
    // DE kernels (de430, de432, de440, etc.) contain planetary barycenter data,
    // not planet centers. For interplanetary trajectory design, the barycenter
    // is the appropriate approximation since it represents the system's center of mass.
    // The difference between barycenter and planet center is:
    //   - Mercury/Venus: ~0 (no significant moons)
    //   - Earth: ~4670 km from Earth center (due to Moon)
    //   - Mars: ~few km (small moons)
    //   - Jupiter/Saturn: varies by thousands of km (large moon systems)
    // For mission design, using barycenters is standard practice.
    if ( lowerName == "mercury" ) return 1;   // Mercury Barycenter (no moons)
    if ( lowerName == "venus" ) return 2;     // Venus Barycenter (no moons)
    if ( lowerName == "earth" ) return 3;     // Earth-Moon Barycenter
    if ( lowerName == "mars" ) return 4;      // Mars Barycenter
    if ( lowerName == "jupiter" ) return 5;   // Jupiter Barycenter
    if ( lowerName == "saturn" ) return 6;    // Saturn Barycenter
    if ( lowerName == "uranus" ) return 7;    // Uranus Barycenter
    if ( lowerName == "neptune" ) return 8;   // Neptune Barycenter
    if ( lowerName == "pluto" ) return 9;     // Pluto Barycenter

    // Moon - requires Earth-Moon relative data from separate kernel segment
    if ( lowerName == "moon" || lowerName == "luna" ) return 301;

    // Planet centers (explicit) - for kernels that have planet center data
    // These use the 1xx/2xx/etc. NAIF IDs for the planet body itself
    if ( lowerName == "mercury center" ) return 199;
    if ( lowerName == "venus center" ) return 299;
    if ( lowerName == "earth center" ) return 399;
    if ( lowerName == "mars center" ) return 499;
    if ( lowerName == "jupiter center" ) return 599;
    if ( lowerName == "saturn center" ) return 699;
    if ( lowerName == "uranus center" ) return 799;
    if ( lowerName == "neptune center" ) return 899;
    if ( lowerName == "pluto center" ) return 999;

    // Major moons
    if ( lowerName == "phobos" ) return 401;
    if ( lowerName == "deimos" ) return 402;
    if ( lowerName == "io" ) return 501;
    if ( lowerName == "europa" ) return 502;
    if ( lowerName == "ganymede" ) return 503;
    if ( lowerName == "callisto" ) return 504;
    if ( lowerName == "titan" ) return 606;
    if ( lowerName == "triton" ) return 801;
    if ( lowerName == "charon" ) return 901;

    return -1;  // Unknown
}

std::string CalcephEphemerisManager::naifIdToBodyName( int naifId )
{
    switch ( naifId )
    {
        case 0: return "SSB";
        case 1: return "Mercury Barycenter";
        case 2: return "Venus Barycenter";
        case 3: return "Earth-Moon Barycenter";
        case 4: return "Mars Barycenter";
        case 5: return "Jupiter Barycenter";
        case 6: return "Saturn Barycenter";
        case 7: return "Uranus Barycenter";
        case 8: return "Neptune Barycenter";
        case 9: return "Pluto Barycenter";
        case 10: return "Sun";
        case 199: return "Mercury";
        case 299: return "Venus";
        case 301: return "Moon";
        case 399: return "Earth";
        case 401: return "Phobos";
        case 402: return "Deimos";
        case 499: return "Mars";
        case 501: return "Io";
        case 502: return "Europa";
        case 503: return "Ganymede";
        case 504: return "Callisto";
        case 599: return "Jupiter";
        case 606: return "Titan";
        case 699: return "Saturn";
        case 799: return "Uranus";
        case 801: return "Triton";
        case 899: return "Neptune";
        case 901: return "Charon";
        case 999: return "Pluto";
        default: return "Body" + std::to_string( naifId );
    }
}

}  // namespace ephemerides
}  // namespace tudat

#endif  // TUDAT_BUILD_WITH_CALCEPH
