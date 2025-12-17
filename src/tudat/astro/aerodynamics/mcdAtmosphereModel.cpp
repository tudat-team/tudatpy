/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "tudat/astro/aerodynamics/mcdAtmosphereModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/utilities.h"
#include "tudat/io/basicInputOutput.h"
#include "mcd.h"
#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace aerodynamics
{

// Constructor
McdAtmosphereModel::McdAtmosphereModel( const std::string& mcdDataPath,
                                        const int dustScenario,
                                        const int perturbationKey,
                                        const double perturbationSeed,
                                        const double gravityWaveLength,
                                        const int highResolutionMode ):
    mcdDataPath_( mcdDataPath ), dustScenario_( dustScenario ), perturbationKey_( perturbationKey ), perturbationSeed_( perturbationSeed ),
    gravityWaveLength_( gravityWaveLength ), highResolutionMode_( highResolutionMode ), currentDensity_( 0.0 ), currentPressure_( 0.0 ),
    currentTemperature_( 0.0 ), currentZonalWind_( 0.0 ), currentMeridionalWind_( 0.0 ), currentSeedOut_( perturbationSeed ),
    cachedAltitude_( TUDAT_NAN ), cachedLongitude_( TUDAT_NAN ), cachedLatitude_( TUDAT_NAN ), cachedTime_( TUDAT_NAN ),
    lastComputedWithExtras_( false ), currentRadialDistance_( 0.0 ), currentAltitudeAboveAreoid_( 0.0 ),
    currentAltitudeAboveSurface_( 0.0 ), currentOrographicHeight_( 0.0 ), currentGcmOrography_( 0.0 ), currentSlopeInclination_( 0.0 ),
    currentSlopeOrientation_( 0.0 ), currentMarsAU_( 0.0 ), currentLs_( 0.0 ), currentLTST_( 0.0 ), currentLocalMeanTime_( 0.0 ),
    currentUniversalSolarTime_( 0.0 ), currentSolarZenithAngle_( 0.0 )
{
    // Resize vectors to hold MCD outputs
    currentMeanVariables_.resize( 5, 0.0 );
    currentExtraVariables_.resize( 100, 0.0 );

    // Set default MCD data path if not provided
    if( mcdDataPath_.empty( ) )
    {
        // Use the compile-time defined path to the MCD data directory
        // This is set by CMake and points to the source tree location
#ifdef MCD_DATA_PATH
        mcdDataPath_ = MCD_DATA_PATH;
#else
        // Fallback: this should not normally happen if CMake is configured correctly
        throw std::runtime_error(
                "Error in MCD atmosphere model: MCD_DATA_PATH not defined at compile time. "
                "Please ensure MCD is properly configured in CMake." );
#endif
    }

    // Ensure path ends with '/'
    if( !mcdDataPath_.empty( ) && mcdDataPath_.back( ) != '/' )
    {
        mcdDataPath_ += '/';
    }

    // Validate input parameters (matching validation in McdAtmosphereSettings)
    if( ( dustScenario_ < 1 || dustScenario_ > 8 ) && ( dustScenario_ < 24 || dustScenario_ > 35 ) )
    {
        throw std::runtime_error( "McdAtmosphereModel: Invalid dustScenario " + std::to_string( dustScenario_ ) +
                                  ". Must be 1-8 or 24-35." );
    }

    if( perturbationKey_ < 0 || perturbationKey_ > 5 )
    {
        throw std::runtime_error( "McdAtmosphereModel: Invalid perturbationKey " + std::to_string( perturbationKey_ ) + ". Must be 0-5." );
    }

    // Validate perturbationSeed for perturbationKey=5
    if( perturbationKey_ == 5 )
    {
        if( perturbationSeed_ < -4.0 || perturbationSeed_ > 4.0 )
        {
            throw std::runtime_error(
                    "McdAtmosphereModel: For perturbationKey=5, perturbationSeed must be in [-4, 4]. "
                    "Got: " +
                    std::to_string( perturbationSeed_ ) );
        }
    }

    if( highResolutionMode_ != 0 && highResolutionMode_ != 1 )
    {
        throw std::runtime_error( "McdAtmosphereModel: Invalid highResolutionMode " + std::to_string( highResolutionMode_ ) +
                                  ". Must be 0 or 1." );
    }

    // Validate gravity wave length if perturbations are used
    if( ( perturbationKey_ == 3 || perturbationKey_ == 4 ) && gravityWaveLength_ < 0.0 )
    {
        throw std::runtime_error(
                "McdAtmosphereModel: gravityWaveLength must be >= 0.0 when using gravity wave perturbations. "
                "Got: " +
                std::to_string( gravityWaveLength_ ) );
    }
}

// Get density
double McdAtmosphereModel::getDensity( const double altitude, const double longitude, const double latitude, const double time )
{
    computePropertiesIfNeeded( altitude, longitude, latitude, time, false );
    return currentDensity_;
}

// Get pressure
double McdAtmosphereModel::getPressure( const double altitude, const double longitude, const double latitude, const double time )
{
    computePropertiesIfNeeded( altitude, longitude, latitude, time, false );
    return currentPressure_;
}

// Get temperature
double McdAtmosphereModel::getTemperature( const double altitude, const double longitude, const double latitude, const double time )
{
    computePropertiesIfNeeded( altitude, longitude, latitude, time, false );
    return currentTemperature_;
}

// Get speed of sound
double McdAtmosphereModel::getSpeedOfSound( const double altitude, const double longitude, const double latitude, const double time )
{
    // Speed of sound needs extvar(60) and extvar(61) = indices 59, 60 in C++
    computePropertiesIfNeeded( altitude, longitude, latitude, time, true );

    // Get gamma and R from extra variables
    double gamma = currentExtraVariables_[ 59 ];  // extvar(60): gamma
    double R = currentExtraVariables_[ 60 ];      // extvar(61): R (J/kg/K)

    if( gamma > 0.0 && R > 0.0 )
    {
        return std::sqrt( gamma * R * currentTemperature_ );
    }
    else
    {
        // Fallback to Mars defaults
        const double defaultGamma = 1.3;
        const double defaultR = 192.0;  // J/kg/K for Mars CO2 atmosphere
        return std::sqrt( defaultGamma * defaultR * currentTemperature_ );
    }
}

bool McdAtmosphereModel::isCached( const double altitude, const double longitude, const double latitude, const double time ) const
{
    // MCD optimization guidelines indicate computational cost hierarchy:
    // 1. Time changes (most expensive) - may trigger dataset reloading
    // 2. Horizontal position (moderate) - requires horizontal interpolation
    // 3. Vertical position (least expensive) - only vertical interpolation

    // Time tolerance: Based on MCD's temporal resolution and diurnal variations
    // - MCD uses 12 time steps per sol (every 2 hours)
    // - Martian sol = 88775.245 seconds
    // - One time step = 88775.245/12 ≈ 7398 seconds
    // - Use 1/100th of time step for smooth temporal evolution
    // - Tolerance: ~74 seconds (captures sub-hourly atmospheric changes)
    const double timeTolerance = 74.0;  // seconds

    // Horizontal position tolerance: Based on atmospheric spatial scales
    // - MCD GCM resolution: ~5.6° × 5.6° (low res) or MOLA 32 pixels/degree (high res)
    // - Mars radius ≈ 3396 km, circumference ≈ 21340 km
    // - 1° ≈ 59.3 km at equator
    // - Atmospheric features: ~100-500 km horizontal scale
    // - Use 0.001 rad ≈ 0.057° ≈ 3.4 km (captures mesoscale features)
    const double horizontalTolerance = 0.001;  // radians (~3.4 km at equator)

    // Vertical tolerance: Based on atmospheric scale heights and MCD vertical resolution
    // - Mars atmospheric scale height: ~10-11 km
    // - Density changes by factor e (2.718) per scale height
    // - MCD vertical resolution: varies from ~1 km (low altitudes) to ~5 km (high altitudes)
    // - Use 100 m tolerance: captures vertical structure without excessive recomputation
    // - At H=10 km: Δρ/ρ ≈ (Δz/H) ≈ 100/10000 = 1% density change
    const double verticalTolerance = 100.0;  // meters

    // Check if state is within tolerances
    if( std::isnan( cachedAltitude_ ) || std::isnan( cachedLongitude_ ) || std::isnan( cachedLatitude_ ) || std::isnan( cachedTime_ ) )
    {
        return false;
    }

    // Check each parameter against its physically-motivated tolerance
    bool altitudeMatch = std::abs( altitude - cachedAltitude_ ) < verticalTolerance;
    bool longitudeMatch = std::abs( longitude - cachedLongitude_ ) < horizontalTolerance;
    bool latitudeMatch = std::abs( latitude - cachedLatitude_ ) < horizontalTolerance;
    bool timeMatch = std::abs( time - cachedTime_ ) < timeTolerance;

    return altitudeMatch && longitudeMatch && latitudeMatch && timeMatch;
}

void McdAtmosphereModel::computePropertiesIfNeeded( const double altitude,
                                                    const double longitude,
                                                    const double latitude,
                                                    const double time,
                                                    const bool needExtraVariables )
{
    // Check if we can use cached values
    bool stateChanged = !isCached( altitude, longitude, latitude, time );

    // Check if extra variables were NOT computed before but are needed now
    bool needToRecompute = needExtraVariables && !lastComputedWithExtras_;

    if( !stateChanged && !needToRecompute )
    {
        return;  // Use cached values
    }

    // Update cache state
    cachedAltitude_ = altitude;
    cachedLongitude_ = longitude;
    cachedLatitude_ = latitude;
    cachedTime_ = time;
    lastComputedWithExtras_ = needExtraVariables;

    // Compute new values
    computePropertiesInternal( altitude, longitude, latitude, time, needExtraVariables );
}

void McdAtmosphereModel::computePropertiesInternal( const double altitude,
                                                    const double longitude,
                                                    const double latitude,
                                                    const double time,
                                                    const bool needExtraVariables )
{
    // Convert time to MCD format (Julian date)
    // MCD uses Julian date with dateKey=0 (Earth time)
    // localTime must be 0 when dateKey=0
    int dateKey = 0;
    double xdate = basic_astrodynamics::convertSecondsSinceEpochToJulianDay< double >( time );
    float localTime = 0.0f;

    // Convert coordinates to degrees
    float longitudeDeg = static_cast< float >( unit_conversions::convertRadiansToDegrees( longitude ) );
    float latitudeDeg = static_cast< float >( unit_conversions::convertRadiansToDegrees( latitude ) );

    // CRITICAL: Use zkey=3 for height above surface (matches Tudat convention)
    int zkey = 3;
    float altitudeInput = static_cast< float >( altitude );

    // Prepare extra variable flags
    int extvarkeys[ 100 ];
    for( int i = 0; i < 100; ++i )
    {
        if( i < 13 )
        {
            // Variables 1-13 are always computed (time/space coordinates)
            extvarkeys[ i ] = 1;
        }
        else if( i < 85 && needExtraVariables )
        {
            // Variables 14-85 only if requested
            extvarkeys[ i ] = 1;
        }
        else
        {
            // Variables 86-100 are unused
            extvarkeys[ i ] = 0;
        }
    }

    // Output variables (MCD uses float)
    float pres, dens, temp, zonwind, merwind;
    float meanvar[ 5 ];
    float extvar[ 100 ];
    float seedout;
    int ier;

    // Set localTime to 0 (required for dateKey=0)
    localTime = 0.0f;

    float seedin_f = static_cast< float >( perturbationSeed_ );
    float gwlength_f = static_cast< float >( gravityWaveLength_ );
    int dustScenario = dustScenario_;
    int perturbationKey = perturbationKey_;
    int highResolutionMode = highResolutionMode_;

    // Call MCD Fortran routine (interface declared in mcd.h)
    __mcd_MOD_call_mcd( &zkey,
                        &altitudeInput,
                        &longitudeDeg,
                        &latitudeDeg,
                        &highResolutionMode,
                        &dateKey,
                        &xdate,
                        &localTime,
                        mcdDataPath_.c_str( ),
                        &dustScenario,
                        &perturbationKey,
                        &seedin_f,
                        &gwlength_f,
                        extvarkeys,
                        &pres,
                        &dens,
                        &temp,
                        &zonwind,
                        &merwind,
                        meanvar,
                        extvar,
                        &seedout,
                        &ier,
                        static_cast< int >( mcdDataPath_.length( ) ) );

    // Check for errors
    if( ier != 0 )
    {
        throw std::runtime_error( "McdAtmosphereModel: MCD routine returned error code " + std::to_string( ier ) );
    }

    // Store basic results
    currentPressure_ = static_cast< double >( pres );
    currentDensity_ = static_cast< double >( dens );
    currentTemperature_ = static_cast< double >( temp );
    currentZonalWind_ = static_cast< double >( zonwind );
    currentMeridionalWind_ = static_cast< double >( merwind );
    currentSeedOut_ = static_cast< double >( seedout );

    // Copy mean variables
    for( int i = 0; i < 5; ++i )
    {
        currentMeanVariables_[ i ] = static_cast< double >( meanvar[ i ] );
    }

    // Copy extra variables (Fortran arrays are 1-indexed, C++ are 0-indexed)
    for( int i = 0; i < 100; ++i )
    {
        currentExtraVariables_[ i ] = static_cast< double >( extvar[ i ] );
    }

    // Extract and store always-computed supplementary variables (extvar 1-13 in Fortran)
    // Remember: Fortran arrays are 1-indexed, C++ arrays are 0-indexed
    if( currentExtraVariables_.size( ) >= 13 )
    {
        currentRadialDistance_ = currentExtraVariables_[ 0 ];        // extvar(1)
        currentAltitudeAboveAreoid_ = currentExtraVariables_[ 1 ];   // extvar(2)
        currentAltitudeAboveSurface_ = currentExtraVariables_[ 2 ];  // extvar(3)
        currentOrographicHeight_ = currentExtraVariables_[ 3 ];      // extvar(4)
        currentGcmOrography_ = currentExtraVariables_[ 4 ];          // extvar(5)
        currentSlopeInclination_ = currentExtraVariables_[ 5 ];      // extvar(6)
        currentSlopeOrientation_ = currentExtraVariables_[ 6 ];      // extvar(7)
        currentMarsAU_ = currentExtraVariables_[ 7 ];                // extvar(8)
        currentLs_ = currentExtraVariables_[ 8 ];                    // extvar(9)
        currentLTST_ = currentExtraVariables_[ 9 ];                  // extvar(10)
        currentLocalMeanTime_ = currentExtraVariables_[ 10 ];        // extvar(11)
        currentUniversalSolarTime_ = currentExtraVariables_[ 11 ];   // extvar(12)
        currentSolarZenithAngle_ = currentExtraVariables_[ 12 ];     // extvar(13)
    }
}
}  // namespace aerodynamics
}  // namespace tudat