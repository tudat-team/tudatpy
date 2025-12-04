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
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/io/basicInputOutput.h"
#include <cmath>
#include <stdexcept>

// Include MCD header
#include "mcd.h"

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
    currentLs_( 0.0 ), currentLTST_( 0.0 ), currentMarsAU_( 0.0 )
{
    // Initialize mean and extra variables vectors
    currentMeanVariables_.resize( 5, 0.0 );
    currentExtraVariables_.resize( 100, 0.0 );

    // Set default MCD data path if not provided by user
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

    // Ensure the path ends with a slash
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
    computeProperties( altitude, longitude, latitude, time );
    return currentDensity_;
}

// Get pressure
double McdAtmosphereModel::getPressure( const double altitude, const double longitude, const double latitude, const double time )
{
    computeProperties( altitude, longitude, latitude, time );
    return currentPressure_;
}

// Get temperature
double McdAtmosphereModel::getTemperature( const double altitude, const double longitude, const double latitude, const double time )
{
    computeProperties( altitude, longitude, latitude, time );
    return currentTemperature_;
}

// Get speed of sound
double McdAtmosphereModel::getSpeedOfSound( const double altitude, const double longitude, const double latitude, const double time )
{
    computeProperties( altitude, longitude, latitude, time );

    // Get gamma and R from extra variables if available
    // extvar(60) = gamma, extvar(61) = R
    double gamma = currentExtraVariables_[ 59 ];  // Index 59 = extvar(60)
    double R = currentExtraVariables_[ 60 ];      // Index 60 = extvar(61)

    if( gamma > 0.0 && R > 0.0 )
    {
        return std::sqrt( gamma * R * currentTemperature_ );
    }
    else
    {
        // Fallback values for Mars atmosphere
        const double defaultGamma = 1.3;
        const double defaultR = 192.0;  // J/(kg·K)
        return std::sqrt( defaultGamma * defaultR * currentTemperature_ );
    }
}

// Compute properties (main computation function)
void McdAtmosphereModel::computeProperties( const double altitude, const double longitude, const double latitude, const double time )
{
    // Convert time to MCD format (Julian date)
    int dateKey = 0;  // Use Earth date (Julian date)
    double xdate;
    float localTime;

    convertTimeToMcdFormat( time, longitude, dateKey, xdate, localTime );

    // Convert coordinates to degrees
    float longitudeDeg = static_cast< float >( unit_conversions::convertRadiansToDegrees( longitude ) );
    float latitudeDeg = static_cast< float >( unit_conversions::convertRadiansToDegrees( latitude ) );

    // Use zkey=2: height above areoid (meters)
    // Tudat's altitude is typically height above the reference ellipsoid (oblate spheroid).
    // The Areoid is the closest approximation to this reference surface in MCD.
    // Using zkey=3 would incorrectly add local topography height to the input altitude.
    int zkey = 2;  // height above areoid
    float altitudeAboveSurface = static_cast< float >( altitude );

    // Prepare extra variable flags (enable all for comprehensive output)
    int extvarkeys[ 100 ];
    for( int i = 0; i < 100; ++i )
    {
        extvarkeys[ i ] = 1;
    }

    // Output variables (MCD uses float)
    float pres, dens, temp, zonwind, merwind;
    float meanvar[ 5 ];
    float extvar[ 100 ];
    float seedout;
    int ier;

    // localTime must be 0 when using dateKey=0 (compulsory)
    localTime = 0.0f;

    // Prepare other inputs as float (remove const to match Fortran interface)
    float seedin_f = static_cast< float >( perturbationSeed_ );
    float gwlength_f = static_cast< float >( gravityWaveLength_ );
    int dustScenario = dustScenario_;
    int perturbationKey = perturbationKey_;
    int highResolutionMode = highResolutionMode_;

    // Call MCD Fortran routine directly
    __mcd_MOD_call_mcd( &zkey,
                        &altitudeAboveSurface,
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

    // Store results (convert float to double)
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

    // Copy extra variables
    for( int i = 0; i < 100; ++i )
    {
        currentExtraVariables_[ i ] = static_cast< double >( extvar[ i ] );
    }

    // Extract Ls and other orbital parameters from extra variables
    // According to MCD documentation:
    // extvar(9) = Ls (solar longitude in degrees)
    // extvar(10) = LST (local true solar time in hours)
    // Note: array indices in C++ are 0-based, so extvar[8] = extvar(9) in Fortran
    if( currentExtraVariables_.size( ) >= 10 )
    {
        currentLs_ = currentExtraVariables_[ 8 ];    // extvar(9) = Ls
        currentLTST_ = currentExtraVariables_[ 9 ];  // extvar(10) = LST
    }
}

// Convert time to MCD format (simplified - only Julian date needed)
void McdAtmosphereModel::convertTimeToMcdFormat( const double time, const double longitude, int& dateKey, double& xdate, float& localTime )
{
    // Convert time from seconds since J2000 to Julian date
    const double J2000_JD = 2451545.0;  // Julian date of J2000 epoch
    const double secondsPerDay = 86400.0;

    double julianDate = J2000_JD + ( time / secondsPerDay );

    // Use Earth date (dateKey = 0)
    dateKey = 0;
    xdate = julianDate;

    // localTime must be set to 0 when using dateKey=0 (compulsory)
    localTime = 0.0f;
}

}  // namespace aerodynamics
}  // namespace tudat