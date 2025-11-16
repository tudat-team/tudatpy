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
    currentTemperature_( 0.0 ), currentZonalWind_( 0.0 ), currentMeridionalWind_( 0.0 ), currentSeedOut_( perturbationSeed )
{
    // Initialize mean and extra variables vectors
    currentMeanVariables_.resize( 5, 0.0 );
    currentExtraVariables_.resize( 100, 0.0 );

    // Set default MCD data path if not provided
    if( mcdDataPath_.empty( ) )
    {
        mcdDataPath_ = paths::getAtmosphereTablesPath( ) + "/MCD_DATA";
    }

    // Validate input parameters
    if( dustScenario_ < 1 || dustScenario_ > 8 )
    {
        throw std::runtime_error( "McdAtmosphereModel: Invalid dust scenario. Must be between 1 and 8." );
    }

    if( perturbationKey_ < 0 || perturbationKey_ > 5 )
    {
        throw std::runtime_error( "McdAtmosphereModel: Invalid perturbation key. Must be between 0 and 5." );
    }

    if( highResolutionMode_ != 0 && highResolutionMode_ != 1 )
    {
        throw std::runtime_error( "McdAtmosphereModel: Invalid high resolution mode. Must be 0 or 1." );
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

    // Calculate speed of sound using: a = sqrt(gamma * R * T)
    // For Mars atmosphere: gamma ≈ 1.3, R ≈ 192 J/(kg·K)
    const double gamma = 1.3;
    const double specificGasConstant = 192.0;  // J/(kg·K)

    return std::sqrt( gamma * specificGasConstant * currentTemperature_ );
}

// Compute properties (main computation function)
void McdAtmosphereModel::computeProperties( const double altitude, const double longitude, const double latitude, const double time )
{
    // Convert time to MCD format
    int dateKey;
    double xdate;
    double localTime;
    convertTimeToMcdFormat( time, longitude, dateKey, xdate, localTime );

    // Convert coordinates to degrees
    double longitudeDeg = unit_conversions::convertRadiansToDegrees( longitude );
    double latitudeDeg = unit_conversions::convertRadiansToDegrees( latitude );

    // Set vertical coordinate type: 3 = height above surface (meters)
    const int zkey = 3;

    // Prepare extra variable flags (all off for now)
    std::vector< int > extvarkeys( 100, 0 );

    // Output variables
    double pres, dens, temp, zonwind, merwind;
    int ier;

    // Call MCD Fortran routine
    callMcdFortran( zkey,
                    altitude,
                    longitudeDeg,
                    latitudeDeg,
                    highResolutionMode_,
                    dateKey,
                    xdate,
                    localTime,
                    mcdDataPath_,
                    dustScenario_,
                    perturbationKey_,
                    perturbationSeed_,
                    gravityWaveLength_,
                    extvarkeys,
                    pres,
                    dens,
                    temp,
                    zonwind,
                    merwind,
                    currentMeanVariables_,
                    currentExtraVariables_,
                    currentSeedOut_,
                    ier );

    // Check for errors
    if( ier != 0 )
    {
        throw std::runtime_error( "McdAtmosphereModel: MCD routine returned error code " + std::to_string( ier ) );
    }

    // Store results
    currentPressure_ = pres;
    currentDensity_ = dens;
    currentTemperature_ = temp;
    currentZonalWind_ = zonwind;
    currentMeridionalWind_ = merwind;
}

// Convert time to MCD format
void McdAtmosphereModel::convertTimeToMcdFormat( const double time, const double longitude, int& dateKey, double& xdate, double& localTime )
{
    // Convert time from seconds since J2000 to calendar date
    basic_astrodynamics::DateTime currentDateTime = basic_astrodynamics::DateTime::fromTime( time );

    // Use Earth date (dateKey = 0) and convert to Julian date
    dateKey = 0;

    // Convert to Julian date
    // PLACEHOLDER: This will call the actual MCD julian conversion routine
    double julianDate;
    int ier;

    // For now, use Tudat's conversion
    julianDate = basic_astrodynamics::convertCalendarDateToJulianDay( currentDateTime.getYear( ),
                                                                      currentDateTime.getMonth( ),
                                                                      currentDateTime.getDay( ),
                                                                      currentDateTime.getHour( ),
                                                                      currentDateTime.getMinute( ),
                                                                      static_cast< double >( currentDateTime.getSeconds( ) ) );

    xdate = julianDate;

    // Local time is set to 0 when using Earth date (compulsory with dateKey=0)
    localTime = 0.0;
}

// Call MCD Fortran routine (PLACEHOLDER)
void McdAtmosphereModel::callMcdFortran( const int zkey,
                                         const double xz,
                                         const double xlon,
                                         const double xlat,
                                         const int hireskey,
                                         const int datekey,
                                         const double xdate,
                                         const double localtime,
                                         const std::string& dset,
                                         const int dust,
                                         const int perturkey,
                                         const double seedin,
                                         const double gwlength,
                                         const std::vector< int >& extvarkeys,
                                         double& pres,
                                         double& dens,
                                         double& temp,
                                         double& zonwind,
                                         double& merwind,
                                         std::vector< double >& meanvar,
                                         std::vector< double >& extvar,
                                         double& seedout,
                                         int& ier )
{
    // PLACEHOLDER: This function will interface with the actual MCD Fortran routine
    // For now, provide dummy values for testing

    // TODO: Replace this with actual call to MCD Fortran routine:
    // __mcd_MOD_call_mcd(&zkey, &xz, &xlon, &xlat, &hireskey,
    //                    &datekey, &xdate, &localtime, dset.c_str(),
    //                    &dust, &perturkey, &seedin, &gwlength,
    //                    extvarkeys.data(), &pres, &dens, &temp,
    //                    &zonwind, &merwind, meanvar.data(),
    //                    extvar.data(), &seedout, &ier,
    //                    dset.length());

    // Dummy implementation for testing purposes
    // Using simple exponential atmosphere approximation
    const double scaleHeight = 11100.0;       // m
    const double surfacePressure = 610.0;     // Pa
    const double surfaceTemperature = 210.0;  // K
    const double surfaceDensity = 0.020;      // kg/m^3

    // Exponential decay
    double tempDecay = std::exp( -xz / ( 2.5 * scaleHeight ) );  // Slower decay
    temp = surfaceTemperature * ( 0.7 + 0.3 * tempDecay );       // Between 147K and 210K

    // Simple wind model (latitude dependent)
    zonwind = 10.0 * std::sin( xlat * mathematical_constants::PI / 180.0 );
    merwind = 5.0 * std::cos( xlat * mathematical_constants::PI / 180.0 );

    // No errors
    ier = 0;
    seedout = seedin;

    // Note: meanvar and extvar remain as initialized (zeros)

    std::cout << "McdAtmosphereModel: Using PLACEHOLDER MCD routine. "
              << "Replace with actual MCD interface." << std::endl;
}

}  // namespace aerodynamics
}  // namespace tudat