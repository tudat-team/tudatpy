/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MCD_ATMOSPHERE_MODEL_H
#define TUDAT_MCD_ATMOSPHERE_MODEL_H

#include <iostream>
#include <vector>
#include <map>
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/dateTime.h"

namespace tudat
{
namespace aerodynamics
{

//! MCD Atmosphere Model class
/*!
 * Class for Mars Climate Database atmosphere model.
 * This class interfaces with the MCD Fortran routines to provide
 * atmospheric properties (density, temperature, pressure, winds) for Mars.
 *
 * ALTITUDE CONVENTION:
 * --------------------
 * Input altitude is "height above local surface" (matching Tudat convention).
 * Internally uses MCD's zkey=3 mode, which allows MCD to handle the conversion
 * to radial distance using its own areoid and topography models.
 *
 * When highResolutionMode=1, MCD uses MOLA topography for accurate surface height.
 * When highResolutionMode=0, MCD uses GCM resolution topography.
 *
 * PARAMETERS:
 * -----------
 * See constructor documentation for detailed parameter descriptions.
 *
 * THREAD SAFETY:
 * --------------
 * This class is NOT thread-safe due to cached internal state.
 * Each thread should use its own instance.
 */
class McdAtmosphereModel : public AtmosphereModel
{
public:
    //! Constructor
    /*!
     * Constructor for MCD atmosphere model.
     * \param mcdDataPath Path to MCD data files (default: "" = use compile-time default)
     * \param dustScenario Dust and solar EUV scenario (1-8 or 24-35, default: 1)
     * \param perturbationKey Perturbation type (0-5, default: 0 = none)
     * \param perturbationSeed Random seed or scaling factor (default: 0.0)
     * \param gravityWaveLength Gravity wave wavelength in meters (default: 0.0 = use MCD default)
     * \param highResolutionMode High resolution topography flag (0 or 1, default: 0)
     */
    McdAtmosphereModel( const std::string& mcdDataPath = "",
                        const int dustScenario = 1,
                        const int perturbationKey = 0,
                        const double perturbationSeed = 0.0,
                        const double gravityWaveLength = 0.0,
                        const int highResolutionMode = 0 );

    //! Destructor
    virtual ~McdAtmosphereModel( ) {}

    //! Get local density
    /*!
     * Returns the local density of the atmosphere in kg/m^3.
     * \param altitude Altitude above local surface (m) - as computed by Tudat
     * \param longitude East longitude (radians)
     * \param latitude Latitude (radians)
     * \param time Time since J2000 (seconds)
     * \return Atmospheric density (kg/m^3)
     */
    virtual double getDensity( const double altitude, const double longitude, const double latitude, const double time );

    //! Get local pressure
    /*!
     * Returns the local pressure of the atmosphere in Pa.
     * \param altitude Altitude above local surface (m) - as computed by Tudat
     * \param longitude East longitude (radians)
     * \param latitude Latitude (radians)
     * \param time Time since J2000 (seconds)
     * \return Atmospheric pressure (Pa)
     */
    virtual double getPressure( const double altitude, const double longitude, const double latitude, const double time );

    //! Get local temperature
    /*!
     * Returns the local temperature of the atmosphere in K.
     * \param altitude Altitude above local surface (m) - as computed by Tudat
     * \param longitude East longitude (radians)
     * \param latitude Latitude (radians)
     * \param time Time since J2000 (seconds)
     * \return Atmospheric temperature (K)
     */
    virtual double getTemperature( const double altitude, const double longitude, const double latitude, const double time );

    //! Get local speed of sound
    /*!
     * Returns the local speed of sound of the atmosphere in m/s.
     * \param altitude Altitude above local surface (m) - as computed by Tudat
     * \param longitude East longitude (radians)
     * \param latitude Latitude (radians)
     * \param time Time since J2000 (seconds)
     * \return Atmospheric speed of sound (m/s)
     */
    virtual double getSpeedOfSound( const double altitude, const double longitude, const double latitude, const double time );

    //! Get zonal wind component (u)
    /*!
     * Returns the zonal (east-west) wind component.
     * \return Zonal wind (m/s)
     */
    double getZonalWind( ) const
    {
        return currentZonalWind_;
    }

    //! Get meridional wind component (v)
    /*!
     * Returns the meridional (north-south) wind component.
     * \return Meridional wind (m/s)
     */
    double getMeridionalWind( ) const
    {
        return currentMeridionalWind_;
    }

    //! Get mean variables (pressure, density, temperature, zonal wind, meridional wind)
    const std::vector< double >& getMeanVariables( ) const
    {
        return currentMeanVariables_;
    }

    //! Get extra variables (100 variables as per MCD specification)
    const std::vector< double >& getExtraVariables( ) const
    {
        return currentExtraVariables_;
    }

    //! Get dust scenario
    int getDustScenario( ) const
    {
        return dustScenario_;
    }

    //! Set dust scenario
    void setDustScenario( const int dustScenario )
    {
        dustScenario_ = dustScenario;
    }

    //! Get high resolution mode
    int getHighResolutionMode( ) const
    {
        return highResolutionMode_;
    }

    //! Set high resolution mode
    void setHighResolutionMode( const int highResolutionMode )
    {
        highResolutionMode_ = highResolutionMode;
    }

    //! Get solar longitude (Ls) in degrees
    double getSolarLongitude( ) const
    {
        return currentLs_;
    }

    //! Get local true solar time in hours
    double getLocalTrueSolarTime( ) const
    {
        return currentLTST_;
    }

    //! Get radial distance from planet center (always available)
    double getRadialDistance( ) const
    {
        return currentRadialDistance_;
    }

    //! Get altitude above areoid (always available)
    double getAltitudeAboveAreoid( ) const
    {
        return currentAltitudeAboveAreoid_;
    }

    //! Get altitude above local surface (always available)
    double getAltitudeAboveSurface( ) const
    {
        return currentAltitudeAboveSurface_;
    }

    //! Get orographic height (always available)
    double getOrographicHeight( ) const
    {
        return currentOrographicHeight_;
    }

    //! Get GCM orography (always available)
    double getGcmOrography( ) const
    {
        return currentGcmOrography_;
    }

    //! Get local slope inclination in degrees (always available if hireskey=1)
    double getSlopeInclination( ) const
    {
        return currentSlopeInclination_;
    }

    //! Get local slope orientation in degrees (always available if hireskey=1)
    double getSlopeOrientation( ) const
    {
        return currentSlopeOrientation_;
    }

    //! Get Sun-Mars distance in AU (always available)
    double getSunMarsDistance( ) const
    {
        return currentMarsAU_;
    }

    //! Get local mean time in hours (always available)
    double getLocalMeanTime( ) const
    {
        return currentLocalMeanTime_;
    }

    //! Get universal solar time in hours (always available)
    double getUniversalSolarTime( ) const
    {
        return currentUniversalSolarTime_;
    }

    //! Get solar zenith angle in degrees (always available)
    double getSolarZenithAngle( ) const
    {
        return currentSolarZenithAngle_;
    }

protected:
    //! Check if current state matches cached state
    /*!
     * Checks if the provided state parameters match the cached state within tolerance.
     * \param altitude Altitude to check (m)
     * \param longitude Longitude to check (rad)
     * \param latitude Latitude to check (rad)
     * \param time Time to check (seconds since J2000)
     * \return True if state matches cached state
     */
    bool isCached( const double altitude, const double longitude, const double latitude, const double time ) const;

    //! Compute properties only if state changed or extras needed
    /*!
     * Checks cached state and only recomputes if necessary.
     * \param altitude Altitude above surface (m) - matches Tudat convention
     * \param longitude Longitude (rad)
     * \param latitude Latitude (rad)
     * \param time Time (seconds since J2000)
     * \param needExtraVariables Whether extra variables 14-100 are needed
     */
    void computePropertiesIfNeeded( const double altitude,
                                    const double longitude,
                                    const double latitude,
                                    const double time,
                                    const bool needExtraVariables = false );

private:
    //! Compute atmospheric properties (internal implementation)
    /*!
     * Internal method that actually calls MCD.
     * Uses zkey=3 (height above surface) to match Tudat's altitude convention.
     * \param altitude Altitude above surface (m)
     * \param longitude Longitude (rad)
     * \param latitude Latitude (rad)
     * \param time Time (seconds since J2000)
     * \param needExtraVariables Whether to compute extra variables 14-100
     */
    void computePropertiesInternal( const double altitude,
                                    const double longitude,
                                    const double latitude,
                                    const double time,
                                    const bool needExtraVariables );

    //! Convert time to MCD format (Julian date only)
    /*!
     * Converts time from seconds since J2000 to Julian date.
     * \param time Time since J2000 (seconds)
     * \param longitude East longitude (radians) - not used with dateKey=0
     * \param dateKey Date key (always 0 for Julian date)
     * \param xdate Output Julian date
     * \param localTime Output local time (always 0 with dateKey=0)
     */
    void convertTimeToMcdFormat( const double time, const double longitude, int& dateKey, double& xdate, float& localTime );

    //! Call MCD Fortran routine via C interface
    void callMcdFortran( const int zkey,
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
                         int& ier );

    //! Path to MCD data files
    std::string mcdDataPath_;

    //! Dust and solar EUV scenario (1-8 or 24-35)
    int dustScenario_;

    //! Perturbation type
    int perturbationKey_;

    //! Perturbation seed
    double perturbationSeed_;

    //! Gravity wave wavelength
    double gravityWaveLength_;

    //! High resolution mode flag
    int highResolutionMode_;

    //! Cached atmospheric density (kg/m^3)
    double currentDensity_;

    //! Cached atmospheric pressure (Pa)
    double currentPressure_;

    //! Cached atmospheric temperature (K)
    double currentTemperature_;

    //! Cached zonal wind (m/s)
    double currentZonalWind_;

    //! Cached meridional wind (m/s)
    double currentMeridionalWind_;

    //! Cached mean variables
    std::vector< double > currentMeanVariables_;

    //! Cached extra variables
    std::vector< double > currentExtraVariables_;

    //! Current random seed value
    double currentSeedOut_;

    //! Cached state variables to avoid redundant MCD calls
    double cachedAltitude_ = TUDAT_NAN;
    double cachedLongitude_ = TUDAT_NAN;
    double cachedLatitude_ = TUDAT_NAN;
    double cachedTime_ = TUDAT_NAN;
    bool lastComputedWithExtras_ = false;

    //! Always-computed supplementary variables (extvar 1-13 in Fortran)
    double currentRadialDistance_;        // extvar(1): Radial distance from planet center (m)
    double currentAltitudeAboveAreoid_;   // extvar(2): Altitude above areoid (m)
    double currentAltitudeAboveSurface_;  // extvar(3): Altitude above local surface (m)
    double currentOrographicHeight_;      // extvar(4): Orographic height (m)
    double currentGcmOrography_;          // extvar(5): GCM orography (m)
    double currentSlopeInclination_;      // extvar(6): Local slope inclination (deg)
    double currentSlopeOrientation_;      // extvar(7): Local slope orientation (deg)
    double currentMarsAU_;                // extvar(8): Sun-Mars distance (AU)
    double currentLs_;                    // extvar(9): Solar longitude Ls (deg)
    double currentLTST_;                  // extvar(10): Local true solar time (hrs)
    double currentLocalMeanTime_;         // extvar(11): Local mean time (hrs)
    double currentUniversalSolarTime_;    // extvar(12): Universal solar time (hrs)
    double currentSolarZenithAngle_;      // extvar(13): Solar zenith angle (deg)
};

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_MCD_ATMOSPHERE_MODEL_H