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
 * Input altitude is expected as "height above local surface" (matching Tudat convention).
 * Internally converted to radial distance: radial_distance = MARS_MEAN_RADIUS + altitude
 * where MARS_MEAN_RADIUS = 3396200.0 m (IAU 2015).
 *
 * LIMITATION: This simplified conversion does not account for:
 *   - Local areoid variations (Mars oblate shape)
 *   - Local MOLA topography (when highResolutionMode=1)
 * This causes ~10-15% systematic differences vs. MCD reference outputs.
 *
 * TODO: Enhance to use proper Body shape model + MOLA topography corrections.
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
     * Note: Internally uses zkey=1 (radial distance from center). The input altitude
     *       is converted to radial distance using a fixed Mars mean radius (3396.2 km).
     *       TODO: Should be improved to use the actual Body shape model.
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
     * Note: Internally uses zkey=1 (radial distance from center). The input altitude
     *       is converted to radial distance using a fixed Mars mean radius (3396.2 km).
     *       TODO: Should be improved to use the actual Body shape model.
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
     * Note: Internally uses zkey=1 (radial distance from center). The input altitude
     *       is converted to radial distance using a fixed Mars mean radius (3396.2 km).
     *       TODO: Should be improved to use the actual Body shape model.
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
     * Note: Internally uses zkey=1 (radial distance from center). The input altitude
     *       is converted to radial distance using a fixed Mars mean radius (3396.2 km).
     *       TODO: Should be improved to use the actual Body shape model.
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

protected:
private:
    //! Compute atmospheric properties
    void computeProperties( const double altitude, const double longitude, const double latitude, const double time );

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

    //! Cached solar longitude (degrees)
    double currentLs_;

    //! Cached local true solar time (hours)
    double currentLTST_;

    //! Cached Sun-Mars distance (AU)
    double currentMarsAU_;
};

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_MCD_ATMOSPHERE_MODEL_H