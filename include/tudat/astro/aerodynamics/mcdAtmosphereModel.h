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
 */
class McdAtmosphereModel : public AtmosphereModel
{
public:
    //! Constructor
    /*!
     * Constructor for MCD atmosphere model.
     * \param mcdDataPath Path to MCD data files
     * \param dustScenario Dust and solar EUV scenario (1-8)
     * \param perturbationKey Perturbation type (0: none, 1-4: various perturbations, 5: n-sigma)
     * \param perturbationSeed Random seed for perturbations (if perturbationKey != 0)
     * \param gravityWaveLength Gravity wave wavelength (if perturbationKey == 3 or 4)
     * \param highResolutionMode Flag for high resolution topography (0: off, 1: on)
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
     * \param altitude Altitude above surface (m)
     * \param longitude East longitude (radians)
     * \param latitude Latitude (radians)
     * \param time Time since J2000 (seconds)
     * \return Atmospheric density (kg/m^3)
     */
    virtual double getDensity( const double altitude, const double longitude, const double latitude, const double time );

    //! Get local pressure
    /*!
     * Returns the local pressure of the atmosphere in Pa.
     * \param altitude Altitude above surface (m)
     * \param longitude East longitude (radians)
     * \param latitude Latitude (radians)
     * \param time Time since J2000 (seconds)
     * \return Atmospheric pressure (Pa)
     */
    virtual double getPressure( const double altitude, const double longitude, const double latitude, const double time );

    //! Get local temperature
    /*!
     * Returns the local temperature of the atmosphere in K.
     * \param altitude Altitude above surface (m)
     * \param longitude East longitude (radians)
     * \param latitude Latitude (radians)
     * \param time Time since J2000 (seconds)
     * \return Atmospheric temperature (K)
     */
    virtual double getTemperature( const double altitude, const double longitude, const double latitude, const double time );

    //! Get local speed of sound
    /*!
     * Returns the local speed of sound of the atmosphere in m/s.
     * \param altitude Altitude above surface (m)
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

protected:
private:
    //! Compute atmospheric properties
    /*!
     * Internal function to compute atmospheric properties at given conditions.
     * This function calls the MCD Fortran routine and caches results.
     * \param altitude Altitude above surface (m)
     * \param longitude East longitude (radians)
     * \param latitude Latitude (radians)
     * \param time Time since J2000 (seconds)
     */
    void computeProperties( const double altitude, const double longitude, const double latitude, const double time );

    //! Convert time to MCD format
    /*!
     * Converts time from seconds since J2000 to MCD date format.
     * \param time Time since J2000 (seconds)
     * \param dateKey Date key (0: Earth date, 1: Mars date with Ls)
     * \param xdate Output date (Julian date if dateKey=0, Ls if dateKey=1)
     * \param localTime Output local time (hours)
     */
    void convertTimeToMcdFormat( const double time, const double longitude, int& dateKey, double& xdate, double& localTime );

    //! Call MCD Fortran routine
    /*!
     * Placeholder function to call the MCD Fortran routine.
     * This will be implemented with actual MCD interface later.
     * \param zkey Vertical coordinate type (1-4)
     * \param xz Vertical coordinate value
     * \param xlon East longitude (degrees)
     * \param xlat Latitude (degrees)
     * \param hireskey High resolution flag
     * \param datekey Date flag (0: Earth date, 1: Mars date)
     * \param xdate Date value
     * \param localtime Local time (hours)
     * \param dset Path to datasets
     * \param dust Dust scenario
     * \param perturkey Perturbation type
     * \param seedin Perturbation seed
     * \param gwlength Gravity wave wavelength
     * \param extvarkeys Extra variable flags
     * \param pres Output pressure
     * \param dens Output density
     * \param temp Output temperature
     * \param zonwind Output zonal wind
     * \param merwind Output meridional wind
     * \param meanvar Output mean variables
     * \param extvar Output extra variables
     * \param seedout Output seed
     * \param ier Output error code
     */
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

    //! Dust and solar EUV scenario (1-8)
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
};

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_MCD_ATMOSPHERE_MODEL_H