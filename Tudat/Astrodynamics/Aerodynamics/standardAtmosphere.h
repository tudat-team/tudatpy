/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_STANDARD_ATMOSPHERE_H
#define TUDAT_STANDARD_ATMOSPHERE_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"

namespace tudat
{
namespace aerodynamics
{

//! Standard atmosphere model class.
/*!
 * Base class for all standard atmosphere models. These models are all independent of time and
 * position, and therefore can return the density, pressure and temperature as a function of
 * altitude only.
 */
class StandardAtmosphere : public AtmosphereModel
{
public:

    //! Get local density.
    /*!
    * Returns the local density parameter of the atmosphere in kg-per-meter^3.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric density.
    */
    virtual double getDensity( const double altitude, const double longitude = 0.0,
                               const double latitude = 0.0, const double time = 0.0 ) = 0;

    //! Get local pressure.
    /*!
    * Returns the local pressure of the atmosphere parameter in Newton-per-meter^2.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric pressure.
    */
    virtual double getPressure( const double altitude, const double longitude = 0.0,
                                const double latitude = 0.0, const double time = 0.0 ) = 0;

    //! Get local temperature.
    /*!
    * Returns the local temperature of the atmosphere parameter in Kelvin.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric temperature.
    */
    virtual double getTemperature( const double altitude, const double longitude = 0.0,
                                   const double latitude = 0.0, const double time = 0.0 ) = 0;

    //! Get local speed of sound.
    /*!
    * Returns the local speed of sound of the atmosphere in m/s.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric speed of sound.
    */
    virtual double getSpeedOfSound( const double altitude, const double longitude = 0.0,
                                    const double latitude = 0.0, const double time = 0.0 ) = 0;
};

//! Typedef for shared-pointer to StandardAtmosphere object.
typedef boost::shared_ptr< StandardAtmosphere > StandardAtmospherePointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_STANDARD_ATMOSPHERE_H
