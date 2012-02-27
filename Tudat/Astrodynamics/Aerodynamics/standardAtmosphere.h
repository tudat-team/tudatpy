/*    Copyright (c) 2010 - 2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110224    B. Tong Minh      File created.
 *      111211    K. Kumar          Minor corrections; corrected const-reference to by-value.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// Baseclass for standard Atmospheres.
// 

#ifndef TUDAT_STANDARD_ATMOSPHERE_H
#define TUDAT_STANDARD_ATMOSPHERE_H

#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"

namespace tudat
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
    virtual getDensity( double altitude, double longitude = 0.0,
                        double latitude = 0.0, double time = 0.0 ) = 0;

    //! Get local pressure.
    /*!
    * Returns the local pressure of the atmosphere parameter in Newton-per-meter^2.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric pressure.
    */
    virtual double getPressure( double altitude, double longitude = 0.0,
                                double latitude = 0.0, double time = 0.0 ) = 0;

    //! Get local temperature.
    /*!
    * Returns the local temperature of the atmosphere parameter in Kelvin.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric temperature.
    */
    virtual double getTemperature( double altitude, double longitude = 0.0,
                                   double latitude = 0.0, double time = 0.0 ) = 0;
};

} // namespace tudat

#endif // TUDAT_STANDARD_ATMOSPHERE_H
