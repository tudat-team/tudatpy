/*! \file atmosphereModel.h
 *    Header file that defines the baseclass atmosphere model included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/AtmosphereModel/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 10 March, 2011
 *    Last modified     : 05 July, 2011
 *
 *    References
 *
 *    Notes
 *      Baseclass for standard and reference Atmospheres.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110224    F.M. Engelen      File created.
 *      110324    J. Melman         Set time to zero by default.
 *      110412    F.M. Engelen      Included flag boolean isAltitudeModelDefault_.
 *      110427    F.M. Engelen      Changed input arguments to altitude, longitude, latitude.
 *      110705    F.M. Engelen      Changed to passing by reference.
 */

#ifndef ATMOSPHEREMODEL_H
#define ATMOSPHEREMODEL_H

// Include statements.
#include "environmentModel.h"
#include "physicalConstants.h"

//! Atmosphere model class.
/*! Base class for all atmosphere models.
 * To keep the function generic for both reference atmospheres and standard
 * atmospheres, altitude, longitude, latitude, and time are inputs for all functions.
 */

class AtmosphereModel : public EnvironmentModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    AtmosphereModel( );

    //! Default destructor.
    /*!
    * Default destructor.
    */
    ~AtmosphereModel( );

    //! Get local density.
    /*!
    * Return the local density parameter of the atmosphere in kg per meter^3.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric density.
    */
    virtual double getDensity( const double& altitude,
                               const double& longitude,
                               const double& latitude,
                               const double& time = 0.0 ) = 0;

    //! Get local pressure.
    /*!
    * Return the local pressure of the atmosphere parameter in Newton per meter^2.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric pressure.
    */
    virtual double getPressure( const double& altitude,
                                const double& longitude,
                                const double& latitude,
                                const double& time = 0.0 ) = 0;

    //! Get local temperature.
    /*!
    * Return the local temperature of the atmosphere parameter in Kelvin.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric temperature.
    */
    virtual double getTemperature( const double& altitude,
                                   const double& longitude,
                                   const double& latitude,
                                   const double& time = 0.0 ) = 0;

protected:

private:
};

#endif // ATMOSPHEREMODEL_H

// End of file.
