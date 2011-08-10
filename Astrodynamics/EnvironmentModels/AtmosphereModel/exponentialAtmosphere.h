/*! \file exponentialAtmosphere.h
 *    Header file that defines the exponential atmosphere model included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/AtmosphereModel/
 *    Version           : 3
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
 *    Last modified     : 10 August, 2011
 *
 *    References
 *
 *    Notes
 *    Baseclass for standard and reference Atmospheres.
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110324    J. Melman         Added overloaded get functions.
 *      110427    F.M. Engelen      Changed input parameter to altitude, longitude and latitude.
 *      110629    F.M. Engelen      Added predefined function.
 *      110705    F.M. Engelen      Changed to passing by reference.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef EXPONENTIALATMOSPHERE_H
#define EXPONENTIALATMOSPHERE_H

// Include statements.
#include <cmath>
#include "atmosphereModel.h"
//#include "simpleAltitudeModel.h"

//! Exponential atmosphere class.
/*! Class with an exponential atmosphere. The user has to initialize this class
 * with a scale altitude, a density at zero meters altitude and a constant
 * temperature.
 */
class ExponentialAtmosphere : public AtmosphereModel
{
public:
    //! Bodies with predefined exponential atmospheres.
    /*!
     *  Bodies with predefined exponential atmospheres.
     */
    enum BodiesWithPredefinedExponentialAtmospheres
    {
        earth
    };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ExponentialAtmosphere( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ExponentialAtmosphere( );

    //! Set predefined exponential atmosphere settings.
    /*!
     * Set predefined exponential atmosphere settings
     * \param bodyWithPredefinedExponentialAtmosphere Body with a predefined exponential
     * atmosphere.
     */
    void setPredefinedExponentialAtmosphere( BodiesWithPredefinedExponentialAtmospheres
                                            bodyWithPredefinedExponentialAtmosphere );

    //! Set scale height.
    /*!
     * Set scale height (property of exponential atmosphere) in meters.
     */
    void setScaleHeight( const double& scaleHeight );

    //! Get scale height.
    /*!
     * Get scale height (property of exponential atmosphere) in meters.
     */
    double getScaleHeight( );

    //! Set density at zero altitude.
    /*!
     * Set density at zero altitude (property of exponential atmosphere) in kg per meter^3.
     */
    void setDensityAtZeroAltitude( const double& densityAtZeroAltitude );

    //! Get density at zero altitude.
    /*!
     * Get density at zero altitude (property of exponential atmosphere) in kg per meter^3.
     */
    double getDensityAtZeroAltitude( );

    //! Set constant temperature.
    /*!
     * Set the atmospheric temperature (constant, property of exponential atmosphere)
     * in Kelvin.
     */
    void setConstantTemperature( const double& constantTemperature );

    //! Get constant temperature.
    /*!
     * Get the atmospheric temperature (constant, property of exponential atmosphere)
     * in Kelvin.
     */
    double getConstantTemperature( );

    //! Get local density.
    /*!
     * Return the local density of the atmosphere in kg per meter^3.
     * It is determined by \f$ \rho = \rho_0 e^{-h/H} \f$, with \f$ \rho_0 \f$ the
     * density at zero altitude, \f$ h \f$ the altitude, and \f$ H \f$ the scaling
     * altitude.
     * \param altitude Altitude.
     * \return Atmospheric density.
     */
    double getDensity( const double& altitude );

    //! Get local density in the general way.
    /*!
     * Return the local density of the atmosphere in kg per meter^3.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return Atmospheric density.
     */
    double getDensity( const double& altitude,
                       const double& longitude,
                       const double& latitude,
                       const double& time = 0.0 );

    //! Get local pressure.
    /*!
     * Return the local pressure of the atmosphere in Newton per meter^2.
     * It is determined by the perfect gas law: \f$ p = \rho R T \f$, with
     * \f$ \rho \f$ the density, \f$ R \f$ the specific gas constant of air,
     * and \f$ T \f$ the temperature.
     * \param altitude Altitude.
     * \return Atmospheric pressure.
     */
    double getPressure( const double& altitude );

    //! Get local pressure in the general way.
    /*!
     * Return the local pressure of the atmosphere in Newton per meter^2.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return Atmospheric pressure.
     */
    double getPressure( const double& altitude,
                        const double& longitude,
                        const double& latitude,
                        const double& time = 0.0 );

    //! Get local temperature.
    /*!
     * Return the local temperature of the atmosphere in Kelvin.
     * \param altitude Altitude.
     * \return Atmospheric temperature.
     */
    double getTemperature( const double& altitude);

    //! Get local temperature.
    /*!
     * Return the local temperature of the atmosphere in Kelvin.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return Atmospheric temperature.
     */
    double getTemperature( const double& altitude,
                           const double& longitude,
                           const double& latitude,
                           const double& time = 0.0 );

protected:

private:

    //! Scale altitude.
    /*!
     * Scale altitude (property of exponential atmosphere) in meters.
     */
    double scaleHeight_;

    //! Constant temperature.
    /*!
     * The atmospheric temperature (constant, property of exponential atmosphere)
     * in Kelvin.
     */
    double constantTemperature_;

    //! density at zero altitude
    /*!
     * Density at zero altitude (property of exponential atmosphere) in kg per meter^3.
     */
    double densityAtZeroAltitude_;
};

#endif // EXPONENTIALATMOSPHERE_H

// End of file.
