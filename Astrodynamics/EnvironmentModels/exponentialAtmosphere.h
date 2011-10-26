/*! \file exponentialAtmosphere.h
 *    Header file that defines the exponential atmosphere model included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
 *    Version           : 6
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
 *      110324    J. Melman         Added overloaded get functions.
 *      110427    F.M. Engelen      Changed input parameter to altitude, longitude and latitude.
 *      110629    F.M. Engelen      Added predefined function.
 *      110705    F.M. Engelen      Changed to passing by reference.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef EXPONENTIALATMOSPHERE_H
#define EXPONENTIALATMOSPHERE_H

// Macros.
#define TUDAT_UNUSED_PARAMETER( unusedParameter ) { ( void ) unusedParameter; }

// Include statements.
#include <cmath>
#include "Astrodynamics/EnvironmentModels/atmosphereModel.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Exponential atmosphere class.
/*!
 * Class with an exponential atmosphere. The user has to initialize this class
 * with a scale altitude, a density at zero meters altitude and a constant
 * temperature. The density is determined by \f$ \rho = \rho_0 e^{-h/H} \f$, with \f$ \rho_0 \f$
 * the density at zero altitude, \f$ h \f$ the altitude, and \f$ H \f$ the scaling
 * altitude. The temperature is taken as constant and the pressure follows from the universal
 * gas law \f$ p = \rho RT$
 */
class ExponentialAtmosphere : public AtmosphereModel
{
public:

    //! Bodies with predefined exponential atmospheres.
    /*!
     *  Bodies with predefined exponential atmospheres.
     */
    enum BodiesWithPredefinedExponentialAtmospheres { earth };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ExponentialAtmosphere( ) : scaleHeight_( -0.0 ), constantTemperature_( -0.0 ),
        densityAtZeroAltitude_( -0.0 ), specificGasConstant_( -0.0 ) { }

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
     * Sets the scale height (property of exponential atmosphere) in meters.
     */
    void setScaleHeight( const double& scaleHeight ) { scaleHeight_ = scaleHeight; }

    //! Get scale height.
    /*!
     * Returns the scale height (property of exponential atmosphere) in meters.
     */
    double getScaleHeight( ) { return scaleHeight_; }

    //! Set density at zero altitude.
    /*!
     * Sets the density at zero altitude (property of exponential atmosphere) in kg per meter^3.
     */
    void setDensityAtZeroAltitude( const double& densityAtZeroAltitude )
    { densityAtZeroAltitude_ = densityAtZeroAltitude; }

    //! Get density at zero altitude.
    /*!
     * Returns the density at zero altitude (property of exponential atmosphere) in kg per meter^3.
     */
    double getDensityAtZeroAltitude( ) { return densityAtZeroAltitude_; }

    //! Set constant temperature.
    /*!
     * Sets the atmospheric temperature (constant, property of exponential atmosphere) in Kelvin.
     */
    void setConstantTemperature( const double& constantTemperature )
    { constantTemperature_ = constantTemperature; }

    //! Get constant temperature.
    /*!
     * Returns the atmospheric temperature (constant, property of exponential atmosphere) in
     * Kelvin.
     */
    double getConstantTemperature( ) { return constantTemperature_; }

    //! Set specific gas constant.
    /*!
     * Sets the specific gas constant of the air in J/(kg K), its value is assumed constant,
     * due to the assumption of constant atmospheric composition.
     */
    void setSpecificGasConstant( const double& specificGasConstant )
    { specificGasConstant_ = specificGasConstant; }

    //! Get specific gas constant.
    /*!
     * Returns the specific gas constant of the air in J/(kg K), its value is assumed constant,
     * due to the assumption of constant atmospheric composition.
     */
    double getSpecificGasConstant( ) { return specificGasConstant_; }

    //! Get local density.
    /*!
     * Returns the local density of the atmosphere in kg per meter^3.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return Atmospheric density.
     */
    double getDensity( const double& altitude, const double& longitude = 0.0,
                       const double& latitude = 0.0, const double& time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return densityAtZeroAltitude_ * std::exp( -altitude / scaleHeight_ );
    }

    //! Get local pressure.
    /*!
     * Returns the local pressure of the atmosphere in Newton per meter^2.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return Atmospheric pressure.
     */
    double getPressure( const double& altitude, const double& longitude = 0.0,
                        const double& latitude = 0.0, const double& time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return getDensity( altitude ) * specificGasConstant_ * constantTemperature_;
    }

    //! Get local temperature.
    /*!
     * Returns the local temperature of the atmosphere in Kelvin.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return Atmospheric temperature.
     */
    double getTemperature( const double& altitude, const double& longitude = 0.0,
                           const double& latitude = 0.0, const double& time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( altitude );
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return constantTemperature_;
    }

protected:

private:

    //! Scale altitude.
    /*!
     * Scale altitude (property of exponential atmosphere) in meters.
     */
    double scaleHeight_;

    //! Constant temperature.
    /*!
     * The atmospheric temperature (constant, property of exponential atmosphere) in Kelvin.
     */
    double constantTemperature_;

    //! density at zero altitude
    /*!
     * Density at zero altitude (property of exponential atmosphere) in kg per meter^3.
     */
    double densityAtZeroAltitude_;

    //! Specific gas constant.
    /*!
     * Specific gas constant of the air, its value is assumed constant, due to the assumption of
     * constant atmospheric composition.
     */
    double specificGasConstant_;
};

}

#endif // EXPONENTIALATMOSPHERE_H

// End of file.
