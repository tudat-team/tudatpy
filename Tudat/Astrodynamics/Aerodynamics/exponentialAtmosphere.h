/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110224    F.M. Engelen      File created.
 *      110324    J. Melman         Added overloaded get functions.
 *      110427    F.M. Engelen      Changed input parameter to altitude, longitude and latitude.
 *      110629    F.M. Engelen      Added predefined function.
 *      110705    F.M. Engelen      Changed to passing by reference.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_EXPONENTIAL_ATMOSPHERE_H
#define TUDAT_EXPONENTIAL_ATMOSPHERE_H

#include <cmath>

#include <TudatCore/Basics/utilityMacros.h>

#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"

namespace tudat
{
namespace aerodynamics
{

//! Exponential atmosphere class.
/*!
 * Class with an exponential atmosphere. The user has to initialize this class
 * with a scale altitude, a density at zero meters altitude and a constant
 * temperature. The density is determined by \f$ \rho = \rho_0 e^{-h/H} \f$, with \f$ \rho_0 \f$
 * the density at zero altitude, \f$ h \f$ the altitude, and \f$ H \f$ the scaling
 * altitude. The temperature is taken as constant and the pressure follows from the universal
 * gas law \f$ p = \rho RT \f$.
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
    ExponentialAtmosphere( )
        : scaleHeight_( -0.0 ),
          constantTemperature_( -0.0 ),
          densityAtZeroAltitude_( -0.0 ),
          specificGasConstant_( -0.0 )
    { }

    //! Set predefined exponential atmosphere settings.
    /*!
     * Sets predefined exponential atmosphere settings.
     * \param bodyWithPredefinedExponentialAtmosphere Body with a predefined exponential
     *          atmosphere.
     */
    void setPredefinedExponentialAtmosphere( BodiesWithPredefinedExponentialAtmospheres
                                             bodyWithPredefinedExponentialAtmosphere );

    //! Set scale height.
    /*!
     * Sets the scale height (property of exponential atmosphere) in meters.
     */
    void setScaleHeight( double scaleHeight ) { scaleHeight_ = scaleHeight; }

    //! Get scale height.
    /*!
     * Returns the scale height (property of exponential atmosphere) in meters.
     */
    double getScaleHeight( ) { return scaleHeight_; }

    //! Set density at zero altitude.
    /*!
     * Sets the density at zero altitude (property of exponential atmosphere) in kg per meter^3.
     */
    void setDensityAtZeroAltitude( const double densityAtZeroAltitude )
    {
        densityAtZeroAltitude_ = densityAtZeroAltitude;
    }

    //! Get density at zero altitude.
    /*!
     * Returns the density at zero altitude (property of exponential atmosphere) in kg per meter^3.
     */
    double getDensityAtZeroAltitude( ) { return densityAtZeroAltitude_; }

    //! Set constant temperature.
    /*!
     * Sets the atmospheric temperature (constant, property of exponential atmosphere) in Kelvin.
     */
    void setConstantTemperature( const double constantTemperature )
    {
        constantTemperature_ = constantTemperature;
    }

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
    void setSpecificGasConstant( const double specificGasConstant )
    {
        specificGasConstant_ = specificGasConstant;
    }

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
    double getDensity( double altitude, double longitude = 0.0,
                       double latitude = 0.0, double time = 0.0 )
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
    double getPressure( double altitude, double longitude = 0.0,
                        double latitude = 0.0, double time = 0.0 )
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
    double getTemperature( double altitude, double longitude = 0.0,
                           double latitude = 0.0, double time = 0.0 )
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

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_EXPONENTIAL_ATMOSPHERE_H
