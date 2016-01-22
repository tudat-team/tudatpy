/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      130120    K. Kumar          Made function calls const-correct; added shared-pointer
 *                                  typedef.
 *      140129    D. Dirkx          Changed Doxygen descriptions
 *      140130    T. Roegiers       Changed Doxygen descriptions
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_EXPONENTIAL_ATMOSPHERE_H
#define TUDAT_EXPONENTIAL_ATMOSPHERE_H

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include <cmath>

#include "Tudat/Basics/utilityMacros.h"

#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

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
     *  Default constructor setting all parameters manually.
     *  \param scaleHeight Scale height of atmosphere model.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param densityAtZeroAltitude Atmospheric density at zero altitude.
     *  \param specificGasConstant The constant specific gas constant of the air
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the air
     */
    ExponentialAtmosphere(
            const double scaleHeight,
            const double constantTemperature,
            const double densityAtZeroAltitude,
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4 )
        : scaleHeight_( scaleHeight ),
          constantTemperature_( constantTemperature ),
          densityAtZeroAltitude_( densityAtZeroAltitude ),
          specificGasConstant_( specificGasConstant ),
          ratioOfSpecificHeats_( ratioOfSpecificHeats )
    { }

    //! Constructor from default atmospheric settings.
    /*!
     *  Constructor from default atmospheric settings.
     *  \param bodyWithPredefinedExponentialAtmosphere Identifier of body for which the
     *  atmosphere is to be created.
     */
    ExponentialAtmosphere(
         const BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere );

    //! Get scale height.
    /*!
     * Returns the scale height (property of exponential atmosphere) in meters.
     * \return scaleHeight Scale height of exponential atmosphere.
     */
    double getScaleHeight( ) { return scaleHeight_; }

    //! Get density at zero altitude.
    /*!
     * Returns the density at zero altitude (property of exponential atmosphere) in kg per meter^3.
     * \return densityAtZeroAltitude Atmospheric density at zero altitude.
     */
    double getDensityAtZeroAltitude( ) { return densityAtZeroAltitude_; }

    //! Get constant temperature.
    /*!
     * Returns the atmospheric temperature (constant, property of exponential atmosphere) in
     * Kelvin.
     * \return constantTemperature Constant atmospheric temperature in exponential atmosphere.
     */
    double getConstantTemperature( ) { return constantTemperature_; }

    //! Get specific gas constant.
    /*!
     * Returns the specific gas constant of the air in J/(kg K), its value is assumed constant,
     * due to the assumption of constant atmospheric composition.
     * \return Specific gas constant in exponential atmosphere.
     */
    double getSpecificGasConstant( ) { return specificGasConstant_; }

    //! Get ratio of specific heats.
    /*!
     * Returns the ratio of specific hears of the air, its value is assumed constant,
     * due to the assumption of constant atmospheric composition.
     * \return Ratio of specific heats exponential atmosphere.
     */
    double getRatioOfSpecificHeats( ) { return ratioOfSpecificHeats_; }

    //! Get local density.
    /*!
     * Returns the local density of the atmosphere in kg per meter^3.
     * \param altitude Altitude at which density is to be computed.
     * \param longitude Longitude at which density is to be computed (not used but included for
     * consistency with base class interface).
     * \param latitude Latitude at which density is to be computed (not used but included for
     * consistency with base class interface).
     * \param time Time at which density is to be computed (not used but included for
     * consistency with base class interface).
     * \return Atmospheric density at specified altitude.
     */
    double getDensity( const double altitude, const double longitude = 0.0,
                       const double latitude = 0.0, const double time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return densityAtZeroAltitude_ * std::exp( -altitude / scaleHeight_ );
    }

    //! Get local pressure.
    /*!
     * Returns the local pressure of the atmosphere in Newton per meter^2.
     * \param altitude Altitude  at which pressure is to be computed.
     * \param longitude Longitude at which pressure is to be computed (not used but included for
     * consistency with base class interface).
     * \param latitude Latitude at which pressure is to be computed (not used but included for
     * consistency with base class interface).
     * \param time Time at which pressure is to be computed (not used but included for
     * consistency with base class interface).
     * \return Atmospheric pressure at specified altitude.
     */
    double getPressure( const double altitude, const double longitude = 0.0,
                        const double latitude = 0.0, const double time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return getDensity( altitude ) * specificGasConstant_ * constantTemperature_;
    }

    //! Get local temperature.
    /*!
     * Returns the local temperature of the atmosphere in Kelvin.
     * \param altitude Altitude at which temperature is to be computed
     * \param longitude Longitude at which temperature is to be computed (not used but included for
     * consistency with base class interface).
     * \param latitude Latitude at which temperature is to be computed (not used but included for
     * consistency with base class interface).
     * \param time Time at which temperature is to be computed (not used but included for
     * consistency with base class interface).
     * \return constantTemperature Atmospheric temperature at specified altitude.
     */
    double getTemperature( const double altitude, const double longitude = 0.0,
                           const double latitude = 0.0, const double time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( altitude );
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return constantTemperature_;
    }

    //! Get local speed of sound in the atmosphere.
    /*!
     * Returns the speed of sound in the atmosphere in m/s.
     * \param altitude Altitude at which speed of sounds is to be computed.
     * \param longitude Longitude at which speed of sounds is to be computed (not used but included
     * for consistency with base class interface).
     * \param latitude Latitude at which speed of sounds is to be computed (not used but included
     * for consistency with base class interface).
     * \param time Time at which speed of sounds is to be computed (not used but included for
     * consistency with base class interface).
     * \return Atmospheric speed of sounds at specified altitude.
     */
    double getSpeedOfSound( const double altitude, const double longitude = 0.0,
                                    const double latitude = 0.0, const double time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( altitude );
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return computeSpeedOfSound(
                    getTemperature( altitude, longitude, latitude, time ), ratioOfSpecificHeats_,
                    specificGasConstant_ );
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

    //! Ratio of specific heats at constant pressure and constant volume.
    /*!
     *  Ratio of specific heats of the atmosphrer at constant pressure and constant volume.
     *  This value is set to a constant, implying constant atmospheric composition.
     */
    double ratioOfSpecificHeats_;

};

//! Typedef for shared-pointer to ExponentialAtmosphere object.
typedef boost::shared_ptr< ExponentialAtmosphere > ExponentialAtmospherePointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_EXPONENTIAL_ATMOSPHERE_H
