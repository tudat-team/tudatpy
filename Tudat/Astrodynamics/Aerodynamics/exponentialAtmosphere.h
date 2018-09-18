/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      The accuracy of this model could be increased by implementing different values for the
 *      scale height and temperature for different altitudes (e.g., lower, middle and upper
 *      atmosphere).
 *
 */

#ifndef TUDAT_EXPONENTIAL_ATMOSPHERE_H
#define TUDAT_EXPONENTIAL_ATMOSPHERE_H

#include <memory>

#include <cmath>

#include "Tudat/Basics/utilityMacros.h"

#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Aerodynamics/standardAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{

namespace aerodynamics
{

//! Bodies with predefined exponential atmospheres.
/*!
 *  Bodies with predefined exponential atmospheres.
 */
enum BodiesWithPredefinedExponentialAtmospheres
{
    undefined_body = -1,
    earth = 0,
    mars = 1
};

//! Exponential atmosphere class.
/*!
 * Class with an exponential atmosphere. The user has to initialize this class
 * with a scale altitude, a density at zero meters altitude and a constant
 * temperature. The density is determined by \f$ \rho = \rho_0 e^{-h/H} \f$, with \f$ \rho_0 \f$
 * the density at zero altitude, \f$ h \f$ the altitude, and \f$ H \f$ the scaling
 * altitude. The temperature is taken as constant and the pressure follows from the universal
 * gas law \f$ p = \rho RT \f$.
 */
class ExponentialAtmosphere : public StandardAtmosphere
{
public:

    //! Default constructor.
    /*!
     *  Default constructor setting all parameters manually.
     *  \param scaleHeight Scale height of atmosphere model.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param densityAtZeroAltitude Atmospheric density at zero altitude.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     */
    ExponentialAtmosphere(
            const double scaleHeight,
            const double constantTemperature,
            const double densityAtZeroAltitude,
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4 ):
        scaleHeight_( scaleHeight ),
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
     * Returns the specific gas constant of the atmosphere in J/(kg K), its value is assumed constant,
     * due to the assumption of constant atmospheric composition.
     * \return Specific gas constant in exponential atmosphere.
     */
    double getSpecificGasConstant( ) { return specificGasConstant_; }

    //! Get ratio of specific heats.
    /*!
     * Returns the ratio of specific hears of the atmosphere, its value is assumed constant,
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
     * \param altitude Altitude at which temperature is to be computed (not used since
     * temperature is assumed to be constant).
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
     * Specific gas constant of the atmosphere, its value is assumed constant, due to the assumption of
     * constant atmospheric composition.    
     */
    double specificGasConstant_;

    //! Ratio of specific heats at constant pressure and constant volume.
    /*!
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     *  This value is set to a constant, implying constant atmospheric composition.
     */
    double ratioOfSpecificHeats_;

};

//! Typedef for shared-pointer to ExponentialAtmosphere object.
typedef std::shared_ptr< ExponentialAtmosphere > ExponentialAtmospherePointer;

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_EXPONENTIAL_ATMOSPHERE_H
