/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Justus, C., Duvall, A., and Keller, V., “Atmospheric Models for Aerocapture,” in
 *          40th AIAA/ASME/SAE/ASEE Joint propulsion Conference, Fort Lauderdale, Florida,
 *          United States, July 2004.
 *      Jah, M., Lisano, M., Born, G., and Axelrad, P., “Mars Aerobraking Spacecraft State
 *          estimation By Processing Inertial Measurement Unit Data,” Journal of Guidance,
 *          Control, and Dynamics, vol. 31, no. 6, pp. 1802–1812, November–December 2008.
 *
 */

#ifndef TUDAT_CUSTOM_CONSTANT_TEMPERATURE_ATMOSPHERE_H
#define TUDAT_CUSTOM_CONSTANT_TEMPERATURE_ATMOSPHERE_H

#include <boost/bind/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <cmath>

#include "tudat/basics/utilityMacros.h"

#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/astro/aerodynamics/standardAtmosphere.h"
#include "tudat/astro/basic_astro/physicalConstants.h"

using namespace boost::placeholders;

namespace tudat
{

namespace aerodynamics
{

enum AvailableConstantTemperatureAtmosphereModels
{
    exponential_atmosphere_model = 0,
    three_wave_atmosphere_model = 1,
    three_term_atmosphere_model = 2
};

//! First atmosphere model, based on exponential atmosphere.
/*!
 *  First atmosphere model, based on exponential atmosphere.
 *  \param altitude Current altitude.
 *  \param longitude Current longitude (unused).
 *  \param latitude Current latitude (unused).
 *  \param time Current time (unused).
 *  \param referenceAltitude Reference altitude.
 *  \param densityAtReferenceAltitude Density at reference altitude condition.
 *  \param scaleHeight Scale height of the atmosphere.
 *  \return Density at current conditions, based on the exponential model with given characteristics.
 */
double exponentialAtmosphereModel( const double altitude, const double longitude, const double latitude, const double time,
                                   const double referenceAltitude, const double densityAtReferenceAltitude, const double scaleHeight );

//! Second atmosphere model, based on a three longitudinal waves model.
/*!
 *  Second atmosphere model, based on a three longitudinal waves model. The parameters for the three
 *  longitudinal waves are hard coded in the function. Based on (Justus, et al., 2004).
 *  \param altitude Current altitude.
 *  \param longitude Current longitude.
 *  \param latitude Current latitude (unused).
 *  \param time Current time (unused).
 *  \param referenceAltitude Reference altitude.
 *  \param densityAtReferenceAltitude Density at reference altitude condition.
 *  \param scaleHeight Scale height of the atmosphere.
 *  \param uncertaintyFactor Factor representing uncertainty in the atmosphere. Can be a
 *      random variable.
 *  \param dustStormFactor Factor representing presence of planet-wide dust storm (particularly
 *      useful for Mars applications).
 *  \return Density at current conditions, based on the three-waves model with given characteristics.
 */
double threeWaveAtmosphereModel( const double altitude, const double longitude, const double latitude, const double time,
                                 const double referenceAltitude, const double densityAtReferenceAltitude, const double scaleHeight,
                                 const double uncertaintyFactor, const double dustStormFactor );

//! Third atmosphere model, based on three constant scale height atmospheres.
/*!
 *  Third atmosphere model, based on three constant scale height atmospheres. The first term is
 *  the same as for the exponential model, whereas the other two are a sine and a cosine term.
 *  Based on (Jah, et al., 2008).
 *  \param altitude Current altitude.
 *  \param longitude Current longitude (unused).
 *  \param latitude Current latitude (unused).
 *  \param time Current time (unused).
 *  \param referenceAltitude Reference altitude.
 *  \param densityAtReferenceAltitude Density at reference altitude condition.
 *  \param scaleHeight Scale height of the atmosphere.
 *  \param modelWeights Weights for each of the three models.
 *  \return Density at current conditions, based on the exponential model with given characteristics.
 */
double threeTermAtmosphereModel( const double altitude, const double longitude, const double latitude, const double time,
                                 const double referenceAltitude, const double densityAtReferenceAltitude, const double scaleHeight,
                                 const std::vector< double >& modelWeights );

//! Custom constant temperature atmosphere class.
/*!
 *  Custom constant temperature atmosphere class. This atmosphere model is initialized
 *  by giving as input a function for the density which may depend on altitude, longitude,
 *  latitude and time. The pressure is then determined with the ideal gas law, and the
 *  temperature is assumed to be constant throughtout the atmosphere domain.
 */
class CustomConstantTemperatureAtmosphere : public AtmosphereModel
{
public:

    //! Typedef for density function.
    typedef std::function< double( const double, const double,
                                   const double, const double ) > DensityFunction;

    //! Default constructor.
    /*!
     *  Default constructor setting all parameters manually.
     *  \param densityFunction Function to retireve the density at the current altitude,
     *      longitude, latitude and time.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     */
    CustomConstantTemperatureAtmosphere(
            const DensityFunction& densityFunction,
            const double constantTemperature,
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4 ) :
        densityFunction_( densityFunction ),
        constantTemperature_( constantTemperature ),
        specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats )
    { }

    //! Constructor.
    /*!
     *  Constructor which uses one of the built-in density functions as input.
     *  \param densityFunctionType Enumeration denoting which density function to implement.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param modelSpecificParameters Vector of parameters to be used to set up the density
     *      function. Both meaning and number of parameters depends on the model.
     */
    CustomConstantTemperatureAtmosphere(
            const AvailableConstantTemperatureAtmosphereModels densityFunctionType,
            const double constantTemperature,
            const double specificGasConstant,
            const double ratioOfSpecificHeats,
            const std::vector< double >& modelSpecificParameters );

    //! Get the function to compute the density at the current conditions.
    /*!
     *  Function to return the function to compute the density at the current conditions.
     *  \return Function to compute the density at the current conditions.
     */
    DensityFunction getDensityFunction( ) { return densityFunction_; }

    //! Set the function to compute the density at the current conditions.
    /*!
     *  Function to reset the function to compute the density at the current conditions.
     *  \param newDensityFunction New function to compute the density at the current conditions.
     */
    void setDensityFunction( DensityFunction& newDensityFunction )
    {
        densityFunction_ = newDensityFunction;
    }

    //! Get constant temperature.
    /*!
     *  Returns the atmospheric temperature (constant, property of exponential atmosphere) in
     *  Kelvin.
     *  \return Constant atmospheric temperature in exponential atmosphere.
     */
    double getConstantTemperature( ) { return constantTemperature_; }

    //! Get specific gas constant.
    /*!
     *  Returns the specific gas constant of the atmosphere in J/(kg K), its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Specific gas constant in exponential atmosphere.
     */
    double getSpecificGasConstant( ) { return specificGasConstant_; }

    //! Get ratio of specific heats.
    /*!
     *  Returns the ratio of specific hears of the atmosphere, its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Ratio of specific heats exponential atmosphere.
     */
    double getRatioOfSpecificHeats( ) { return ratioOfSpecificHeats_; }

    //! Get local density.
    /*!
     *  Returns the local density of the atmosphere in kg per meter^3.
     *  \param altitude Altitude at which density is to be computed.
     *  \param longitude Longitude at which density is to be computed.
     *  \param latitude Latitude at which density is to be computed.
     *  \param time Time at which density is to be computed.
     *  \return Atmospheric density at specified altitude.
     */
    double getDensity( const double altitude, const double longitude = 0.0,
                       const double latitude = 0.0, const double time = 0.0 )
    {
        return densityFunction_( altitude, longitude, latitude, time );
    }

    //! Get local pressure.
    /*!
     *  Returns the local pressure of the atmosphere in Newton per meter^2.
     *  \param altitude Altitude  at which pressure is to be computed.
     *  \param longitude Longitude at which pressure is to be computed.
     *  \param latitude Latitude at which pressure is to be computed.
     *  \param time Time at which pressure is to be computed.
     *  \return Atmospheric pressure at specified altitude.
     */
    double getPressure( const double altitude, const double longitude = 0.0,
                        const double latitude = 0.0, const double time = 0.0 )
    {
        return getDensity( altitude, longitude, latitude, time ) * specificGasConstant_ * constantTemperature_;
    }

    //! Get local temperature.
    /*!
     *  Returns the local temperature of the atmosphere in Kelvin.
     *  \param altitude Altitude at which temperature is to be computed (not used since
     *      temperature is assumed to be constant).
     *  \param longitude Longitude at which temperature is to be computed (not used but included for
     *      consistency with base class interface).
     *  \param latitude Latitude at which temperature is to be computed (not used but included for
     *      consistency with base class interface).
     *  \param time Time at which temperature is to be computed (not used but included for
     *      consistency with base class interface).
     *  \return constantTemperature Atmospheric temperature at specified altitude.
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
     *  Returns the speed of sound in the atmosphere in m/s.
     *  \param altitude Altitude at which speed of sounds is to be computed.
     *  \param longitude Longitude at which speed of sounds is to be computed (not used but included
     *      for consistency with base class interface).
     *  \param latitude Latitude at which speed of sounds is to be computed (not used but included
     *      for consistency with base class interface).
     *  \param time Time at which speed of sounds is to be computed (not used but included for
     *      consistency with base class interface).
     *  \return Atmospheric speed of sounds at specified altitude.
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

    //! Function to compute the density at the current conditions.
    /*!
     *  Function to compute the density at the current conditions. Note that the independent
     *  variables need to be altitude, longitude, latitude and time, in this precise order.
     */
    DensityFunction densityFunction_;

    //! Constant temperature.
    /*!
     *  The atmospheric temperature (constant, property of exponential atmosphere) in Kelvin.
     */
    const double constantTemperature_;

    //! Specific gas constant.
    /*!
     *  Specific gas constant of the atmosphere, its value is assumed constant, due to the
     *  assumption of constant atmospheric composition.
     */
    const double specificGasConstant_;

    //! Ratio of specific heats at constant pressure and constant volume.
    /*!
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     *  This value is set to a constant, implying constant atmospheric composition.
     */
    const double ratioOfSpecificHeats_;

};

//! Typedef for shared-pointer to CustomConstantTemperatureAtmosphere object.
typedef std::shared_ptr< CustomConstantTemperatureAtmosphere > CustomConstantTemperatureAtmospherePointer;

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_CUSTOM_CONSTANT_TEMPERATURE_ATMOSPHERE_H
