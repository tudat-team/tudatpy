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
 *      The provided USSA1976 table file, generated with the pascal file, has a small error which
 *      can be observed at the pressure at sea level. This in 101320 in the file but should be
 *      101325. If this error is not acceptable, another table file should be used.
 *
 */

#ifndef TUDAT_TABULATED_ATMOSPHERE_H
#define TUDAT_TABULATED_ATMOSPHERE_H

#include <string>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/utilityMacros.h"

#include "Tudat/Astrodynamics/Aerodynamics/standardAtmosphere.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"

namespace tudat
{
namespace aerodynamics
{

//! Tabulated atmosphere class.
/*!
 * Tabulated atmospheres class, for example US1976. The default path from which the files are
 * obtained is: /External/AtmosphereTables
 * NOTE: for the moment it only works for tables with 4 columns: altitude, density, pressure and
 * temperature.
 */
class TabulatedAtmosphere : public StandardAtmosphere
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param atmosphereTableFile File containing atmospheric properties.
     *  The file name of the atmosphere table. The file should contain four columns of data,
     *  containing altitude (first column), and the associated density, pressure and density values
     *  in the second, third and fourth columns
     *  \param specificGasConstant The constant specific gas constant of the air
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the air
     */
    TabulatedAtmosphere( const std::string& atmosphereTableFile,
                         const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                         const double ratioOfSpecificHeats = 1.4 )
        : atmosphereTableFile_( atmosphereTableFile ), specificGasConstant_( specificGasConstant ),
          ratioOfSpecificHeats_( ratioOfSpecificHeats )
    {
        initialize( atmosphereTableFile_ );
    }

    //! Get atmosphere table file name.
    /*!
     * Returns atmosphere table file name.
     * \return The atmosphere table file.
     */
    std::string getAtmosphereTableFile( ) { return atmosphereTableFile_; }

    //! Get specific gas constant.
    /*!
     * Returns the specific gas constant of the air in J/(kg K), its value is assumed constant.
     * \return specificGasConstant Specific gas constant in exponential atmosphere.
     */
    double getSpecificGasConstant( ) { return specificGasConstant_; }

    //! Get ratio of specific heats.
    /*!
     * Returns the ratio of specific hears of the air, its value is assumed constant,.
     * \return Ratio of specific heats exponential atmosphere.
     */
    double getRatioOfSpecificHeats( ) { return ratioOfSpecificHeats_; }

    //! Get local density.
    /*!
     * Returns the local density parameter of the atmosphere in kg per meter^3.
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
        return cubicSplineInterpolationForDensity_->interpolate( altitude );
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
        return cubicSplineInterpolationForPressure_->interpolate( altitude );
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
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return cubicSplineInterpolationForTemperature_->interpolate( altitude );
    }

    //! Get local speed of sound in the atmosphere.
    /*!
     * Returns the speed of sound in the atmosphere in m/s.
     * \param altitude Altitude at which speed of sound is to be computed.
     * \param longitude Longitude at which speed of sound is to be computed (not used but included
     * for consistency with base class interface).
     * \param latitude Latitude at which speed of sound is to be computed (not used but included
     * for consistency with base class interface).
     * \param time Time at which speed of sound is to be computed (not used but included for
     * consistency with base class interface).
     * \return Atmospheric speed of sound at specified altitude.
     */
    double getSpeedOfSound( const double altitude, const double longitude = 0.0,
                            const double latitude = 0.0, const double time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return computeSpeedOfSound(
                    getTemperature( altitude, longitude, latitude, time ), ratioOfSpecificHeats_,
                    specificGasConstant_ );
    }

protected:

private:

    //! Initialize atmosphere table reader.
    /*!
     * Initializes the atmosphere table reader.
     * \param atmosphereTableFile The name of the atmosphere table.
     */
    void initialize( const std::string& atmosphereTableFile );

    //! The file name of the atmosphere table.
    /*!
     *  The file name of the atmosphere table. The file should contain four columns of data,
     *  containing altitude (first column), and the associated density, pressure and density values
     *  in the second, third and fourth columns.
     */
    std::string atmosphereTableFile_;

    //! Vector containing the altitude.
    /*!
     *  Vector containing the altitude.
     */
    std::vector< double > altitudeData_;

    //! Vector containing the density data as a function of the altitude.
    /*!
     *  Vector containing the density data as a function of the altitude.
     */
    std::vector< double > densityData_;

    //! Vector containing the pressure data as a function of the altitude.
    /*!
     *  Vector containing the pressure data as a function of the altitude.
     */
    std::vector< double > pressureData_;

    //! Vector containing the temperature data as a function of the altitude.
    /*!
     *  Vector containing the temperature data as a function of the altitude.
     */
    std::vector< double > temperatureData_;

    //! Cubic spline interpolation for density.
    /*!
     *  Cubic spline interpolation for density.
     */
    interpolators::CubicSplineInterpolatorDoublePointer cubicSplineInterpolationForDensity_;

    //! Cubic spline interpolation for pressure.
    /*!
     *  Cubic spline interpolation for pressure.
     */
    interpolators::CubicSplineInterpolatorDoublePointer cubicSplineInterpolationForPressure_;

    //! Cubic spline interpolation for temperature.
    /*!
     *  Cubic spline interpolation for temperature.
     */
    interpolators::CubicSplineInterpolatorDoublePointer cubicSplineInterpolationForTemperature_;

    //! Specific gas constant.
    /*!
     * Specific gas constant of the air, its value is assumed constant, due to the assumption of
     * constant atmospheric composition.
     */
    double specificGasConstant_;

    /*!
     *  Ratio of specific heats of the atmosphrer at constant pressure and constant volume.
     *  This value is set to a constant, implying constant atmospheric composition.
     */
    double ratioOfSpecificHeats_;
};

//! Typedef for shared-pointer to TabulatedAtmosphere object.
typedef boost::shared_ptr< TabulatedAtmosphere > TabulatedAtmospherePointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_TABULATED_ATMOSPHERE_H
