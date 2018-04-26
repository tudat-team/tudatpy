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
#include "Tudat/Mathematics/Interpolators/multiLinearInterpolator.h"
#include "Tudat/InputOutput/tabulatedAtmosphereReader.h"

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
template< int NumberOfIndependentVariables >
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
    TabulatedAtmosphere( const std::vector< std::string >& atmosphereTableFile,
                         const std::vector< AtmosphereDependentVariables > dependentVariables =
    { density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere },
                         const std::vector< AtmosphereIndependentVariables > independentVariables =
    { altitude_dependent_atmosphere },
                         const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                         const double ratioOfSpecificHeats = 1.4 ):
        atmosphereTableFile_( atmosphereTableFile ), independentVariables_( independentVariables ),
        numberOfIndependentVariables_( NumberOfIndependentVariables ), dependentVariables_( dependentVariables ),
        specificGasConstant_( specificGasConstant ), ratioOfSpecificHeats_( ratioOfSpecificHeats )
    {
        initialize( atmosphereTableFile_ );
    }

    //! Destructor
    ~TabulatedAtmosphere( ){ }

    //! Get atmosphere table file name.
    /*!
     * Returns atmosphere table file name.
     * \return The atmosphere table file.
     */
    std::vector< std::string > getAtmosphereTableFile( ) { return atmosphereTableFile_; }

    //! Get local density.
    /*!
     * Returns the local density parameter of the atmosphere in kg per meter^3.
     * \param altitude Altitude at which density is to be computed.
     * \param longitude Longitude at which density is to be computed.
     * \param latitude Latitude at which density is to be computed.
     * \param time Time at which density is to be computed.
     * \return Atmospheric density at specified conditions.
     */
    double getDensity( const double altitude, const double longitude = 0.0,
                       const double latitude = 0.0, const double time = 0.0 )
    {
        if ( numberOfIndependentVariables_ == 1 )
        {
            // Get only independent variable
            double independentVariableData;
            switch ( independentVariables_.at( 0 ) )
            {
            case altitude_dependent_atmosphere:
                independentVariableData = altitude;
                break;
            case latitude_dependent_atmosphere:
                independentVariableData = latitude;
                break;
            case longitude_dependent_atmosphere:
                independentVariableData = longitude;
                break;
            case time_dependent_atmosphere:
                independentVariableData = time;
                break;
            }

            // Give output
            return cubicSplineInterpolationForDensity_->interpolate( independentVariableData );
        }
        else
        {
            // Get list of independent variables
            std::vector< double > independentVariableData;
            for ( int i = 0; i < numberOfIndependentVariables_; i++ )
            {
                switch ( independentVariables_.at( i ) )
                {
                case altitude_dependent_atmosphere:
                    independentVariableData.push_back( altitude );
                    break;
                case latitude_dependent_atmosphere:
                    independentVariableData.push_back( latitude );
                    break;
                case longitude_dependent_atmosphere:
                    independentVariableData.push_back( longitude );
                    break;
                case time_dependent_atmosphere:
                    independentVariableData.push_back( time );
                    break;
                }
            }

            // Give output
            return multiLinearInterpolationForDensity_->interpolate( independentVariableData );
        }
    }

    //! Get local pressure.
    /*!
     * Returns the local pressure of the atmosphere in Newton per meter^2.
     * \param altitude Altitude  at which pressure is to be computed.
     * \param longitude Longitude at which pressure is to be computed.
     * \param latitude Latitude at which pressure is to be computed.
     * \param time Time at which pressure is to be computed.
     * \return Atmospheric pressure at specified conditions.
     */
    double getPressure( const double altitude, const double longitude = 0.0,
                        const double latitude = 0.0, const double time = 0.0 )
    {
        if ( numberOfIndependentVariables_ == 1 )
        {
            // Get only independent variable
            double independentVariableData;
            switch ( independentVariables_.at( 0 ) )
            {
            case altitude_dependent_atmosphere:
                independentVariableData = altitude;
                break;
            case latitude_dependent_atmosphere:
                independentVariableData = latitude;
                break;
            case longitude_dependent_atmosphere:
                independentVariableData = longitude;
                break;
            case time_dependent_atmosphere:
                independentVariableData = time;
                break;
            }

            // Give output
            return cubicSplineInterpolationForPressure_->interpolate( independentVariableData );
        }
        else
        {
            // Get list of independent variables
            std::vector< double > independentVariableData;
            for ( int i = 0; i < numberOfIndependentVariables_; i++ )
            {
                switch ( independentVariables_.at( i ) )
                {
                case altitude_dependent_atmosphere:
                    independentVariableData.push_back( altitude );
                    break;
                case latitude_dependent_atmosphere:
                    independentVariableData.push_back( latitude );
                    break;
                case longitude_dependent_atmosphere:
                    independentVariableData.push_back( longitude );
                    break;
                case time_dependent_atmosphere:
                    independentVariableData.push_back( time );
                    break;
                }
            }

            // Give output
            return multiLinearInterpolationForPressure_->interpolate( independentVariableData );
        }
    }

    //! Get local temperature.
    /*!
     * Returns the local temperature of the atmosphere in Kelvin.
     * \param altitude Altitude at which temperature is to be computed
     * \param longitude Longitude at which temperature is to be computed.
     * \param latitude Latitude at which temperature is to be computed.
     * \param time Time at which temperature is to be computed.
     * \return constantTemperature Atmospheric temperature at specified conditions.
     */
    double getTemperature( const double altitude, const double longitude = 0.0,
                           const double latitude = 0.0, const double time = 0.0 )
    {
        if ( numberOfIndependentVariables_ == 1 )
        {
            // Get only independent variable
            double independentVariableData;
            switch ( independentVariables_.at( 0 ) )
            {
            case altitude_dependent_atmosphere:
                independentVariableData = altitude;
                break;
            case latitude_dependent_atmosphere:
                independentVariableData = latitude;
                break;
            case longitude_dependent_atmosphere:
                independentVariableData = longitude;
                break;
            case time_dependent_atmosphere:
                independentVariableData = time;
                break;
            }

            // Give output
            return cubicSplineInterpolationForTemperature_->interpolate( independentVariableData );
        }
        else
        {
            // Get list of independent variables
            std::vector< double > independentVariableData;
            for ( int i = 0; i < numberOfIndependentVariables_; i++ )
            {
                switch ( independentVariables_.at( i ) )
                {
                case altitude_dependent_atmosphere:
                    independentVariableData.push_back( altitude );
                    break;
                case latitude_dependent_atmosphere:
                    independentVariableData.push_back( latitude );
                    break;
                case longitude_dependent_atmosphere:
                    independentVariableData.push_back( longitude );
                    break;
                case time_dependent_atmosphere:
                    independentVariableData.push_back( time );
                    break;
                }
            }

            // Give output
            return multiLinearInterpolationForTemperature_->interpolate( independentVariableData );
        }
    }

    //! Get specific gas constant.
    /*!
     * Returns the specific gas constant of the atmosphere in J/(kg K), its value is assumed constant.
     * \param altitude Altitude at which specific gas constant is to be computed.
     * \param longitude Longitude at which specific gas constant is to be computed.
     * \param latitude Latitude at which specific gas constant is to be computed.
     * \param time Time at which specific gas constant is to be computed.
     * \return specificGasConstant Specific gas constant at specified conditions.
     */
    double getSpecificGasConstant( const double altitude, const double longitude = 0.0,
                                   const double latitude = 0.0, const double time = 0.0 )
    {
        if ( dependentVariablesDependency_.at( gas_constant_dependent_atmosphere ) )
        {
            if ( numberOfIndependentVariables_ == 1 )
            {
                // Get only independent variable
                double independentVariableData;
                switch ( independentVariables_.at( 0 ) )
                {
                case altitude_dependent_atmosphere:
                    independentVariableData = altitude;
                    break;
                case latitude_dependent_atmosphere:
                    independentVariableData = latitude;
                    break;
                case longitude_dependent_atmosphere:
                    independentVariableData = longitude;
                    break;
                case time_dependent_atmosphere:
                    independentVariableData = time;
                    break;
                }

                // Give output
                return cubicSplineInterpolationForGasConstant_->interpolate( independentVariableData );
            }
            else
            {
                // Get list of independent variables
                std::vector< double > independentVariableData;
                for ( int i = 0; i < numberOfIndependentVariables_; i++ )
                {
                    switch ( independentVariables_.at( i ) )
                    {
                    case altitude_dependent_atmosphere:
                        independentVariableData.push_back( altitude );
                        break;
                    case latitude_dependent_atmosphere:
                        independentVariableData.push_back( latitude );
                        break;
                    case longitude_dependent_atmosphere:
                        independentVariableData.push_back( longitude );
                        break;
                    case time_dependent_atmosphere:
                        independentVariableData.push_back( time );
                        break;
                    }
                }

                // Give output
                return multiLinearInterpolationForGasConstant_->interpolate( independentVariableData );
            }
        }
        else
        {
            return specificGasConstant_;
        }
    }

    //! Get ratio of specific heats.
    /*!
     * Returns the ratio of specific heats of the air, its value is assumed constant,.
     * \param altitude Altitude at which ratio of specific heats is to be computed
     * \param longitude Longitude at which ratio of specific heats is to be computed.
     * \param latitude Latitude at which ratio of specific heats is to be computed.
     * \param time Time at which ratio of specific heats is to be computed.
     * \return Ratio of specific heats at specified conditions.
     */
    double getRatioOfSpecificHeats( const double altitude, const double longitude = 0.0,
                                    const double latitude = 0.0, const double time = 0.0 )
    {
        if ( dependentVariablesDependency_.at( specific_heat_ratio_dependent_atmosphere ) )
        {
            if ( numberOfIndependentVariables_ == 1 )
            {
                // Get only independent variable
                double independentVariableData;
                switch ( independentVariables_.at( 0 ) )
                {
                case altitude_dependent_atmosphere:
                    independentVariableData = altitude;
                    break;
                case latitude_dependent_atmosphere:
                    independentVariableData = latitude;
                    break;
                case longitude_dependent_atmosphere:
                    independentVariableData = longitude;
                    break;
                case time_dependent_atmosphere:
                    independentVariableData = time;
                    break;
                }

                // Give output
                return cubicSplineInterpolationForSpecificHeatRatio_->interpolate( independentVariableData );
            }
            else
            {
                // Get list of independent variables
                std::vector< double > independentVariableData;
                for ( int i = 0; i < numberOfIndependentVariables_; i++ )
                {
                    switch ( independentVariables_.at( i ) )
                    {
                    case altitude_dependent_atmosphere:
                        independentVariableData.push_back( altitude );
                        break;
                    case latitude_dependent_atmosphere:
                        independentVariableData.push_back( latitude );
                        break;
                    case longitude_dependent_atmosphere:
                        independentVariableData.push_back( longitude );
                        break;
                    case time_dependent_atmosphere:
                        independentVariableData.push_back( time );
                        break;
                    }
                }

                // Give output
                return multiLinearInterpolationForSpecificHeatRatio_->interpolate( independentVariableData );
            }
        }
        else
        {
            return ratioOfSpecificHeats_;
        }
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
     * \return Atmospheric speed of sound at specified conditions.
     */
    double getSpeedOfSound( const double altitude, const double longitude,
                            const double latitude, const double time )
    {
        return computeSpeedOfSound(
                    getTemperature( altitude, longitude, latitude, time ),
                    getSpecificGasConstant( altitude, longitude, latitude, time ),
                    getRatioOfSpecificHeats( altitude, longitude, latitude, time ) );
    }

protected:

private:

    //! Initialize atmosphere table reader.
    /*!
     * Initializes the atmosphere table reader.
     * \param atmosphereTableFile The name of the atmosphere table.
     */
    void initialize( const std::vector< std::string >& atmosphereTableFile );

    //! The file name of the atmosphere table.
    /*!
     *  The file name of the atmosphere table. The file should contain four columns of data,
     *  containing altitude (first column), and the associated density, pressure and density values
     *  in the second, third and fourth columns.
     */
    std::vector< std::string > atmosphereTableFile_;

    //! A vector of strings containing the names of the independent variables contained in the atmosphere file
    /*!
     * A vector of strings containing the names of the independent variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereIndependentVariables > independentVariables_;

    //! Vector containing the independent variables.
    /*!
     *  Vector containing the independent variables.
     */
    std::vector< std::vector< double > > independentVariablesData_;

    //! Integer specifying number of independent variables.
    /*!
     *  Integer specifying number of independent variables.
     */
    int numberOfIndependentVariables_;

    //! Boolean that determines if the atmosphere is dependent on altitude, latitude, longitude and/or time.
    /*!
     *  Boolean that determines if the atmosphere is dependent on altitude, latitude, longitude and/or time.
     *  Note that the order is important.
     */
    std::vector< bool > independentVariablesDependency_ = { false, false, false, false };

    //! A vector of strings containing the names of the variables contained in the atmosphere file
    /*!
     * A vector of strings containing the names of the variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereDependentVariables > dependentVariables_;

    //! Boolean that determines if the atmosphere file contains dentity, pressure, temperature, gas constant and/or
    //! ratio of specific heats.
    /*!
     *  Boolean that determines if the atmosphere file contains dentity, pressure, temperature, gas constant and/or
     *  ratio of specific heats. Note that the order is important.
     */
    std::vector< bool > dependentVariablesDependency_ = { false, false, false, false, false };

    //! Vector containing the density data as a function of the altitude.
    /*!
     *  Vector containing the density data as a function of the altitude.
     */
    boost::multi_array< double, static_cast< size_t >( NumberOfIndependentVariables ) > densityData_;

    //! Vector containing the pressure data as a function of the altitude.
    /*!
     *  Vector containing the pressure data as a function of the altitude.
     */
    boost::multi_array< double, static_cast< size_t >( NumberOfIndependentVariables ) > pressureData_;

    //! Vector containing the temperature data as a function of the altitude.
    /*!
     *  Vector containing the temperature data as a function of the altitude.
     */
    boost::multi_array< double, static_cast< size_t >( NumberOfIndependentVariables ) > temperatureData_;

    //! Vector containing the gas constant data as a function of the altitude.
    /*!
     *  Vector containing the gas constant data as a function of the altitude.
     */
    boost::multi_array< double, static_cast< size_t >( NumberOfIndependentVariables ) > gasConstantData_;

    //! Vector containing the specific heat ratio data as a function of the altitude.
    /*!
     *  Vector containing the specific heat ratio data as a function of the altitude.
     */
    boost::multi_array< double, static_cast< size_t >( NumberOfIndependentVariables ) > specificHeatRatioData_;

    //! Specific gas constant of the atmosphere.
    /*!
     * Specific gas constant of the atmosphere.
     */
    double specificGasConstant_;

    //! Ratio of specific heats of the atmosphere at constant pressure and constant volume.
    /*!
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double ratioOfSpecificHeats_;

    //! Cubic spline interpolation for density.
    /*!
     *  Cubic spline interpolation for density.
     */
    interpolators::CubicSplineInterpolatorDoublePointer cubicSplineInterpolationForDensity_;

    //! N-dimensional linear interpolation for density.
    /*!
     *  N-dimensional linear interpolation for density.
     */
    boost::shared_ptr< interpolators::MultiLinearInterpolator
            < double, double, NumberOfIndependentVariables > > multiLinearInterpolationForDensity_;

    //! Cubic spline interpolation for pressure.
    /*!
     *  Cubic spline interpolation for pressure.
     */
    interpolators::CubicSplineInterpolatorDoublePointer cubicSplineInterpolationForPressure_;

    //! N-dimensional linear interpolation for pressure.
    /*!
     *  N-dimensional linear interpolation for pressure.
     */
    boost::shared_ptr< interpolators::MultiLinearInterpolator
            < double, double, NumberOfIndependentVariables > > multiLinearInterpolationForPressure_;

    //! Cubic spline interpolation for temperature.
    /*!
     *  Cubic spline interpolation for temperature.
     */
    interpolators::CubicSplineInterpolatorDoublePointer cubicSplineInterpolationForTemperature_;

    //! N-dimensional linear interpolation for temperature.
    /*!
     *  N-dimensional linear interpolation for temperature.
     */
    boost::shared_ptr< interpolators::MultiLinearInterpolator
            < double, double, NumberOfIndependentVariables > > multiLinearInterpolationForTemperature_;

    //! Cubic spline interpolation for specific gas constant.
    /*!
     *  Cubic spline interpolation for specific gas constant.
     */
    interpolators::CubicSplineInterpolatorDoublePointer cubicSplineInterpolationForGasConstant_;

    //! N-dimensional linear interpolation for specific gas constant.
    /*!
     *  N-dimensional linear interpolation for specific gas constant.
     */
    boost::shared_ptr< interpolators::MultiLinearInterpolator
            < double, double, NumberOfIndependentVariables > > multiLinearInterpolationForGasConstant_;

    //! Cubic spline interpolation for ratio of specific heats.
    /*!
     *  Cubic spline interpolation for ratio of specific heats.
     */
    interpolators::CubicSplineInterpolatorDoublePointer cubicSplineInterpolationForSpecificHeatRatio_;

    //! N-dimensional linear interpolation for ratio of specific heats.
    /*!
     *  N-dimensional linear interpolation for ratio of specific heats.
     */
    boost::shared_ptr< interpolators::MultiLinearInterpolator
            < double, double, NumberOfIndependentVariables > > multiLinearInterpolationForSpecificHeatRatio_;
};

////! Typedef for shared-pointer to TabulatedAtmosphere object.
//template< int NumberOfIndependentVariables = 1 >
//typedef boost::shared_ptr< TabulatedAtmosphere > TabulatedAtmospherePointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_TABULATED_ATMOSPHERE_H
