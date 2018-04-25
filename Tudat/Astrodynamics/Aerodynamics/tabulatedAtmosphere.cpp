/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <boost/make_shared.hpp>
#include <iostream>
#include "Tudat/InputOutput/matrixTextFileReader.h"

#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"

namespace tudat
{
namespace aerodynamics
{

//! Initialize atmosphere table reader.
void TabulatedAtmosphere::initialize( const std::vector< std::string >& atmosphereTableFile,
                                      const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
                                      independentVariableNames )
{
    // Locally store the atmosphere table file name.
    atmosphereTableFile_ = atmosphereTableFile;

    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( atmosphereTableFile_.begin( )->second );

    // Check input consistency
    if( independentVariableNames.size( ) != numberOfIndependentVariables )
    {
        throw std::runtime_error( "Error when creating tabulated atmosphere from file, input sizes are inconsistent" );
    }

    // Call approriate file reading function for N independent variables
    std::pair< boost::multi_array< Eigen::Vector1d, numberOfIndependentVariables >,
            std::vector< std::vector< double > > > containerOfAtmosphereTableFileData;
    if ( ( numberOfIndependentVariables > 0 ) && ( numberOfIndependentVariables < 4 ) )
    {
        containerOfAtmosphereTableFileData =
                input_output::readTabulatedAtmosphere< numberOfIndependentVariables >( atmosphereTableFile_ );
    }
    else
    {
        throw std::runtime_error( "Error when reading tabulated atmosphere from file, found " +
                                  std::to_string( numberOfIndependentVariables ) +
                                  " independent variables, up to 3 currently supported" );
    }

    // Initialize vectors.
    altitudeData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    densityData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    pressureData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    temperatureData_.resize( containerOfAtmosphereTableFileData.rows( ) );

    // Get order of dependent variables
    int densityIndex = 0;
    int pressureIndex = 0;
    int temperatureIndex = 0;
    int specificHeatRatioIndex = 0;
    int gasConstantIndex = 0;
    containsSpecificHeatRatio_ = false;
    containsGasConstant_ = false;
    for( unsigned int i = 0; i < dependentVariables_.size( ); i++ )
    {
        switch( dependentVariables_[ i ] )
        {
        case density_dependent_atmosphere:
            densityIndex = i + 1;
            break;
        case pressure_dependent_atmosphere:
            pressureIndex = i + 1;
            break;
        case temperature_dependent_atmosphere:
            temperatureIndex = i + 1;
            break;
        case specific_heat_ratio_dependent_atmosphere:
            containsSpecificHeatRatio_ = true;
            specificHeatRatioIndex = i + 1;
            break;
        case gas_constant_dependent_atmosphere:
            containsGasConstant_ = true;
            gasConstantIndex = i + 1;
            break;
        default:
            std::string errorMessage = "Error, dependent variable " +
                    std::to_string( dependentVariables_[i] ) +
                    " not found in tabulated atmosphere";
            throw std::runtime_error( errorMessage );
        }
    }

    if( densityIndex == 0 || pressureIndex == 0 || temperatureIndex == 0 )
    {
        throw std::runtime_error(
                    "Error, tabulated atmosphere must be initialized with at least temperature, pressure and density" );
    }

    // Loop through all the strings stored in the container and store the data
    // in the right Eigen::VectorXd.
    for ( int i = 0; i < containerOfAtmosphereTableFileData.rows( ); i++  )
    {
        altitudeData_[ i ] = containerOfAtmosphereTableFileData( i, 0 );
        densityData_[ i ] = containerOfAtmosphereTableFileData( i, densityIndex );
        pressureData_[ i ] = containerOfAtmosphereTableFileData( i, pressureIndex );
        temperatureData_[ i ] = containerOfAtmosphereTableFileData( i, temperatureIndex );
        if( containsSpecificHeatRatio_ )
        {
            specificHeatRatioData_[ i ] = containerOfAtmosphereTableFileData( i, specificHeatRatioIndex );
        }
        if( containsGasConstant_ )
        {
            gasConstantData_[ i ] = containerOfAtmosphereTableFileData( i, gasConstantIndex );
        }
    }


    using namespace interpolators;

    cubicSplineInterpolationForDensity_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, densityData_ );
    cubicSplineInterpolationForPressure_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, pressureData_ );
    cubicSplineInterpolationForTemperature_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, temperatureData_ );

    if( containsSpecificHeatRatio_ )
    {
        cubicSplineInterpolationForSpecificHeatRatio_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, specificHeatRatioData_);
    }

    if( containsGasConstant_ )
    {
        cubicSplineInterpolationForGasConstant_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, gasConstantData_);
    }
}

} // namespace aerodynamics
} // namespace tudat
