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
void TabulatedAtmosphere::initialize( const std::vector< std::string >& atmosphereTableFile )
{
    // Locally store the atmosphere table file name.
    atmosphereTableFile_ = atmosphereTableFile;

    // Retrieve number of independent variables from file.
    numberOfIndependentVariables_ = input_output::getNumberOfIndependentVariablesInCoefficientFile(
                atmosphereTableFile_.at( 0 ) );

    // Check input consistency
    if ( static_cast< int >( independentVariables_.size( ) ) != numberOfIndependentVariables_ )
    {
        throw std::runtime_error( "Error when creating tabulated atmosphere from file, "
                                  "number of specified independent variables differs from file." );
    }

    // Create interpolators
    createAtmosphereInterpolators< 1 >( atmosphereTableFile );
}

//! Initialize atmosphere table reader.
template< int NumberOfIndependentVariables >
void TabulatedAtmosphere::createAtmosphereInterpolators( const std::vector< std::string >& atmosphereTableFile )
{
    using namespace interpolators;

    // Get order of dependent variables
    std::vector< int > dependentVariableIndices = { 0, 0, 0, 0, 0 }; // only 5 dependent variables supported
    for ( unsigned int i = 0; i < dependentVariables_.size( ); i++ )
    {
        if ( i <= dependentVariableIndices.size( ) )
        {
            dependentVariableIndices.at( dependentVariables_.at( i ) ) = i;
            dependentVariablesDependency_.at( dependentVariables_.at( i ) ) = true;
        }
        else
        {
            std::string errorMessage = "Error, dependent variable " + std::to_string( dependentVariables_.at( i ) ) +
                    " not found in tabulated atmosphere.";
            throw std::runtime_error( errorMessage );
        }
    }

    // Check that density, pressure and temperature are present
    if ( !( dependentVariablesDependency_.at( 0 ) || dependentVariablesDependency_.at( 1 ) ||
            dependentVariablesDependency_.at( 2 ) ) )
    {
        throw std::runtime_error( "Error, tabulated atmosphere must be initialized with at least "
                                  "density, pressure and temperature." );
    }

    // Create interpolators for variables requested by users, depending on the number of variables
    switch ( NumberOfIndependentVariables )
    {
    case 1:
    {
        // Create interpolators for density, pressure and temperature
        interpolationForDensity_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                    independentVariablesData_.at( 0 ), tabulatedAtmosphereData.at( dependentVariableIndices.at( 0 ) ) );
        interpolationForPressure_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                    independentVariablesData_.at( 0 ), tabulatedAtmosphereData.at( dependentVariableIndices.at( 1 ) ) );
        interpolationForTemperature_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                    independentVariablesData_.at( 0 ), tabulatedAtmosphereData.at( dependentVariableIndices.at( 2 ) ) );

        // Create remaining interpolators, if requested by user
        if ( dependentVariablesDependency_.at( 3 ) )
        {
            interpolationForGasConstant_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                        independentVariablesData_.at( 0 ), tabulatedAtmosphereData.at( dependentVariableIndices.at( 3 ) ) );
        }
        if ( dependentVariablesDependency_.at( 4 ) )
        {
            interpolationForSpecificHeatRatio_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                        independentVariablesData_.at( 0 ), tabulatedAtmosphereData.at( dependentVariableIndices.at( 4 ) ) );
        }
        break;
    }
    case 2:
    case 3:
    case 4:
    {
        // Call approriate file reading function for N independent variables
        std::pair< std::vector< boost::multi_array< double, static_cast< size_t >( NumberOfIndependentVariables ) > >,
                std::vector< std::vector< double > > > tabulatedAtmosphereData;
        if ( ( NumberOfIndependentVariables > 0 ) && ( NumberOfIndependentVariables < 5 ) )
        {
            // Retrieve number of dependent variables from user.
            int numberOfDependentVariables = dependentVariables_.size( );
            // consistency with number of files is checked in readTabulatedAtmosphere function

            // Extract data
            tabulatedAtmosphereData = input_output::readTabulatedAtmosphere<
                    numberOfDependentVariables, NumberOfIndependentVariables >( atmosphereTableFile );
        }
        else
        {
            throw std::runtime_error( "Error when reading tabulated atmosphere from file, found " +
                                      std::to_string( NumberOfIndependentVariables ) +
                                      " independent variables, up to 4 currently supported." );
        }

        // Assign independent variables
        independentVariablesData_ = tabulatedAtmosphereData.second;

        // Assign dependent variables
        interpolationForDensity_ =
                boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                    independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 0 ) ) );

        interpolationForPressure_ =
                boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                    independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 1 ) ) );
        interpolationForTemperature_ =
                boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                    independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 2 ) ) );
        if ( dependentVariablesDependency_.at( 3 ) )
        {
            interpolationForGasConstant_ =
                    boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                        independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 3 ) ) );
        }
        if ( dependentVariablesDependency_.at( 4 ) )
        {
            interpolationForSpecificHeatRatio_ =
                    boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                        independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 4 ) ) );
        }
        break;
    }
    }
    // Error in case of more than 4 (or less than 1) independent variables has already been given.
}

} // namespace aerodynamics
} // namespace tudat
