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
void TabulatedAtmosphere::initialize( const std::map< int, std::string >& atmosphereTableFile )
{
    // Locally store the atmosphere table file name.
    atmosphereTableFile_ = atmosphereTableFile;

    // Check input consistency
    if ( atmosphereTableFile_.size( ) != dependentVariables_.size( ) )
    {
        throw std::runtime_error( "Error when creating tabulated atmosphere from file, "
                                  "number of specified independent variables differs from file." );
    }

    // Retrieve number of independent variables from file.
    numberOfIndependentVariables_ = input_output::getNumberOfIndependentVariablesInCoefficientFile(
                atmosphereTableFile_.at( 0 ) );

    // Check number of independent variables
    if ( ( numberOfIndependentVariables_ < 1 ) || ( numberOfIndependentVariables_ > 4 ) )
    {
        throw std::runtime_error( "Error when reading tabulated atmosphere from file, found " +
                                  std::to_string( numberOfIndependentVariables_ ) +
                                  " independent variables, up to 4 currently supported." );
    }
    // Could also check to make sure that no duplicate (in)dependent variables are input

    // Check input consistency
    if ( static_cast< int >( independentVariables_.size( ) ) != numberOfIndependentVariables_ )
    {
        throw std::runtime_error( "Error when creating tabulated atmosphere from file, "
                                  "number of specified independent variables differs from file." );
    }

    // Create interpolators
    createAtmosphereInterpolators( atmosphereTableFile );
}

//! Initialize atmosphere table reader.
void TabulatedAtmosphere::createAtmosphereInterpolators( const std::map< int, std::string >& atmosphereTableFile )
{
    using namespace interpolators;

    // Get order of dependent variables
    for ( unsigned int i = 0; i < dependentVariables_.size( ); i++ )
    {
        if ( i <= dependentVariableIndices_.size( ) )
        {
            dependentVariableIndices_.at( dependentVariables_.at( i ) ) = i;
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

    // Retrieve number of dependent variables from user.
    int numberOfDependentVariables = dependentVariables_.size( );
    // consistency with number of files is checked in readTabulatedAtmosphere function

    // Create interpolators for variables requested by users, depending on the number of variables
    switch ( numberOfIndependentVariables_ )
    {
    case 1:
    {
        // Call approriate file reading function for 1 independent variables
        Eigen::MatrixXd tabulatedAtmosphereData = input_output::readMatrixFromFile(
                    atmosphereTableFile_.at( 0 ), " \t", "%" );

        // Check whether data is present in the file.
        if ( tabulatedAtmosphereData.rows( ) < 1 || tabulatedAtmosphereData.cols( ) < 1 )
        {
            std::string errorMessage = "The atmosphere table file " + atmosphereTableFile_.at( 0 ) + " is empty";
            throw std::runtime_error( errorMessage );
        }

        // Check
        if ( numberOfDependentVariables != ( tabulatedAtmosphereData.cols( ) - 1 ) )
        {
            throw std::runtime_error( "Number of specified dependent variables does not match file." );
        }

        // Extract variables from file
        std::vector< std::vector< double > > dependentVariablesData;
        for ( unsigned int i = 0; i < tabulatedAtmosphereData.rows( ); i++ )
        {
            independentVariablesData_.at( 0 ).at( i ) = tabulatedAtmosphereData( i, 0 );
            for ( int j = 0; j < numberOfDependentVariables; j++ )
            {
                if ( dependentVariablesDependency_.at( j ) )
                {
                    dependentVariablesData.at( j ).at( i ) = tabulatedAtmosphereData( i, dependentVariableIndices_.at( j ) );
                }
            }
        }

        // Create interpolators for density, pressure and temperature
        interpolationForDensity_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                    independentVariablesData_.at( 0 ), dependentVariablesData.at( 0 ) );
        interpolationForPressure_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                    independentVariablesData_.at( 0 ), dependentVariablesData.at( 1 ) );
        interpolationForTemperature_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                    independentVariablesData_.at( 0 ), dependentVariablesData.at( 2 ) );

        // Create remaining interpolators, if requested by user
        if ( dependentVariablesDependency_.at( 3 ) )
        {
            interpolationForGasConstant_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                        independentVariablesData_.at( 0 ), dependentVariablesData.at( 3 ) );
        }
        if ( dependentVariablesDependency_.at( 4 ) )
        {
            interpolationForSpecificHeatRatio_ = boost::make_shared< CubicSplineInterpolatorDouble >(
                        independentVariablesData_.at( 0 ), dependentVariablesData.at( 4 ) );
        }
        break;
    }
    case 2:
    {
        createMultiDimensionalAtmosphereInterpolators< 2 >( );
    }
    case 3:
    {
        createMultiDimensionalAtmosphereInterpolators< 3 >( );
    }
    case 4:
    {
        throw std::runtime_error( "Currently, only three independent variables are supported." );
    }
    }
    // Error in case of more than 4 (or less than 1) independent variables has already been given.
}

template< int NumberOfIndependentVariables >
void TabulatedAtmosphere::createMultiDimensionalAtmosphereInterpolators( )
{
    using namespace interpolators;

    // Call approriate file reading function for N independent variables
    std::pair< std::vector< boost::multi_array< double, static_cast< size_t >( NumberOfIndependentVariables ) > >,
            std::vector< std::vector< double > > > tabulatedAtmosphereData;

    // Extract data
    tabulatedAtmosphereData = input_output::readTabulatedAtmosphere< NumberOfIndependentVariables >( atmosphereTableFile_ );

    // Assign independent variables
    independentVariablesData_ = tabulatedAtmosphereData.second;

    // Assign dependent variables
    interpolationForDensity_ =
            boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 0 ) ) );
    interpolationForPressure_ =
            boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 1 ) ) );
    interpolationForTemperature_ =
            boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 2 ) ) );
    if ( dependentVariablesDependency_.at( 3 ) )
    {
        interpolationForGasConstant_ =
                boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                    independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 3 ) ) );
    }
    if ( dependentVariablesDependency_.at( 4 ) )
    {
        interpolationForSpecificHeatRatio_ =
                boost::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                    independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 4 ) ) );
    }
}

} // namespace aerodynamics
} // namespace tudat
