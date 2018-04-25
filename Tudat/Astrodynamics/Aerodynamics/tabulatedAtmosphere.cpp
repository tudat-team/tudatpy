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
template< int NumberOfIndependentVariables >
void TabulatedAtmosphere< NumberOfIndependentVariables >::initialize( const std::vector< std::string >& atmosphereTableFile )
{
    // Locally store the atmosphere table file name.
    atmosphereTableFile_ = atmosphereTableFile;

    // Retrieve number of independent variables from file.
    int numberOfIndependentVariablesInFile =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( atmosphereTableFile_.at( 0 ) );

    // Check input consistency
    if ( independentVariables_.size( ) != NumberOfIndependentVariables )
    {
        throw std::runtime_error( "Error when creating tabulated atmosphere from file, "
                                  "number of specified independent variables, differs from file." );
    }

    // Retrieve number of dependent variables from user.
    const int numberOfDependentVariables = dependentVariables_.size( );
    // consistency with number of files is checked in readTabulatedAtmosphere function

    // Call approriate file reading function for N independent variables
    std::pair< std::vector< boost::multi_array< double, static_cast< size_t >( NumberOfIndependentVariables ) > >,
            std::vector< std::vector< double > > > tabulatedAtmosphereData;
    if ( ( NumberOfIndependentVariables > 0 ) && ( NumberOfIndependentVariables < 5 ) )
    {
        tabulatedAtmosphereData =
                input_output::readTabulatedAtmosphere< numberOfDependentVariables, NumberOfIndependentVariables >(
                    atmosphereTableFile_ );
    }
    else
    {
        throw std::runtime_error( "Error when reading tabulated atmosphere from file, found " +
                                  std::to_string( NumberOfIndependentVariables ) +
                                  " independent variables, up to 4 currently supported." );
    }

    // Get order of independent variables
    std::vector< int > independentVariableIndices = { 0, 0, 0, 0 }; // only 4 independent variables supported
    for ( unsigned int i = 0; i < independentVariables_.size( ); i++ )
    {
        if ( i < 5 )
        {
            independentVariableIndices.at( independentVariables_.at( i ) ) = i;
            independentVariablesDependency_.at( independentVariables_.at( i ) ) = true;
        }
        else
        {
            std::string errorMessage = "Error, independent variable " +
                    std::to_string( independentVariables_.at( i ) ) +
                    " not found in tabulated atmosphere.";
            throw std::runtime_error( errorMessage );
        }
    }

    // Get order of dependent variables
    std::vector< int > dependentVariableIndices = { 0, 0, 0, 0, 0 }; // only 5 dependent variables supported
    for ( unsigned int i = 0; i < dependentVariables_.size( ); i++ )
    {
        if ( i < 6 )
        {
            dependentVariableIndices.at( dependentVariables_.at( i ) ) = i;
            dependentVariablesDependency_.at( dependentVariables_.at( i ) ) = true;
        }
        else
        {
            std::string errorMessage = "Error, dependent variable " +
                    std::to_string( dependentVariables_.at( i ) ) +
                    " not found in tabulated atmosphere.";
            throw std::runtime_error( errorMessage );
        }
    }

    // Check that density, pressure and temperature are present
    if ( dependentVariablesDependency_.at( 0 ) || dependentVariablesDependency_.at( 1 ) ||
         dependentVariablesDependency_.at( 2 ) )
    {
        throw std::runtime_error( "Error, tabulated atmosphere must be initialized with at least "
                                  "density, pressure and temperature." );
    }

    // Assign independent variables
    independentVariablesData_ = tabulatedAtmosphereData.second;

    // Assign dependent variables
    if ( dependentVariablesDependency_.at( 0 ) )
    {
        densityData_ = tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 0 ) );
    }
    if ( dependentVariablesDependency_.at( 1 ) )
    {
        pressureData_ = tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 1 ) );
    }
    if ( dependentVariablesDependency_.at( 2 ) )
    {
        temperatureData_ = tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 2 ) );
    }
    if ( dependentVariablesDependency_.at( 3 ) )
    {
        gasConstantData_ = tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 3 ) );
    }
    if ( dependentVariablesDependency_.at( 4 ) )
    {
        specificHeatRatioData_ = tabulatedAtmosphereData.first.at( dependentVariableIndices.at( 4 ) );
    }

    using namespace interpolators;

    // Use cubic spline if only one variable is used
    if ( NumberOfIndependentVariables == 1 )
    {
        // Get only independent variable
        boost::shared_ptr< std::vector< double > > independentVariableDataPointer =
                boost::make_shared< std::vector< double > >( );
        switch ( independentVariables_.at( 0 ) )
        {
        case altitude_dependent_atmosphere:
            independentVariableDataPointer = boost::make_shared< std::vector< double > >(
                        tabulatedAtmosphereData.second.at( independentVariableIndices.at( 0 ) ) );
        case latitude_dependent_atmosphere:
            independentVariableDataPointer = boost::make_shared< std::vector< double > >(
                        tabulatedAtmosphereData.second.at( independentVariableIndices.at( 1 ) ) );
        case longitude_dependent_atmosphere:
            independentVariableDataPointer = boost::make_shared< std::vector< double > >(
                        tabulatedAtmosphereData.second.at( independentVariableIndices.at( 2 ) ) );
        case time_dependent_atmosphere:
            independentVariableDataPointer = boost::make_shared< std::vector< double > >(
                        tabulatedAtmosphereData.second.at( independentVariableIndices.at( 3 ) ) );
        }

        if ( dependentVariablesDependency_.at( 0 ) )
        {
            cubicSplineInterpolationForDensity_ =
                    boost::make_shared< CubicSplineInterpolatorDouble >( independentVariableDataPointer,
                                                                         densityData_ );
        }
        if ( dependentVariablesDependency_.at( 1 ) )
        {
            cubicSplineInterpolationForPressure_ =
                    boost::make_shared< CubicSplineInterpolatorDouble >( independentVariableDataPointer,
                                                                         pressureData_ );
        }
        if ( dependentVariablesDependency_.at( 2 ) )
        {
            cubicSplineInterpolationForTemperature_ =
                    boost::make_shared< CubicSplineInterpolatorDouble >( independentVariableDataPointer,
                                                                         temperatureData_ );
        }
        if ( dependentVariablesDependency_.at( 3 ) )
        {
            cubicSplineInterpolationForGasConstant_ =
                    boost::make_shared< CubicSplineInterpolatorDouble >( independentVariableDataPointer,
                                                                         gasConstantData_);
        }
        if ( dependentVariablesDependency_.at( 4 ) )
        {
            cubicSplineInterpolationForSpecificHeatRatio_ =
                    boost::make_shared< CubicSplineInterpolatorDouble >( independentVariableDataPointer,
                                                                         specificHeatRatioData_);
        }
    }
    else
    {
        if ( dependentVariablesDependency_.at( 0 ) )
        {
            multiLinearInterpolationForDensity_ = boost::make_shared< interpolators::MultiLinearInterpolator
                    < double, double, NumberOfIndependentVariables > >(
                        independentVariablesData_, densityData_ );
        }
        if ( dependentVariablesDependency_.at( 1 ) )
        {
            multiLinearInterpolationForPressure_ = boost::make_shared< interpolators::MultiLinearInterpolator
                    < double, double, NumberOfIndependentVariables > >(
                        independentVariablesData_, pressureData_ );
        }
        if ( dependentVariablesDependency_.at( 2 ) )
        {
            multiLinearInterpolationForTemperature_ = boost::make_shared< interpolators::MultiLinearInterpolator
                    < double, double, NumberOfIndependentVariables > >(
                        independentVariablesData_, temperatureData_ );
        }
        if ( dependentVariablesDependency_.at( 3 ) )
        {
            multiLinearInterpolationForGasConstant_ = boost::make_shared< interpolators::MultiLinearInterpolator
                    < double, double, NumberOfIndependentVariables > >(
                        independentVariablesData_, gasConstantData_ );
        }
        if ( dependentVariablesDependency_.at( 4 ) )
        {
            multiLinearInterpolationForSpecificHeatRatio_ = boost::make_shared< interpolators::MultiLinearInterpolator
                    < double, double, NumberOfIndependentVariables > >(
                        independentVariablesData_, specificHeatRatioData_ );
        }
    }
}

} // namespace aerodynamics
} // namespace tudat
