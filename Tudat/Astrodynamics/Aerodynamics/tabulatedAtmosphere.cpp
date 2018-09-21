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

#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"

#include <iostream>
#include <boost/make_shared.hpp>
#include "Tudat/InputOutput/matrixTextFileReader.h"

namespace tudat
{

namespace aerodynamics
{

using namespace interpolators;

//! Check uniqueness of input.
template< typename VariableType >
void checkVariableUniqueness( std::vector< VariableType > variables )
{
    // Sort variables
    std::sort( variables.begin( ), variables.end( ) );

    // Check uniqueness
    unsigned int numberOfUniqueElements = std::distance( variables.begin( ),
                                                         std::unique( variables.begin( ), variables.end( ) ) );

    // Give error in case of non-unique variables
    if ( numberOfUniqueElements != variables.size( ) )
    {
        throw std::runtime_error( "Error, in tabulated atmosphere. Duplicate entry in (in)dependent variables." );
    }
}

//! Function to create the interpolators based on the tabulated atmosphere files.
void TabulatedAtmosphere::createAtmosphereInterpolators( )
{
    // Check uniqueness
    checkVariableUniqueness< AtmosphereDependentVariables >( dependentVariables_ );
    checkVariableUniqueness< AtmosphereIndependentVariables >( independentVariables_ );

    // Retrieve number of dependent variables from user.
    unsigned int numberOfDependentVariables = dependentVariables_.size( );

    // Check that number of dependent variables does not exceed limit
    if ( numberOfDependentVariables > dependentVariablesDependency_.size( ) )
    {
        throw std::runtime_error( "Error, in tabulated atmosphere. Number of dependent variables exceeds current limit." );
    }

    // Check input consistency
    if ( independentVariables_.size( ) != 1 )
    {
        if ( atmosphereTableFile_.size( ) != numberOfDependentVariables )
        {
            throw std::runtime_error( "Error when creating tabulated atmosphere from file, "
                                      "number of specified dependent variables differs from file." );
        }

        // Retrieve number of independent variables from file.
        numberOfIndependentVariables_ = input_output::getNumberOfIndependentVariablesInCoefficientFile(
                    atmosphereTableFile_.at( 0 ) );

        // Check that number of independent variables does not exceed limit
        if ( ( numberOfIndependentVariables_ < 1 ) || ( numberOfIndependentVariables_ > 4 ) )
        {
            throw std::runtime_error( "Error when reading tabulated atmosphere from file, found " +
                                      std::to_string( numberOfIndependentVariables_ ) +
                                      " independent variables, up to 4 currently supported." );
        }

        // Check input consistency
        if ( independentVariables_.size( ) != numberOfIndependentVariables_ )
        {
            throw std::runtime_error( "Error when creating tabulated atmosphere from file, "
                                      "number of specified independent variables differs from file." );
        }
    }
    else
    {
        numberOfIndependentVariables_ = 1; // if only one independent variable is specified, only one file will
                                           // be provided, and it cannot be opened with the same function
    }

    // Get order of dependent variables
    for ( unsigned int i = 0; i < numberOfDependentVariables; i++ )
    {
        dependentVariablesDependency_.at( dependentVariables_.at( i ) ) = true;
        dependentVariableIndices_.at( dependentVariables_.at( i ) ) = i;
    }

    // Check that density, pressure and temperature are present
    if ( !( dependentVariablesDependency_.at( 0 ) || dependentVariablesDependency_.at( 1 ) ||
            dependentVariablesDependency_.at( 2 ) ) )
    {
        throw std::runtime_error( "Error, tabulated atmosphere must be initialized with at least "
                                  "density, pressure and temperature." );
    }

    // Assign values to default boundary handling methods
    if ( boundaryHandling_.empty( ) )
    {
        boundaryHandling_ = std::vector< BoundaryInterpolationType >( numberOfIndependentVariables_,
                                                                      use_boundary_value );
    }
    else
    {
        if ( boundaryHandling_.size( ) != numberOfIndependentVariables_ )
        {
            throw std::runtime_error( "Error, in tabulated atmosphere. Number of boundary handling methods provided does "
                                      "not match number of independent variables." );
        }
    }

    // Assign values to default extrapolation values
    if ( defaultExtrapolationValue_.empty( ) )
    {
        defaultExtrapolationValue_ = std::vector< std::vector< std::pair< double, double > > >(
                    numberOfDependentVariables, std::vector< std::pair< double, double > >(
                        numberOfIndependentVariables_, std::make_pair(
                            IdentityElement::getAdditionIdentity< double >( ),
                            IdentityElement::getAdditionIdentity< double >( ) ) ) );
    }
    else
    {
        if ( defaultExtrapolationValue_.size( ) != numberOfDependentVariables )
        {
            throw std::runtime_error( "Error, in tabulated atmosphere. Number of default extrapolation values provided "
                                      "does not match number of dependent variables." );
        }
    }

    // Create interpolators for variables requested by users, depending on the number of variables
    switch ( numberOfIndependentVariables_ )
    {
    case 1:
    {
        // Call approriate file reading function for 1 independent variables
        Eigen::MatrixXd tabulatedAtmosphereData = input_output::readMatrixFromFile( atmosphereTableFile_.at( 0 ), " \t", "%" );

        // Extract information on file size
        unsigned int numberOfColumnsInFile = tabulatedAtmosphereData.cols( );
        unsigned int numberOfRowsInFile = tabulatedAtmosphereData.rows( );


        // Check whether data is present in the file.
        if ( numberOfRowsInFile < 1 || numberOfColumnsInFile < 1 )
        {
            std::string errorMessage = "Error, in tabulated atmosphere. The atmosphere table file " +
                    atmosphereTableFile_.at( 0 ) + " is empty";
            throw std::runtime_error( errorMessage );
        }

        // Check whether number of dependent variables matches number of columns
        if ( numberOfDependentVariables != ( numberOfColumnsInFile - 1 ) )
        {
            throw std::runtime_error( "Error, in tabulated atmosphere. "
                                      "Number of specified dependent variables does not match file." );
        }

        // Assign sizes to vectors
        independentVariablesData_.resize( numberOfIndependentVariables_ );
        std::vector< std::vector< double > > dependentVariablesData;
        dependentVariablesData.resize( numberOfDependentVariables );

        // Extract variables from file
        for ( unsigned int i = 0; i < numberOfRowsInFile; i++ )
        {
            independentVariablesData_.at( 0 ).push_back( tabulatedAtmosphereData( i, 0 ) );
            for ( unsigned int j = 0; j < dependentVariablesDependency_.size( ); j++ )
            {
                if ( dependentVariablesDependency_.at( j ) )
                {
                    dependentVariablesData.at( dependentVariableIndices_.at( j ) ).push_back(
                                tabulatedAtmosphereData( i, dependentVariableIndices_.at( j ) + 1 ) );
                }
            }
        }

        // Create interpolators for density, pressure and temperature
        interpolatorForDensity_ = std::make_shared< CubicSplineInterpolatorDouble >(
                    independentVariablesData_.at( 0 ), dependentVariablesData.at( dependentVariableIndices_.at( 0 ) ),
                    huntingAlgorithm, boundaryHandling_.at( 0 ), defaultExtrapolationValue_.at( dependentVariableIndices_.at( 0 ) ).at( 0 ) );
        interpolatorForPressure_ = std::make_shared< CubicSplineInterpolatorDouble >(
                    independentVariablesData_.at( 0 ), dependentVariablesData.at( dependentVariableIndices_.at( 1 ) ),
                    huntingAlgorithm, boundaryHandling_.at( 0 ), defaultExtrapolationValue_.at( dependentVariableIndices_.at( 1 ) ).at( 0 ) );
        interpolatorForTemperature_ = std::make_shared< CubicSplineInterpolatorDouble >(
                    independentVariablesData_.at( 0 ), dependentVariablesData.at( dependentVariableIndices_.at( 2 ) ),
                    huntingAlgorithm, boundaryHandling_.at( 0 ), defaultExtrapolationValue_.at( dependentVariableIndices_.at( 2 ) ).at( 0 ) );

        // Create remaining interpolators, if requested by user
        if ( dependentVariablesDependency_.at( 3 ) )
        {
            interpolatorForGasConstant_ = std::make_shared< CubicSplineInterpolatorDouble >(
                        independentVariablesData_.at( 0 ), dependentVariablesData.at( dependentVariableIndices_.at( 3 ) ),
                        huntingAlgorithm, boundaryHandling_.at( 0 ), defaultExtrapolationValue_.at( dependentVariableIndices_.at( 3 ) ).at( 0 ) );
        }
        if ( dependentVariablesDependency_.at( 4 ) )
        {
            interpolatorForSpecificHeatRatio_ = std::make_shared< CubicSplineInterpolatorDouble >(
                        independentVariablesData_.at( 0 ), dependentVariablesData.at( dependentVariableIndices_.at( 4 ) ),
                        huntingAlgorithm, boundaryHandling_.at( 0 ), defaultExtrapolationValue_.at( dependentVariableIndices_.at( 4 ) ).at( 0 ) );
        }
        if ( dependentVariablesDependency_.at( 5 ) )
        {
            interpolatorForMolarMass_ = std::make_shared< CubicSplineInterpolatorDouble >(
                        independentVariablesData_.at( 0 ), dependentVariablesData.at( dependentVariableIndices_.at( 5 ) ),
                        huntingAlgorithm, boundaryHandling_.at( 0 ), defaultExtrapolationValue_.at( dependentVariableIndices_.at( 5 ) ).at( 0 ) );
        }
        break;
    }
    case 2:
    {
        createMultiDimensionalAtmosphereInterpolators< 2 >( );
        break;
    }
    case 3:
    {
        createMultiDimensionalAtmosphereInterpolators< 3 >( );
        break;
    }
    case 4:
    {
        createMultiDimensionalAtmosphereInterpolators< 4 >( );
        break;
    }
    default:
        throw std::runtime_error( "Error, multi-dimensional atmosphere with N>4 not yet implemented" );
    }
}

//! Create interpolators for specified dependent variables, taking into consideration the number
//! of independent variables (which is greater than one).
template< unsigned int NumberOfIndependentVariables >
void TabulatedAtmosphere::createMultiDimensionalAtmosphereInterpolators( )
{
    // Call approriate file reading function for N independent variables
    std::pair< std::vector< boost::multi_array< double, static_cast< size_t >( NumberOfIndependentVariables ) > >,
            std::vector< std::vector< double > > > tabulatedAtmosphereData;

    // Extract data
    tabulatedAtmosphereData = input_output::readTabulatedAtmosphere< NumberOfIndependentVariables >( atmosphereTableFile_ );

    // Assign independent variables
    independentVariablesData_ = tabulatedAtmosphereData.second;

    // Create interpolators for density, pressure and temperature
    interpolatorForDensity_ =
            std::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 0 ) ),
                huntingAlgorithm, boundaryHandling_, defaultExtrapolationValue_.at( dependentVariableIndices_.at( 0 ) ) );
    interpolatorForPressure_ =
            std::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 1 ) ),
                huntingAlgorithm, boundaryHandling_, defaultExtrapolationValue_.at( dependentVariableIndices_.at( 1 ) ) );
    interpolatorForTemperature_ =
            std::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 2 ) ),
                huntingAlgorithm, boundaryHandling_, defaultExtrapolationValue_.at( dependentVariableIndices_.at( 2 ) ) );

    // Create remaining interpolators, if requested by user
    if ( dependentVariablesDependency_.at( 3 ) )
    {
        interpolatorForGasConstant_ =
                std::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                    independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 3 ) ),
                    huntingAlgorithm, boundaryHandling_, defaultExtrapolationValue_.at( dependentVariableIndices_.at( 3 ) ) );
    }
    if ( dependentVariablesDependency_.at( 4 ) )
    {
        interpolatorForSpecificHeatRatio_ =
                std::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                    independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 4 ) ),
                    huntingAlgorithm, boundaryHandling_, defaultExtrapolationValue_.at( dependentVariableIndices_.at( 4 ) ) );
    }
    if ( dependentVariablesDependency_.at( 5 ) )
    {
        interpolatorForMolarMass_ =
                std::make_shared< MultiLinearInterpolator< double, double, NumberOfIndependentVariables > >(
                    independentVariablesData_, tabulatedAtmosphereData.first.at( dependentVariableIndices_.at( 5 ) ),
                    huntingAlgorithm, boundaryHandling_, defaultExtrapolationValue_.at( dependentVariableIndices_.at( 5 ) ) );
    }
}

} // namespace aerodynamics

} // namespace tudat
