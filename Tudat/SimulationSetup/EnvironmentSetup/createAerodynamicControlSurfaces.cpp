/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EnvironmentSetup/createAerodynamicControlSurfaces.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files
std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings >
readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::map< int, std::string > momentCoefficientFiles,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames )
{
    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( forceCoefficientFiles.begin( )->second );

    // Call approriate file reading function for N independent variables
    std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > coefficientSettings;
    if( numberOfIndependentVariables == 1 )
    {
        coefficientSettings = readGivenSizeTabulatedControlIncrementAerodynamicCoefficientsFromFiles< 1 >(
                    forceCoefficientFiles, momentCoefficientFiles, independentVariableNames );
    }
    else if( numberOfIndependentVariables == 2 )
    {
        coefficientSettings = readGivenSizeTabulatedControlIncrementAerodynamicCoefficientsFromFiles< 2 >(
                    forceCoefficientFiles, momentCoefficientFiles, independentVariableNames );
    }
    else if( numberOfIndependentVariables == 3 )
    {
        coefficientSettings = readGivenSizeTabulatedControlIncrementAerodynamicCoefficientsFromFiles< 3 >(
                    forceCoefficientFiles, momentCoefficientFiles, independentVariableNames );

    }
    else
    {
        throw std::runtime_error( "Error when reading aerodynamic control increment coefficient settings from file, found " +
                                  std::to_string( numberOfIndependentVariables ) +
                                  " independent variables, up to 3 currently supported" );
    }

    return coefficientSettings;
}

//! Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files
std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings >
readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames )
{
    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( forceCoefficientFiles.begin( )->second );

    // Call approriate file reading function for N independent variables
    std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > coefficientSettings;
    if( numberOfIndependentVariables == 1 )
    {
        coefficientSettings = readGivenSizeTabulatedControlIncrementAerodynamicCoefficientsFromFiles< 1 >(
                    forceCoefficientFiles, independentVariableNames );
    }
    else if( numberOfIndependentVariables == 2 )
    {
        coefficientSettings = readGivenSizeTabulatedControlIncrementAerodynamicCoefficientsFromFiles< 2 >(
                    forceCoefficientFiles, independentVariableNames );
    }
    else if( numberOfIndependentVariables == 3 )
    {
        coefficientSettings = readGivenSizeTabulatedControlIncrementAerodynamicCoefficientsFromFiles< 3 >(
                    forceCoefficientFiles, independentVariableNames );
    }
    else
    {
        throw std::runtime_error( "Error when reading aerodynamic coefficient settings from file, found " +
                                  std::to_string( numberOfIndependentVariables ) +
                                  " independent variables, up to 3 currently supported" );
    }
    return coefficientSettings;
}


//! Function to tabulated control surface aerodynamic coefficients from associated settings object
std::shared_ptr< aerodynamics::ControlSurfaceIncrementAerodynamicInterface >
createControlSurfaceIncrementAerodynamicCoefficientInterface(
        const std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    using namespace tudat::aerodynamics;

    std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > coefficientInterface;

    // Check type of interface that is to be created.
    switch( coefficientSettings->getAerodynamicCoefficientType( ) )
    {
    case tabulated_coefficients:
    {
        // Check number of dimensions of tabulated coefficients.
        int numberOfDimensions = coefficientSettings->getIndependentVariableNames( ).size( );
        switch( numberOfDimensions )
        {
        case 1:
        {
            coefficientInterface = createTabulatedControlSurfaceIncrementAerodynamicCoefficientInterface< 1 >(
                        coefficientSettings, body );
            break;
        }
        case 2:
        {
            coefficientInterface = createTabulatedControlSurfaceIncrementAerodynamicCoefficientInterface< 2 >(
                        coefficientSettings, body );
            break;
        }
        case 3:
        {
            coefficientInterface = createTabulatedControlSurfaceIncrementAerodynamicCoefficientInterface< 3 >(
                        coefficientSettings, body );
            break;
        }
        case 4:
        {
            coefficientInterface = createTabulatedControlSurfaceIncrementAerodynamicCoefficientInterface< 4 >(
                        coefficientSettings, body );
            break;
        }
        case 5:
        {
            coefficientInterface = createTabulatedControlSurfaceIncrementAerodynamicCoefficientInterface< 5 >(
                        coefficientSettings, body );
            break;
        }
        case 6:
        {
            coefficientInterface = createTabulatedControlSurfaceIncrementAerodynamicCoefficientInterface< 6 >(
                        coefficientSettings, body );
            break;
        }
        default:
            throw std::runtime_error( "Error when making tabulated control surface aerodynamic coefficient interface, " +
                                      std::to_string( numberOfDimensions ) + " dimensions not yet implemented" );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error, do not recognize control surface aerodynamic coefficient settings for " + body );
    }


    return coefficientInterface;
}

}

}
