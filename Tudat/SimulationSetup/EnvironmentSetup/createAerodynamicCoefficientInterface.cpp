/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/lambda/lambda.hpp>

#include "Tudat/SimulationSetup/EnvironmentSetup/createAerodynamicCoefficientInterface.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create aerodynamic coefficient settings from coefficients stored in data files
std::shared_ptr< AerodynamicCoefficientSettings > readTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::map< int, std::string > momentCoefficientFiles,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame,
        const bool areCoefficientsInNegativeAxisDirection,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( forceCoefficientFiles.begin( )->second );

    // Call approriate file reading function for N independent variables
    std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings;
    if( numberOfIndependentVariables == 1 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 1 >(
                    forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea, lateralReferenceLength,
                    momentReferencePoint, independentVariableNames, areCoefficientsInAerodynamicFrame,
                    areCoefficientsInNegativeAxisDirection, interpolatorSettings );
    }
    else if( numberOfIndependentVariables == 2 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 2 >(
                    forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea, lateralReferenceLength,
                    momentReferencePoint, independentVariableNames, areCoefficientsInAerodynamicFrame,
                    areCoefficientsInNegativeAxisDirection, interpolatorSettings );
    }
    else if( numberOfIndependentVariables == 3 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 3 >(
                    forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea, lateralReferenceLength,
                    momentReferencePoint, independentVariableNames, areCoefficientsInAerodynamicFrame,
                    areCoefficientsInNegativeAxisDirection, interpolatorSettings );
    }
    else
    {
        throw std::runtime_error( "Error when reading aerodynamic coefficient settings from file, found " +
                                  std::to_string( numberOfIndependentVariables ) +
                                  " independent variables, up to 3 currently supported" );
    }
    return coefficientSettings;
}

//! Function to create aerodynamic coefficient settings from coefficients stored in data files
std::shared_ptr< AerodynamicCoefficientSettings >
readTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const double referenceArea,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame,
        const bool areCoefficientsInNegativeAxisDirection,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( forceCoefficientFiles.begin( )->second );

    // Call approriate file reading function for N independent variables
    std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings;
    if( numberOfIndependentVariables == 1 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 1 >(
                    forceCoefficientFiles, referenceArea, independentVariableNames,
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection, interpolatorSettings );
    }
    else if( numberOfIndependentVariables == 2 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 2 >(
                    forceCoefficientFiles, referenceArea, independentVariableNames,
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection, interpolatorSettings );
    }
    else if( numberOfIndependentVariables == 3 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 3 >(
                    forceCoefficientFiles, referenceArea, independentVariableNames,
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection, interpolatorSettings );
    }
    else
    {
        throw std::runtime_error( "Error when reading aerodynamic coefficient settings from file, found " +
                                  std::to_string( numberOfIndependentVariables ) +
                                  " independent variables, up to 3 currently supported" );
    }
    return coefficientSettings;
}

//! Function to create an aerodynamic coefficient interface containing constant coefficients.
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createConstantCoefficientAerodynamicCoefficientInterface(
        const Eigen::Vector3d constantForceCoefficient,
        const Eigen::Vector3d constantMomentCoefficient,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool areCoefficientsInAerodynamicFrame,
        const bool areCoefficientsInNegativeAxisDirection )
{
    // Create coefficient interface
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface =
            std::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                [ = ]( const std::vector< double >& ){ return constantForceCoefficient; },
                [ = ]( const std::vector< double >& ){ return constantMomentCoefficient; },
                referenceLength, referenceArea, lateralReferenceLength, momentReferencePoint,
                std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >( ),
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
    coefficientInterface->updateFullCurrentCoefficients( std::vector< double >( ) );
    return coefficientInterface;
}

//! Factory function for tabulated (1-D independent variables) aerodynamic coefficient interface from coefficient settings.
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createUnivariateTabulatedCoefficientAerodynamicCoefficientInterface(
        const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    using namespace tudat::interpolators;

    // Check consistency of type.
    std::shared_ptr< TabulatedAerodynamicCoefficientSettings< 1 > > tabulatedCoefficientSettings =
            std::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< 1 > >( coefficientSettings );
    if( tabulatedCoefficientSettings == nullptr )
    {
        throw std::runtime_error(
                    "Error, expected tabulated aerodynamic coefficients of size " +
                    std::to_string( 1 ) + "for body " + body );
    }
    else
    {

        // Retrieve or generate interpolation settings
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector3d > > forceInterpolator;
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector3d > > momentInterpolator;
        if ( tabulatedCoefficientSettings->getInterpolatorSettings( ) == nullptr )
        {
            forceInterpolator = createOneDimensionalInterpolator( tabulatedCoefficientSettings->getForceCoefficients( ),
                                                                  std::make_shared< InterpolatorSettings >( linear_interpolator ) );
            momentInterpolator = createOneDimensionalInterpolator( tabulatedCoefficientSettings->getForceCoefficients( ),
                                                                   std::make_shared< InterpolatorSettings >( linear_interpolator ) );
        }
        else
        {
            forceInterpolator = createOneDimensionalInterpolator( tabulatedCoefficientSettings->getForceCoefficients( ),
                                                                  std::dynamic_pointer_cast< InterpolatorSettings >(
                                                                      tabulatedCoefficientSettings->getInterpolatorSettings( ) ) );
            momentInterpolator = createOneDimensionalInterpolator( tabulatedCoefficientSettings->getForceCoefficients( ),
                                                                   std::dynamic_pointer_cast< InterpolatorSettings >(
                                                                       tabulatedCoefficientSettings->getInterpolatorSettings( ) ) );
        }

        // Create aerodynamic coefficient interface.
        return  std::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                    std::bind( &Interpolator< double, Eigen::Vector3d >::interpolate, forceInterpolator, std::placeholders::_1 ),
                    std::bind( &Interpolator< double, Eigen::Vector3d >::interpolate, momentInterpolator, std::placeholders::_1 ),
                    tabulatedCoefficientSettings->getReferenceLength( ),
                    tabulatedCoefficientSettings->getReferenceArea( ),
                    tabulatedCoefficientSettings->getReferenceLength( ),
                    tabulatedCoefficientSettings->getMomentReferencePoint( ),
                    tabulatedCoefficientSettings->getIndependentVariableNames( ),
                    tabulatedCoefficientSettings->getAreCoefficientsInAerodynamicFrame( ),
                    tabulatedCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( ) );
    }
}

//! Function to create and aerodynamic coefficient interface.
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createAerodynamicCoefficientInterface(
        const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    using namespace tudat::aerodynamics;

    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface;

    // Check type of interface that is to be created.
    switch( coefficientSettings->getAerodynamicCoefficientType( ) )
    {
    case constant_aerodynamic_coefficients:
    {
        // Check consistency of type.
        std::shared_ptr< ConstantAerodynamicCoefficientSettings > constantCoefficientSettings =
                std::dynamic_pointer_cast< ConstantAerodynamicCoefficientSettings >(
                    coefficientSettings );
        if( constantCoefficientSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected constant aerodynamic coefficients for body " + body );
        }
        else
        {
            // create constant interface.
            coefficientInterface = createConstantCoefficientAerodynamicCoefficientInterface(
                        constantCoefficientSettings->getConstantForceCoefficient( ),
                        constantCoefficientSettings->getConstantMomentCoefficient( ),
                        constantCoefficientSettings->getReferenceLength( ),
                        constantCoefficientSettings->getReferenceArea( ),
                        constantCoefficientSettings->getReferenceLength( ),
                        constantCoefficientSettings->getMomentReferencePoint( ),
                        constantCoefficientSettings->getAreCoefficientsInAerodynamicFrame( ),
                        constantCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( ) );
        }
        break;
    }
    case tabulated_coefficients:
    {
        // Check number of dimensions of tabulated coefficients.
        int numberOfDimensions = coefficientSettings->getIndependentVariableNames( ).size( );
        switch( numberOfDimensions )
        {
        case 1:
        {
            coefficientInterface = createUnivariateTabulatedCoefficientAerodynamicCoefficientInterface(
                        coefficientSettings, body );
            break;
        }
        case 2:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 2 >(
                        coefficientSettings, body );
            break;
        }
        case 3:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 3 >(
                        coefficientSettings, body );
            break;
        }
        case 4:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 4 >(
                        coefficientSettings, body );
            break;
        }
        case 5:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 5 >(
                        coefficientSettings, body );
            break;
        }
        case 6:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 6 >(
                        coefficientSettings, body );
            break;
        }
        default:
            throw std::runtime_error( "Error when making tabulated aerodynamic coefficient interface, " +
                                      std::to_string( numberOfDimensions ) + " dimensions not yet implemented" );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error, do not recognize aerodynamic coefficient settings for " + body );
    }

    // Create and set control surfaces
    if( coefficientSettings->getControlSurfaceSettings( ).size( ) != 0 )
    {
        std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > >
                controlSurfaceIncrementInterfaces;
        std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > >
                controlSurfaceSettings = coefficientSettings->getControlSurfaceSettings( );
        for( std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > >::iterator
             settingIterator = controlSurfaceSettings.begin( ); settingIterator != controlSurfaceSettings.end( );
             settingIterator++ )
        {
            controlSurfaceIncrementInterfaces[ settingIterator->first ] =
                    createControlSurfaceIncrementAerodynamicCoefficientInterface(
                        settingIterator->second, body );
        }
        coefficientInterface->setControlSurfaceIncrements( controlSurfaceIncrementInterfaces );

    }

    return coefficientInterface;
}

} // simulation_setup

} // tudat
