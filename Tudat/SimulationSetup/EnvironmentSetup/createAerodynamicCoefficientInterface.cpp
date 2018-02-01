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


//! Function to create aerodynamic coefficient settings fom coefficients stored in data files
boost::shared_ptr< AerodynamicCoefficientSettings > readTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::map< int, std::string > momentCoefficientFiles,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame ,
        const bool areCoefficientsInNegativeAxisDirection )
{
    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( forceCoefficientFiles.begin( )->second );

    // Call approriate file reading function for N independent variables
    boost::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings;
    if( numberOfIndependentVariables == 1 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 1 >(
                    forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea, lateralReferenceLength,
                    momentReferencePoint, independentVariableNames, areCoefficientsInAerodynamicFrame,
                    areCoefficientsInNegativeAxisDirection );
    }
    else if( numberOfIndependentVariables == 2 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 2 >(
                    forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea, lateralReferenceLength,
                    momentReferencePoint, independentVariableNames, areCoefficientsInAerodynamicFrame,
                    areCoefficientsInNegativeAxisDirection );
    }
    else if( numberOfIndependentVariables == 3 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 3 >(
                    forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea, lateralReferenceLength,
                    momentReferencePoint, independentVariableNames, areCoefficientsInAerodynamicFrame,
                    areCoefficientsInNegativeAxisDirection );
    }
    else
    {
        throw std::runtime_error( "Error when reading aerodynamic coefficient settings from file, found " +
                                  std::to_string( numberOfIndependentVariables ) +
                                  " independent variables, up to 3 currently supported" );
    }
    return coefficientSettings;
}

//! Function to create aerodynamic coefficient settings fom coefficients stored in data files
boost::shared_ptr< AerodynamicCoefficientSettings >
readTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const double referenceArea,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame,
        const bool areCoefficientsInNegativeAxisDirection )
{
    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( forceCoefficientFiles.begin( )->second );

    // Call approriate file reading function for N independent variables
    boost::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings;
    if( numberOfIndependentVariables == 1 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 1 >(
                    forceCoefficientFiles, referenceArea, independentVariableNames,
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
    }
    else if( numberOfIndependentVariables == 2 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 2 >(
                    forceCoefficientFiles, referenceArea, independentVariableNames,
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
    }
    else if( numberOfIndependentVariables == 3 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 3 >(
                    forceCoefficientFiles, referenceArea, independentVariableNames,
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
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
boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createConstantCoefficientAerodynamicCoefficientInterface(
        const Eigen::Vector3d constantForceCoefficient,
        const Eigen::Vector3d constantMomentCoefficient,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool areCoefficientsInAerodynamicFrame,
        const bool areCoefficientsInNegativeAxisDirection  )
{
    // Create coefficient interface
    boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface =
            boost::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                boost::lambda::constant( constantForceCoefficient ),
                boost::lambda::constant( constantMomentCoefficient ),
                referenceLength, referenceArea, lateralReferenceLength, momentReferencePoint,
                std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >( ),
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
    coefficientInterface->updateFullCurrentCoefficients( std::vector< double >( ) );

    return coefficientInterface;
}

//! Factory function for tabulated (1-D independent variables) aerodynamic coefficient interface from coefficient settings.
boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createUnivariateTabulatedCoefficientAerodynamicCoefficientInterface(
        const boost::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    // Check consistency of type.
    boost::shared_ptr< TabulatedAerodynamicCoefficientSettings< 1 > > tabulatedCoefficientSettings =
            boost::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< 1 > >(
                coefficientSettings );
    if( tabulatedCoefficientSettings == NULL )
    {
        throw std::runtime_error(
                    "Error, expected tabulated aerodynamic coefficients of size " +
                    std::to_string( 1 ) + "for body " + body );
    }
    else
    {
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > forceInterpolator =
                interpolators::createOneDimensionalInterpolator(
                    tabulatedCoefficientSettings->getForceCoefficients( ),
                    tabulatedCoefficientSettings->getInterpolationSettings( ) );
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > momentInterpolator =
                interpolators::createOneDimensionalInterpolator(
                    tabulatedCoefficientSettings->getForceCoefficients( ),
                    tabulatedCoefficientSettings->getInterpolationSettings( ) );

        // Create aerodynamic coefficient interface.
        return  boost::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                    boost::bind( &interpolators::Interpolator
                                 < double, Eigen::Vector3d >::interpolate, forceInterpolator, _1 ),
                    boost::bind( &interpolators::Interpolator
                                 < double, Eigen::Vector3d >::interpolate, momentInterpolator, _1 ),
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
boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createAerodynamicCoefficientInterface(
        const boost::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    using namespace tudat::aerodynamics;

    boost::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface;

    // Check type of interface that is to be created.
    switch( coefficientSettings->getAerodynamicCoefficientType( ) )
    {
    case constant_aerodynamic_coefficients:
    {
        // Check consistency of type.
        boost::shared_ptr< ConstantAerodynamicCoefficientSettings > constantCoefficientSettings =
                boost::dynamic_pointer_cast< ConstantAerodynamicCoefficientSettings >(
                    coefficientSettings );
        if( constantCoefficientSettings == NULL )
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
        std::map< std::string, boost::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > >
                controlSurfaceIncrementInterfaces;
        std::map< std::string, boost::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > >
                controlSurfaceSettings = coefficientSettings->getControlSurfaceSettings( );
        for( std::map< std::string, boost::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > >::iterator
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

}

}
