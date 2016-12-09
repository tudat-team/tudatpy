/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< AerodynamicCoefficientSettings >
readTabulatedAerodynamicCoefficientsFromFiles(
        const std::vector< std::string > forceCoefficientFile,
        const double referenceArea,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame,
        const bool areCoefficientsInNegativeAxisDirection )
{
    // Open file and create file stream.
    std::fstream stream( forceCoefficientFile.at( 0 ).c_str( ), std::ios::in );

    // Check if file opened correctly.
    if ( stream.fail( ) )
    {
        boost::throw_exception(
                    std::runtime_error( boost::str(
                                            boost::format( "Data file '%s' could not be opened." ) %
                                            forceCoefficientFile.at( 0 ).c_str( ) ) ) );
    }

    std::string line;
    std::vector< std::string > vectorOfIndividualStrings;

    int numberOfIndependentVariables = -1;

    while ( !stream.fail( ) && !stream.eof( ) && ( numberOfIndependentVariables < 0 ) )
    {
        // Get line from stream
        std::getline( stream, line );
        boost::algorithm::trim( line );

        if( line.size( ) > 0 && !( line.at( 0 ) == '#' ) )
        {
            boost::algorithm::split( vectorOfIndividualStrings,
                                     line,
                                     boost::algorithm::is_any_of( "\t ;, " ),
                                     boost::algorithm::token_compress_on );
            try
            {
                numberOfIndependentVariables = boost::lexical_cast< double >( vectorOfIndividualStrings.at( 0 ).at( 0 ) );
            }
            catch( std::runtime_error )
            {
                throw std::runtime_error( "Error when reading aerodynamic coefficient settings from file, input is inconsistent" );
            }
        }
    }

    boost::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings;
    if( numberOfIndependentVariables == 1 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 1 >(
                    forceCoefficientFile, referenceArea, independentVariableNames,
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
    }
    else if( numberOfIndependentVariables == 2 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 2 >(
                    forceCoefficientFile, referenceArea, independentVariableNames,
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
    }
    else if( numberOfIndependentVariables == 3 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 3 >(
                    forceCoefficientFile, referenceArea, independentVariableNames,
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
    }
    else
    {
        throw std::runtime_error( "Error when reading aerodynamic coefficient settings from file, found " +
                                  boost::lexical_cast< std::string >( numberOfIndependentVariables ) +
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
                    boost::lexical_cast<  std::string >( 1 ) + "for body " + body );
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
                                      boost::lexical_cast< std::string >( numberOfDimensions ) + " dimensions not yet implemented" );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error, do not recognize aerodynamic coefficient settings for " + body );
    }
    return coefficientInterface;
}

//! Function to create a flight conditions object
boost::shared_ptr< aerodynamics::FlightConditions > createFlightConditions(
        const boost::shared_ptr< Body > bodyWithFlightConditions,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::function< double( ) > angleOfAttackFunction,
        const boost::function< double( ) > angleOfSideslipFunction,
        const boost::function< double( ) > bankAngleFunction,
        const boost::function< void( const double ) > angleUpdateFunction )
{
    // Check whether all required environment models are set.
    if( centralBody->getAtmosphereModel( ) == NULL )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, body " + nameOfBodyExertingAcceleration +
                    " has no atmosphere model." );
    }

    if( centralBody->getShapeModel( ) == NULL )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, body " + nameOfBodyExertingAcceleration +
                    " has no shape model." );
    }

    if( centralBody->getRotationalEphemeris( ) == NULL )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, body " + nameOfBodyExertingAcceleration +
                    " has no rotation model." );
    }

    if( bodyWithFlightConditions->getAerodynamicCoefficientInterface( ) == NULL )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, body "+ nameOfBodyUndergoingAcceleration +
                    " has no aerodynamic coefficients." );
    }


    // Create function to rotate state from intertial to body-fixed frame.
    boost::function< Eigen::Quaterniond( ) > rotationToFrameFunction =
            boost::bind( &Body::getCurrentRotationToLocalFrame, centralBody );
    boost::function< Eigen::Matrix3d( ) > rotationMatrixToFrameDerivativeFunction =
            boost::bind( &Body::getCurrentRotationMatrixDerivativeToLocalFrame, centralBody );

    boost::function< Eigen::Matrix< double, 6, 1 >( ) > bodyStateFunction = boost::bind( &Body::getState, bodyWithFlightConditions );
    boost::function< Eigen::Matrix< double, 6, 1 >( ) > centralBodyStateFunction = boost::bind( &Body::getState, centralBody );

    boost::function< Eigen::Matrix< double, 6, 1 >( ) > relativeBodyFixedStateFunction =
            boost::bind( &ephemerides::transformRelativeStateToFrame< double >,
                         bodyStateFunction, centralBodyStateFunction,
                         rotationToFrameFunction,
                         rotationMatrixToFrameDerivativeFunction );

    // Create aerodynamic angles calculator and set in flight conditions.
    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
            boost::make_shared< reference_frames::AerodynamicAngleCalculator >(
                relativeBodyFixedStateFunction,
                boost::bind( &simulation_setup::Body::getCurrentRotationToGlobalFrame, centralBody ),
                nameOfBodyExertingAcceleration, 1,
                angleOfAttackFunction, angleOfSideslipFunction, bankAngleFunction, angleUpdateFunction );


    // Create flight conditions.
    boost::function< double( const std::string& )> controlSurfaceDeflectionFunction;
    if( bodyWithFlightConditions->getVehicleSystems( ) != NULL )
    {
        controlSurfaceDeflectionFunction = boost::bind(
                    &system_models::VehicleSystems::getCurrentControlSurfaceDeflection,
                    bodyWithFlightConditions->getVehicleSystems( ), _1 );
    }
    boost::shared_ptr< aerodynamics::FlightConditions > flightConditions =
            boost::make_shared< aerodynamics::FlightConditions >(
                centralBody->getAtmosphereModel( ), centralBody->getShapeModel( ),
                bodyWithFlightConditions->getAerodynamicCoefficientInterface( ), aerodynamicAngleCalculator,
                controlSurfaceDeflectionFunction );

    return flightConditions;


}


//! Function to set the angle of attack to trimmed conditions.
boost::shared_ptr< aerodynamics::TrimOrientationCalculator > setTrimmedConditions(
        const boost::shared_ptr< aerodynamics::FlightConditions > flightConditions )
{
    // Create trim object.
    boost::shared_ptr< aerodynamics::TrimOrientationCalculator > trimOrientation =
            boost::make_shared< aerodynamics::TrimOrientationCalculator >(
                flightConditions->getAerodynamicCoefficientInterface( ) );

    // Create angle-of-attack function from trim object.
    boost::function< std::vector< double >( ) > untrimmedIndependentVariablesFunction =
            boost::bind( &aerodynamics::FlightConditions::getAerodynamicCoefficientIndependentVariables,
                         flightConditions );
    boost::function< std::map< std::string, std::vector< double > >( ) > untrimmedControlSurfaceIndependentVariableFunction =
            boost::bind( &aerodynamics::FlightConditions::getControlSurfaceAerodynamicCoefficientIndependentVariables,
                         flightConditions );

    flightConditions->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                boost::bind( &aerodynamics::TrimOrientationCalculator::findTrimAngleOfAttackFromFunction, trimOrientation,
                             untrimmedIndependentVariablesFunction, untrimmedControlSurfaceIndependentVariableFunction ) );

    return trimOrientation;
}

//! Function to set the angle of attack to trimmed conditions.
boost::shared_ptr< aerodynamics::TrimOrientationCalculator > setTrimmedConditions(
        const boost::shared_ptr< Body > bodyWithFlightConditions )
{
    if( bodyWithFlightConditions->getFlightConditions( ) == NULL )
    {
        throw std::runtime_error( "Error, body does not have FlightConditions when setting trim conditions." );
    }

    return setTrimmedConditions( bodyWithFlightConditions->getFlightConditions( ) );
}

//! Function that must be called to link the AerodynamicGuidance object to the simulation
void setGuidanceAnglesFunctions(
        const boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator )
{
    angleCalculator->setOrientationAngleFunctions(
                boost::bind( &aerodynamics::AerodynamicGuidance::getCurrentAngleOfAttack, aerodynamicGuidance ),
                boost::bind( &aerodynamics::AerodynamicGuidance::getCurrentAngleOfSideslip, aerodynamicGuidance ),
                boost::bind( &aerodynamics::AerodynamicGuidance::getCurrentBankAngle, aerodynamicGuidance ),
                boost::bind( &aerodynamics::AerodynamicGuidance::updateGuidance, aerodynamicGuidance,_1 ) );
}

//! Function that must be called to link the AerodynamicGuidance object to the simulation
void setGuidanceAnglesFunctions(
        const boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const boost::shared_ptr< simulation_setup::Body > bodyWithAngles )
{
    boost::shared_ptr< reference_frames::DependentOrientationCalculator >  orientationCalculator =
            bodyWithAngles->getDependentOrientationCalculator( );
    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator =
            boost::dynamic_pointer_cast< reference_frames::AerodynamicAngleCalculator >( orientationCalculator );

    if( angleCalculator == NULL )
    {
        throw std::runtime_error( "Error, body does not have AerodynamicAngleCalculator when setting aerodynamic guidance" );
    }
    else
    {
        setGuidanceAnglesFunctions( aerodynamicGuidance, angleCalculator );
    }
}

} // namespace simulation_setup

} // namespace tudat
