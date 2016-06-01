/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/SimulationSetup/createFlightConditions.h"

namespace tudat
{

namespace simulation_setup
{

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
    coefficientInterface->updateCurrentCoefficients( std::vector< double >( ) );

    return coefficientInterface;
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
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 1 >(
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
        const boost::function< double( ) > bankAngleFunction )
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

    // Create function to calculate the altitude from current body-fixed state
    boost::function< double( const Eigen::Vector3d ) > altitudeFunction =
            boost::bind( &basic_astrodynamics::BodyShapeModel::getAltitude,
                         centralBody->getShapeModel( ), _1 );

    // Create function to rotate state from intertial to body-fixed frame.
    boost::function< Eigen::Quaterniond( ) > rotationToFrameFunction =
            boost::bind( &Body::getCurrentRotationToLocalFrame, centralBody );
    boost::function< Eigen::Matrix3d( ) > rotationMatrixToFrameDerivativeFunction =
            boost::bind( &Body::getCurrentRotationMatrixDerivativeToLocalFrame, centralBody );
    boost::function< basic_mathematics::Vector6d( const basic_mathematics::Vector6d& ) >
            transformationToCentralBodyFrame =
            boost::bind(
                static_cast< basic_mathematics::Vector6d(&)(
                    const basic_mathematics::Vector6d&,
                    const boost::function< Eigen::Quaterniond( ) >,
                    const boost::function< Eigen::Matrix3d( ) > ) >(
                    &ephemerides::transformStateToFrame ),
                _1, rotationToFrameFunction,
                rotationMatrixToFrameDerivativeFunction );

    // Create flight conditions.
    boost::shared_ptr< aerodynamics::FlightConditions > flightConditions =
            boost::make_shared< aerodynamics::FlightConditions >(
                centralBody->getAtmosphereModel( ), altitudeFunction,
                boost::bind( &Body::getState, bodyWithFlightConditions ),
                boost::bind( &Body::getState, centralBody ),
                transformationToCentralBodyFrame,
                bodyWithFlightConditions->getAerodynamicCoefficientInterface( ) );

    // Create aerodynamic angles calculator and set in flight conditions.
    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
            boost::make_shared< reference_frames::AerodynamicAngleCalculator >(
                boost::bind( &aerodynamics::FlightConditions::getCurrentBodyCenteredBodyFixedState,
                             flightConditions ),
                angleOfAttackFunction, angleOfSideslipFunction, bankAngleFunction, 1 );
    flightConditions->setAerodynamicAngleCalculator( aerodynamicAngleCalculator );

    return flightConditions;


}

} // namespace simulation_setup

} // namespace tudat
