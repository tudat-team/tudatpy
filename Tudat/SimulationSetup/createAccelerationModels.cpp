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
#include <boost/lexical_cast.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/SimulationSetup/accelerationSettings.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/createFlightConditions.h"

namespace tudat
{

namespace simulation_setup
{

using namespace aerodynamics;
using namespace gravitation;
using namespace basic_astrodynamics;
using namespace electro_magnetism;
using namespace ephemerides;

//! Function to add to double-returning functions.
double evaluateDoubleFunctions(
        const boost::function< double( ) >& function1,
        const boost::function< double( ) >& function2 )
{
    return function1( ) + function2( );
}


//! Function to create central gravity acceleration model.
boost::shared_ptr< CentralGravitationalAccelerationModel3d > createCentralGravityAcceleratioModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const bool useCentralBodyFixedFrame )
{
    // Declare pointer to return object.
    boost::shared_ptr< CentralGravitationalAccelerationModel3d > accelerationModelPointer;

    // Check if body is endowed with a gravity field model (i.e. is capable of exerting
    // gravitation acceleration).
    if( bodyExertingAcceleration->getGravityFieldModel( ) == NULL )
    {
        throw std::runtime_error(
                    std::string( "Error, gravity field model not set when making central ") +
                    " gravitational acceleration of " + nameOfBodyExertingAcceleration + " on " +
                    nameOfBodyUndergoingAcceleration );
    }
    else
    {
        boost::function< double( ) > gravitationalParameterFunction;

        // Set correct value for gravitational parameter.
        if( useCentralBodyFixedFrame == 0  ||
                bodyUndergoingAcceleration->getGravityFieldModel( ) == NULL )
        {
            gravitationalParameterFunction =
                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                 bodyExertingAcceleration->getGravityFieldModel( ) );
        }
        else
        {
            boost::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                 bodyExertingAcceleration->getGravityFieldModel( ) );
            boost::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                 bodyUndergoingAcceleration->getGravityFieldModel( ) );
            gravitationalParameterFunction =
                    boost::bind( &evaluateDoubleFunctions,
                                 gravitationalParameterOfBodyExertingAcceleration,
                                 gravitationalParameterOfBodyUndergoingAcceleration );
        }

        // Create acceleration object.
        accelerationModelPointer =
                boost::make_shared< CentralGravitationalAccelerationModel3d >(
                    boost::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                    gravitationalParameterFunction,
                    boost::bind( &Body::getPosition, bodyExertingAcceleration ),
                    useCentralBodyFixedFrame );
    }


    return accelerationModelPointer;
}

//! Function to create spherical harmonic gravity acceleration model.
boost::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModelXd >
createSphericalHarmonicsGravityAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame )
{
    // Declare pointer to return object
    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModelXd > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    boost::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicsSettings =
            boost::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >(
                accelerationSettings );
    if( sphericalHarmonicsSettings == NULL )
    {
        throw std::runtime_error(
                    std::string( "Error, acceleration settings inconsistent ") +
                    " making sh gravitational acceleration of " + nameOfBodyExertingAcceleration +
                    " on " + nameOfBodyUndergoingAcceleration );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        boost::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );

        boost::shared_ptr< RotationalEphemeris> rotationalEphemeris =
                bodyExertingAcceleration->getRotationalEphemeris( );
        if( sphericalHarmonicsGravityField == NULL )
        {
            throw std::runtime_error(
                        std::string( "Error, spherical harmonic gravity field model not set when ")
                        + " making sh gravitational acceleration of " +
                        nameOfBodyExertingAcceleration +
                        " on " + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            if( rotationalEphemeris == NULL )
            {
                throw std::runtime_error( "Warning when making spherical harmonic acceleration on body " +
                                          nameOfBodyUndergoingAcceleration + ", no rotation model found for " +
                                          nameOfBodyExertingAcceleration );
            }

            if( rotationalEphemeris->getTargetFrameOrientation( ) !=
                    sphericalHarmonicsGravityField->getFixedReferenceFrame( ) )
            {
                throw std::runtime_error( "Warning when making spherical harmonic acceleration on body " +
                                          nameOfBodyUndergoingAcceleration + ", rotation model found for " +
                                          nameOfBodyExertingAcceleration + " is incompatible, frames are: " +
                                          rotationalEphemeris->getTargetFrameOrientation( ) + " and " +
                                          sphericalHarmonicsGravityField->getFixedReferenceFrame( ) );
            }

            boost::function< double( ) > gravitationalParameterFunction;

            // Check if mutual acceleration is to be used.
            if( useCentralBodyFixedFrame == false ||
                    bodyUndergoingAcceleration->getGravityFieldModel( ) == NULL )
            {
                gravitationalParameterFunction =
                        boost::bind( &SphericalHarmonicsGravityField::getGravitationalParameter,
                                     sphericalHarmonicsGravityField );
            }
            else
            {
                // Create function returning summed gravitational parameter of the two bodies.
                boost::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                        boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                     sphericalHarmonicsGravityField );
                boost::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                        boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                     bodyUndergoingAcceleration->getGravityFieldModel( ) );
                gravitationalParameterFunction =
                        boost::bind( &evaluateDoubleFunctions,
                                     gravitationalParameterOfBodyExertingAcceleration,
                                     gravitationalParameterOfBodyUndergoingAcceleration );
            }

            // Create acceleration object.
            accelerationModel =
                    boost::make_shared< SphericalHarmonicsGravitationalAccelerationModelXd >
                    ( boost::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                      gravitationalParameterFunction,
                      sphericalHarmonicsGravityField->getReferenceRadius( ),
                      boost::bind( &SphericalHarmonicsGravityField::getCosineCoefficients,
                                   sphericalHarmonicsGravityField,
                                   sphericalHarmonicsSettings->maximumDegree_,
                                   sphericalHarmonicsSettings->maximumOrder_ ),
                      boost::bind( &SphericalHarmonicsGravityField::getSineCoefficients,
                                   sphericalHarmonicsGravityField,
                                   sphericalHarmonicsSettings->maximumDegree_,
                                   sphericalHarmonicsSettings->maximumOrder_ ),
                      boost::bind( &Body::getPosition, bodyExertingAcceleration ),
                      boost::bind( &Body::getCurrentRotationToGlobalFrame,
                                   bodyExertingAcceleration ), useCentralBodyFixedFrame );
        }
    }
    return accelerationModel;
}


//! Function to create a third body central gravity acceleration model.
boost::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration >
createThirdBodyCentralGravityAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody )
{
    // Declare pointer to return object.
    boost::shared_ptr< ThirdBodyCentralGravityAcceleration > accelerationModelPointer;

    // Create acceleration object.
    accelerationModelPointer =  boost::make_shared< ThirdBodyCentralGravityAcceleration >(
                boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createCentralGravityAcceleratioModel( bodyUndergoingAcceleration,
                                                          bodyExertingAcceleration,
                                                          nameOfBodyUndergoingAcceleration,
                                                          nameOfBodyExertingAcceleration, 0 ) ),
                boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createCentralGravityAcceleratioModel( centralBody, bodyExertingAcceleration,
                                                          nameOfCentralBody,
                                                          nameOfBodyExertingAcceleration, 0 ) ), nameOfCentralBody );

    return accelerationModelPointer;
}


//! Function to create an aerodynamic acceleration model.
boost::shared_ptr< aerodynamics::AerodynamicAcceleration > createAerodynamicAcceleratioModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    // Check existence of required environment models
    if( bodyUndergoingAcceleration->getAerodynamicCoefficientInterface( ) == NULL )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, body " +
                                  nameOfBodyUndergoingAcceleration +
                                  "has no aerodynamic coefficients." );
    }

    if( bodyExertingAcceleration->getAtmosphereModel( ) == NULL )
    {
        throw std::runtime_error(  "Error when making aerodynamic acceleration, central body " +
                                   nameOfBodyExertingAcceleration + " has no atmosphere model.");
    }

    if( bodyExertingAcceleration->getShapeModel( ) == NULL )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, central body " +
                                  nameOfBodyExertingAcceleration + " has no shape model." );
    }

    // Retrieve flight conditions; create object if not yet extant.
    boost::shared_ptr< FlightConditions > bodyFlightConditions =
            bodyUndergoingAcceleration->getFlightConditions( );
    if( bodyFlightConditions == NULL )
    {
        bodyUndergoingAcceleration->setFlightConditions(
                    createFlightConditions( bodyUndergoingAcceleration,
                                            bodyExertingAcceleration,
                                            nameOfBodyUndergoingAcceleration,
                                            nameOfBodyExertingAcceleration ) );
        bodyFlightConditions = bodyUndergoingAcceleration->getFlightConditions( );
    }

    // Retrieve frame in which aerodynamic coefficients are defined.
    boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > aerodynamicCoefficients =
            bodyUndergoingAcceleration->getAerodynamicCoefficientInterface( );
    reference_frames::AerodynamicsReferenceFrames accelerationFrame;
    if( aerodynamicCoefficients->getAreCoefficientsInAerodynamicFrame( ) )
    {
        accelerationFrame = reference_frames::aerodynamic_frame;
    }
    else
    {
        accelerationFrame = reference_frames::body_frame;
    }

    // Create function to transform from frame of aerodynamic coefficienrs to that of propagation.
    boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > toPropagationFrameTransformation;
    toPropagationFrameTransformation =
            reference_frames::getAerodynamicForceTransformationFunction(
                bodyFlightConditions->getAerodynamicAngleCalculator( ),
                accelerationFrame,
                boost::bind( &Body::getCurrentRotationToGlobalFrame, bodyExertingAcceleration ),
                reference_frames::inertial_frame );

    boost::function< Eigen::Vector3d( ) > coefficientFunction =
            boost::bind( &AerodynamicCoefficientInterface::getCurrentForceCoefficients,
                         aerodynamicCoefficients );
    boost::function< Eigen::Vector3d( ) > coefficientInPropagationFrameFunction =
            boost::bind( static_cast< Eigen::Vector3d(&)(
                             const boost::function< Eigen::Vector3d( ) >,
                             const boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > ) >(
                             &reference_frames::transformVector ),
                         coefficientFunction, toPropagationFrameTransformation );

    // Create acceleration model.
    return boost::make_shared< AerodynamicAcceleration >(
                coefficientInPropagationFrameFunction,
                boost::bind( &FlightConditions::getCurrentDensity, bodyFlightConditions ),
                boost::bind( &FlightConditions::getCurrentAirspeed, bodyFlightConditions ),
                boost::bind( &Body::getBodyMass, bodyUndergoingAcceleration ),
                boost::bind( &AerodynamicCoefficientInterface::getReferenceArea,
                             aerodynamicCoefficients ),
                aerodynamicCoefficients->getAreCoefficientsInNegativeAxisDirection( ) );
}

//! Function to create a cannonball radiation pressure acceleration model.
boost::shared_ptr< CannonBallRadiationPressureAcceleration >
createCannonballRadiationPressureAcceleratioModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    // Retrieve radiation pressure interface
    if( bodyUndergoingAcceleration->getRadiationPressureInterfaces( ).count(
                nameOfBodyExertingAcceleration ) == 0 )
    {
        throw std::runtime_error(
                    "Error when making radiation pressure, no radiation pressure interface found  in " +
                    nameOfBodyUndergoingAcceleration +
                    " for body " + nameOfBodyExertingAcceleration );
    }
    boost::shared_ptr< RadiationPressureInterface > radiationPressureInterface =
            bodyUndergoingAcceleration->getRadiationPressureInterfaces( ).at(
                nameOfBodyExertingAcceleration );

    // Create acceleration model.
    return boost::make_shared< CannonBallRadiationPressureAcceleration >(
                boost::bind( &Body::getPosition, bodyExertingAcceleration ),
                boost::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                boost::bind( &RadiationPressureInterface::getCurrentRadiationPressure, radiationPressureInterface ),
                boost::bind( &RadiationPressureInterface::getRadiationPressureCoefficient, radiationPressureInterface ),
                boost::bind( &RadiationPressureInterface::getArea, radiationPressureInterface ),
                boost::bind( &Body::getBodyMass, bodyUndergoingAcceleration ) );

}


//! Function to create acceleration model object.
boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > createAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfCentralBody )
{
    // Declare pointer to return object.
    boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;

    // Switch to call correct acceleration model type factory function.
    switch( accelerationSettings->accelerationType_ )
    {
    case central_gravity:
        // Check if body is a single-body central gravity acceleration (use third-body if not)
        if( nameOfCentralBody == nameOfBodyExertingAcceleration ||
                isFrameInertial( nameOfCentralBody ) )
        {
            // Check if gravitational parameter to use is sum of gravitational paramater of the
            // two bodies.
            bool useCentralBodyFixedFrame = 0;
            if( nameOfCentralBody == nameOfBodyExertingAcceleration )
            {
                useCentralBodyFixedFrame = 1;
            }

            accelerationModelPointer = createCentralGravityAcceleratioModel(
                        bodyUndergoingAcceleration,
                        bodyExertingAcceleration,
                        nameOfBodyUndergoingAcceleration,
                        nameOfBodyExertingAcceleration, useCentralBodyFixedFrame );
        }
        // Create third body central gravity acceleration
        else
        {

            accelerationModelPointer = createThirdBodyCentralGravityAccelerationModel(
                        bodyUndergoingAcceleration,
                        bodyExertingAcceleration,
                        centralBody,
                        nameOfBodyUndergoingAcceleration,
                        nameOfBodyExertingAcceleration,
                        nameOfCentralBody );
        }
        break;
    case spherical_harmonic_gravity:
        if( nameOfCentralBody == nameOfBodyExertingAcceleration ||
                isFrameInertial( nameOfCentralBody ) )
        {
            // Check if gravitational parameter to use is sum of gravitational paramater of the
            // two bodies.
            bool useCentralBodyFixedFrame = 0;
            if( nameOfCentralBody == nameOfBodyExertingAcceleration )
            {
                useCentralBodyFixedFrame = 1;
            }

            accelerationModelPointer = createSphericalHarmonicsGravityAcceleration(
                        bodyUndergoingAcceleration,
                        bodyExertingAcceleration,
                        nameOfBodyUndergoingAcceleration,
                        nameOfBodyExertingAcceleration,
                        accelerationSettings, useCentralBodyFixedFrame );

        }
        else
        {
            throw std::runtime_error(
                        "Error, cannot yet make third body spherical harmonic acceleration." );

        }
        break;
    case aerodynamic:
        accelerationModelPointer = createAerodynamicAcceleratioModel(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration );
        break;
    case cannon_ball_radiation_pressure:
        accelerationModelPointer = createCannonballRadiationPressureAcceleratioModel(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration );
        break;
    default:
        throw std::runtime_error(
                    std::string( "Error, acceleration model ") +
                    boost::lexical_cast< std::string >( accelerationSettings->accelerationType_ ) +
                    " not recognized when making acceleration model of" +
                    nameOfBodyExertingAcceleration + " on " +
                    nameOfBodyUndergoingAcceleration );
        break;
    }
    return accelerationModelPointer;
}

//! Function to create a set of acceleration models from a map of bodies and acceleration model
//! types.
AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::map< std::string, std::string >& centralBodies )
{
    // Declare return map.
    AccelerationMap accelerationModelMap;

    // Iterate over all bodies which are undergoing acceleration
    for( SelectedAccelerationMap::const_iterator bodyIterator =
         selectedAccelerationPerBody.begin( ); bodyIterator != selectedAccelerationPerBody.end( );
         bodyIterator++ )
    {
        boost::shared_ptr< Body > currentCentralBody;

        // Retrieve name of body undergoing acceleration.
        std::string bodyUndergoingAcceleration = bodyIterator->first;

        // Retrieve name of current central body.
        std::string currentCentralBodyName = centralBodies.at( bodyUndergoingAcceleration );

        if( !isFrameInertial( currentCentralBodyName ) )
        {
            if( bodyMap.count( currentCentralBodyName ) == 0 )
            {
                throw std::runtime_error(
                            std::string( "Error, could not find non-inertial central body ") +
                            currentCentralBodyName + " of " + bodyUndergoingAcceleration +
                            " when making acceleration model." );
            }
            else
            {
                currentCentralBody = bodyMap.at( currentCentralBodyName );
            }
        }

        // Check if body undergoing acceleration is included in bodyMap
        if( bodyMap.count( bodyUndergoingAcceleration ) ==  0 )
        {
            throw std::runtime_error(
                        std::string( "Error when making acceleration models, requested forces" ) +
                        "acting on body " + bodyUndergoingAcceleration  +
                        ", but no such body found in map of bodies" );
        }

        // Declare map of acceleration models acting on current body.
        SingleBodyAccelerationMap mapOfAccelerationsForBody;

        // Retrieve list of required acceleration model types and bodies exerting accelerationd on
        // current body.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >
                accelerationsForBody = bodyIterator->second;

        // Iterate over all bodies exerting an acceleration
        for( std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >::
             iterator body2Iterator = accelerationsForBody.begin( );
             body2Iterator != accelerationsForBody.end( ); body2Iterator++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = body2Iterator->first;

            // Check if body exerting acceleration is included in bodyMap
            if( bodyMap.count( bodyExertingAcceleration ) ==  0 )
            {
                throw std::runtime_error(
                            std::string( "Error when making acceleration models, requested forces ")
                            + "acting on body " + bodyUndergoingAcceleration  + " due to body " +
                            bodyExertingAcceleration +
                            ", but no such body found in map of bodies" );
            }

            // Retrieve list of accelerations due to current body.
            std::vector< boost::shared_ptr< AccelerationSettings > > accelerationList =
                    body2Iterator->second;

            for( unsigned int i = 0; i < accelerationList.size( ); i++ )
            {
                // Create acceleration model.
                mapOfAccelerationsForBody[ bodyExertingAcceleration ].push_back(
                            createAccelerationModel( bodyMap.at( bodyUndergoingAcceleration ),
                                                     bodyMap.at( bodyExertingAcceleration ),
                                                     accelerationList.at( i ),
                                                     bodyUndergoingAcceleration,
                                                     bodyExertingAcceleration,
                                                     currentCentralBody,
                                                     currentCentralBodyName ) );
            }
        }

        // Put acceleration models on current body in return map.
        accelerationModelMap[ bodyUndergoingAcceleration ] = mapOfAccelerationsForBody;
    }

    return accelerationModelMap;
}

//! Function to create acceleration models from a map of bodies and acceleration model types.
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::vector< std::string >& propagatedBodies,
        const std::vector< std::string >& centralBodies )
{
    if( centralBodies.size( ) != propagatedBodies.size( ) )
    {
        throw std::runtime_error( "Error, number of propagated bodies must equal number of central bodies" );
    }

    std::map< std::string, std::string > centralBodyMap;
    for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
    {
        centralBodyMap[ propagatedBodies.at( i ) ] = centralBodies.at( i );
    }

    return createAccelerationModelsMap( bodyMap, selectedAccelerationPerBody, centralBodyMap );
}


} // namespace simulation_setup

} // namespace tudat
