/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>

#include <boost/make_shared.hpp>
#include <boost/bind.hpp>
#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Propulsion/thrustMagnitudeWrapper.h"
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Astrodynamics/Relativity/relativisticAccelerationCorrection.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/Basics/utilities.h"
#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"

namespace tudat
{

namespace simulation_setup
{

using namespace aerodynamics;
using namespace gravitation;
using namespace basic_astrodynamics;
using namespace electro_magnetism;
using namespace ephemerides;


//! Function to create a direct (i.e. not third-body) gravitational acceleration (of any type)
boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > createDirectGravitationalAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfCentralBody,
        const bool isCentralBody )
{
    // Check if sum of gravitational parameters (i.e. inertial force w.r.t. central body) should be used.
    bool sumGravitationalParameters = 0;
    if( ( nameOfCentralBody == nameOfBodyExertingAcceleration ) && bodyUndergoingAcceleration != NULL )
    {
        sumGravitationalParameters = 1;
    }


    // Check type of acceleration model and create.
    boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel;
    switch( accelerationSettings->accelerationType_ )
    {
    case central_gravity:
        accelerationModel = createCentralGravityAcceleratioModel(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    sumGravitationalParameters );
        break;
    case spherical_harmonic_gravity:
        accelerationModel = createSphericalHarmonicsGravityAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings,
                    sumGravitationalParameters );
        break;
    case mutual_spherical_harmonic_gravity:
        accelerationModel = createMutualSphericalHarmonicsGravityAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings,
                    sumGravitationalParameters,
                    isCentralBody );
        break;
    default:

        std::string errorMessage = "Error when making gravitional acceleration model, cannot parse type " +
                std::to_string( accelerationSettings->accelerationType_ );
        throw std::runtime_error( errorMessage );
    }
    return accelerationModel;
}

//! Function to create a third-body gravitational acceleration (of any type)
boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > createThirdBodyGravitationalAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Check type of acceleration model and create.
    boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel;
    switch( accelerationSettings->accelerationType_ )
    {
    case central_gravity:
        accelerationModel = boost::make_shared< ThirdBodyCentralGravityAcceleration >(
                    boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                        createDirectGravitationalAcceleration(
                            bodyUndergoingAcceleration, bodyExertingAcceleration,
                            nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 0 ) ),
                    boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                        createDirectGravitationalAcceleration(
                            centralBody, bodyExertingAcceleration,
                            nameOfCentralBody, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 1 ) ), nameOfCentralBody );
        break;
    case spherical_harmonic_gravity:
        accelerationModel = boost::make_shared< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                    boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            bodyUndergoingAcceleration, bodyExertingAcceleration,
                            nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 0 ) ),
                    boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            centralBody, bodyExertingAcceleration, nameOfCentralBody, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 1 ) ), nameOfCentralBody );
        break;
    case mutual_spherical_harmonic_gravity:
        accelerationModel = boost::make_shared< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                    boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            bodyUndergoingAcceleration, bodyExertingAcceleration,
                            nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 0 ) ),
                    boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            centralBody, bodyExertingAcceleration, nameOfCentralBody, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 1 ) ), nameOfCentralBody );
        break;
    default:

        std::string errorMessage = "Error when making third-body gravitional acceleration model, cannot parse type " +
                std::to_string( accelerationSettings->accelerationType_ );
        throw std::runtime_error( errorMessage );
    }
    return accelerationModel;
}

//! Function to create gravitational acceleration (of any type)
boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > createGravitationalAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfCentralBody )
{

    boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;
    if( accelerationSettings->accelerationType_ != central_gravity &&
            accelerationSettings->accelerationType_ != spherical_harmonic_gravity &&
            accelerationSettings->accelerationType_ != mutual_spherical_harmonic_gravity )
    {
        throw std::runtime_error( "Error when making gravitational acceleration, type is inconsistent" );
    }

    if( nameOfCentralBody == nameOfBodyExertingAcceleration || ephemerides::isFrameInertial( nameOfCentralBody ) )
    {
        accelerationModelPointer = createDirectGravitationalAcceleration( bodyUndergoingAcceleration,
                                                                          bodyExertingAcceleration,
                                                                          nameOfBodyUndergoingAcceleration,
                                                                          nameOfBodyExertingAcceleration,
                                                                          accelerationSettings,
                                                                          nameOfCentralBody, false );
    }
    else
    {
        accelerationModelPointer = createThirdBodyGravitationalAcceleration( bodyUndergoingAcceleration,
                                                                             bodyExertingAcceleration,
                                                                             centralBody,
                                                                             nameOfBodyUndergoingAcceleration,
                                                                             nameOfBodyExertingAcceleration,
                                                                             nameOfCentralBody, accelerationSettings );
    }

    return accelerationModelPointer;
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
                    boost::bind( &utilities::sumFunctionReturn< double >,
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
boost::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel >
createSphericalHarmonicsGravityAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame )
{
    // Declare pointer to return object
    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

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
                        boost::bind( &utilities::sumFunctionReturn< double >,
                                     gravitationalParameterOfBodyExertingAcceleration,
                                     gravitationalParameterOfBodyUndergoingAcceleration );
            }

            // Create acceleration object.
            accelerationModel =
                    boost::make_shared< SphericalHarmonicsGravitationalAccelerationModel >
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

//! Function to create mutual spherical harmonic gravity acceleration model.
boost::shared_ptr< gravitation::MutualSphericalHarmonicsGravitationalAccelerationModel >
createMutualSphericalHarmonicsGravityAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame,
        const bool acceleratedBodyIsCentralBody )
{
    using namespace basic_astrodynamics;

    // Declare pointer to return object
    boost::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    boost::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicsSettings =
            boost::dynamic_pointer_cast< MutualSphericalHarmonicAccelerationSettings >( accelerationSettings );
    if( mutualSphericalHarmonicsSettings == NULL )
    {
        std::string errorMessage = "Error, expected mutual spherical harmonics acceleration settings when making acceleration model on " +
                nameOfBodyUndergoingAcceleration + "due to " + nameOfBodyExertingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        boost::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyExertingAcceleration =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        boost::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyUndergoingAcceleration->getGravityFieldModel( ) );

        if( sphericalHarmonicsGravityFieldOfBodyExertingAcceleration == NULL )
        {

            std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );

        }
        else if( sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration == NULL )
        {

            std::string errorMessage = "Error " + nameOfBodyUndergoingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else
        {
            boost::function< double( ) > gravitationalParameterFunction;

            // Create function returning summed gravitational parameter of the two bodies.
            if( useCentralBodyFixedFrame == false )
            {
                gravitationalParameterFunction =
                        boost::bind( &SphericalHarmonicsGravityField::getGravitationalParameter,
                                     sphericalHarmonicsGravityFieldOfBodyExertingAcceleration );
            }
            else
            {
                // Create function returning summed gravitational parameter of the two bodies.
                boost::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                        boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                     sphericalHarmonicsGravityFieldOfBodyExertingAcceleration );
                boost::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                        boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                     sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration );
                gravitationalParameterFunction =
                        boost::bind( &utilities::sumFunctionReturn< double >,
                                     gravitationalParameterOfBodyExertingAcceleration,
                                     gravitationalParameterOfBodyUndergoingAcceleration );
            }

            // Create acceleration object.

            int maximumDegreeOfUndergoingBody, maximumOrderOfUndergoingBody;
            if( !acceleratedBodyIsCentralBody )
            {
                maximumDegreeOfUndergoingBody = mutualSphericalHarmonicsSettings->maximumDegreeOfBodyUndergoingAcceleration_;
                maximumOrderOfUndergoingBody = mutualSphericalHarmonicsSettings->maximumOrderOfBodyUndergoingAcceleration_;
            }
            else
            {
                maximumDegreeOfUndergoingBody = mutualSphericalHarmonicsSettings->maximumDegreeOfCentralBody_;
                maximumOrderOfUndergoingBody = mutualSphericalHarmonicsSettings->maximumOrderOfCentralBody_;
            }

            accelerationModel = boost::make_shared< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        boost::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                        boost::bind( &Body::getPosition, bodyExertingAcceleration ),
                        gravitationalParameterFunction,
                        sphericalHarmonicsGravityFieldOfBodyExertingAcceleration->getReferenceRadius( ),
                        sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration->getReferenceRadius( ),
                        boost::bind( &SphericalHarmonicsGravityField::getCosineCoefficients,
                                     sphericalHarmonicsGravityFieldOfBodyExertingAcceleration,
                                     mutualSphericalHarmonicsSettings->maximumDegreeOfBodyExertingAcceleration_,
                                     mutualSphericalHarmonicsSettings->maximumOrderOfBodyExertingAcceleration_ ),
                        boost::bind( &SphericalHarmonicsGravityField::getSineCoefficients,
                                     sphericalHarmonicsGravityFieldOfBodyExertingAcceleration,
                                     mutualSphericalHarmonicsSettings->maximumDegreeOfBodyExertingAcceleration_,
                                     mutualSphericalHarmonicsSettings->maximumOrderOfBodyExertingAcceleration_ ),
                        boost::bind( &SphericalHarmonicsGravityField::getCosineCoefficients,
                                     sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration,
                                     maximumDegreeOfUndergoingBody,
                                     maximumOrderOfUndergoingBody ),
                        boost::bind( &SphericalHarmonicsGravityField::getSineCoefficients,
                                     sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration,
                                     maximumDegreeOfUndergoingBody,
                                     maximumOrderOfUndergoingBody ),
                        boost::bind( &Body::getCurrentRotationToGlobalFrame,
                                     bodyExertingAcceleration ),
                        boost::bind( &Body::getCurrentRotationToGlobalFrame,
                                     bodyUndergoingAcceleration ),
                        useCentralBodyFixedFrame );
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

//! Function to create a third body spheric harmonic gravity acceleration model.
boost::shared_ptr< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >
createThirdBodySphericalHarmonicGravityAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings )
{
    using namespace basic_astrodynamics;

    // Declare pointer to return object
    boost::shared_ptr< ThirdBodySphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    boost::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicsSettings =
            boost::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >( accelerationSettings );
    if( sphericalHarmonicsSettings == NULL )
    {
        std::string errorMessage = "Error, expected spherical harmonics acceleration settings when making acceleration model on " +
                nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        boost::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        if( sphericalHarmonicsGravityField == NULL )
        {
            std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making third body spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else
        {

            accelerationModel =  boost::make_shared< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                        boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                            createSphericalHarmonicsGravityAcceleration(
                                bodyUndergoingAcceleration, bodyExertingAcceleration, nameOfBodyUndergoingAcceleration,
                                nameOfBodyExertingAcceleration, sphericalHarmonicsSettings, 0 ) ),
                        boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                            createSphericalHarmonicsGravityAcceleration(
                                centralBody, bodyExertingAcceleration, nameOfCentralBody,
                                nameOfBodyExertingAcceleration, sphericalHarmonicsSettings, 0 ) ), nameOfCentralBody );
        }
    }
    return accelerationModel;
}

//! Function to create a third body mutual spheric harmonic gravity acceleration model.
boost::shared_ptr< gravitation::ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >
createThirdBodyMutualSphericalHarmonicGravityAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Declare pointer to return object
    boost::shared_ptr< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    boost::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicsSettings =
            boost::dynamic_pointer_cast< MutualSphericalHarmonicAccelerationSettings >( accelerationSettings );
    if( mutualSphericalHarmonicsSettings == NULL )
    {

        std::string errorMessage = "Error, expected mutual spherical harmonics acceleration settings when making acceleration model on " +
                nameOfBodyUndergoingAcceleration +
                " due to " + nameOfBodyExertingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        boost::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyExertingAcceleration =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        boost::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyUndergoingAcceleration->getGravityFieldModel( ) );
        boost::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfCentralBody =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    centralBody->getGravityFieldModel( ) );

        if( sphericalHarmonicsGravityFieldOfBodyExertingAcceleration == NULL )
        {
            std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else if( sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration == NULL )
        {
            std::string errorMessage = "Error " + nameOfBodyUndergoingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else if( sphericalHarmonicsGravityFieldOfCentralBody == NULL )
        {
            std::string errorMessage = "Error " + nameOfCentralBody + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else
        {
            boost::shared_ptr< MutualSphericalHarmonicAccelerationSettings > accelerationSettingsForCentralBodyAcceleration =
                    boost::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                        mutualSphericalHarmonicsSettings->maximumDegreeOfBodyExertingAcceleration_,
                        mutualSphericalHarmonicsSettings->maximumOrderOfBodyExertingAcceleration_,
                        mutualSphericalHarmonicsSettings->maximumDegreeOfCentralBody_,
                        mutualSphericalHarmonicsSettings->maximumOrderOfCentralBody_ );
            accelerationModel =  boost::make_shared< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                        boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                            createMutualSphericalHarmonicsGravityAcceleration(
                                bodyUndergoingAcceleration, bodyExertingAcceleration, nameOfBodyUndergoingAcceleration,
                                nameOfBodyExertingAcceleration, mutualSphericalHarmonicsSettings, 0, 0 ) ),
                        boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                            createMutualSphericalHarmonicsGravityAcceleration(
                                centralBody, bodyExertingAcceleration, nameOfCentralBody,
                                nameOfBodyExertingAcceleration, accelerationSettingsForCentralBodyAcceleration, 0, 1 ) ),
                        nameOfCentralBody );
        }
    }
    return accelerationModel;
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
            boost::bind( &reference_frames::transformVectorFunctionFromVectorFunctions,
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

//! Function to create an orbiter relativistic correction acceleration model
boost::shared_ptr< relativity::RelativisticAccelerationCorrection > createRelativisticCorrectionAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const NamedBodyMap& bodyMap )
{
    using namespace relativity;

    // Declare pointer to return object
    boost::shared_ptr< RelativisticAccelerationCorrection > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    boost::shared_ptr< RelativisticAccelerationCorrectionSettings > relativisticAccelerationSettings =
            boost::dynamic_pointer_cast< RelativisticAccelerationCorrectionSettings >(
                accelerationSettings );
    if( relativisticAccelerationSettings == NULL )
    {
        throw std::runtime_error( "Error, expected relativistic acceleration settings when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration );
    }
    else
    {

        // Retrieve function pointers for properties of bodies exerting/undergoing acceleration.
        boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingAcceleration =
                boost::bind( &Body::getState, bodyExertingAcceleration );
        boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingAcceleration =
                boost::bind( &Body::getState, bodyUndergoingAcceleration );

        boost::function< double( ) > centralBodyGravitationalParameterFunction;
        boost::shared_ptr< GravityFieldModel > gravityField = bodyExertingAcceleration->getGravityFieldModel( );
        if( gravityField == NULL )
        {
            throw std::runtime_error( "Error " + nameOfBodyExertingAcceleration + " does not have a gravity field " +
                                      "when making relativistic acceleration on" + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            centralBodyGravitationalParameterFunction =
                    boost::bind( &GravityFieldModel::getGravitationalParameter, bodyExertingAcceleration->getGravityFieldModel( ) );
        }

        // Create acceleration model if only schwarzschild term is to be used.
        if( relativisticAccelerationSettings->calculateLenseThirringCorrection_ == false &&
                relativisticAccelerationSettings->calculateDeSitterCorrection_ == false )
        {
            boost::function< double( ) > ppnGammaFunction = boost::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet );
            boost::function< double( ) > ppnBetaFunction = boost::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet );

            // Create acceleration model.
            accelerationModel = boost::make_shared< RelativisticAccelerationCorrection >
                    ( stateFunctionOfBodyUndergoingAcceleration,
                      stateFunctionOfBodyExertingAcceleration,
                      centralBodyGravitationalParameterFunction,
                      ppnGammaFunction, ppnBetaFunction );

        }
        else
        {

            // Retrieve parameters of primary body if de Sitter term is to be used.
            boost::function< Eigen::Vector6d( ) > stateFunctionOfPrimaryBody;
            boost::function< double( ) > primaryBodyGravitationalParameterFunction;
            if( relativisticAccelerationSettings->calculateDeSitterCorrection_ == true )
            {
                if(  bodyMap.count( relativisticAccelerationSettings->primaryBody_ ) == 0 )
                {
                    throw std::runtime_error( "Error, no primary body " + relativisticAccelerationSettings->primaryBody_ +
                                              " found when making de Sitter acceleration correction" );
                }
                stateFunctionOfPrimaryBody =
                        boost::bind( &Body::getState, bodyMap.at( relativisticAccelerationSettings->primaryBody_ ) );

                if(  bodyMap.at( relativisticAccelerationSettings->primaryBody_ )->getGravityFieldModel( ) == NULL )
                {
                    throw std::runtime_error( "Error, primary body " + relativisticAccelerationSettings->primaryBody_ +
                                              " has no gravity field when making de Sitter acceleration correction" );
                }

                primaryBodyGravitationalParameterFunction =
                        boost::bind( &GravityFieldModel::getGravitationalParameter,
                                     bodyMap.at( relativisticAccelerationSettings->primaryBody_ )->getGravityFieldModel( ) );


            }

            // Retrieve angular momentum vector if Lense-Thirring
            boost::function< Eigen::Vector3d( ) > angularMomentumFunction;
            if( relativisticAccelerationSettings->calculateLenseThirringCorrection_ == true  )
            {
                angularMomentumFunction = boost::lambda::constant(
                            relativisticAccelerationSettings->centralBodyAngularMomentum_ );
            }

            if( relativisticAccelerationSettings->calculateDeSitterCorrection_ == true )
            {
                // Create acceleration model with Lense-Thirring and de Sitter terms.
                accelerationModel = boost::make_shared< RelativisticAccelerationCorrection >
                        ( stateFunctionOfBodyUndergoingAcceleration,
                          stateFunctionOfBodyExertingAcceleration,
                          stateFunctionOfPrimaryBody,
                          centralBodyGravitationalParameterFunction,
                          primaryBodyGravitationalParameterFunction,
                          relativisticAccelerationSettings->primaryBody_,
                          angularMomentumFunction,
                          boost::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet ),
                          boost::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet ),
                          relativisticAccelerationSettings->calculateSchwarzschildCorrection_ );
            }
            else
            {
                // Create acceleration model with Lense-Thirring and term.
                accelerationModel = boost::make_shared< RelativisticAccelerationCorrection >
                        ( stateFunctionOfBodyUndergoingAcceleration,
                          stateFunctionOfBodyExertingAcceleration,
                          centralBodyGravitationalParameterFunction,
                          angularMomentumFunction,
                          boost::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet ),
                          boost::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet ),
                          relativisticAccelerationSettings->calculateSchwarzschildCorrection_ );
            }
        }
    }
    return accelerationModel;
}


//! Function to create empirical acceleration model.
boost::shared_ptr< EmpiricalAcceleration > createEmpiricalAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  boost::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Declare pointer to return object
    boost::shared_ptr< EmpiricalAcceleration > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    boost::shared_ptr< EmpiricalAccelerationSettings > empiricalSettings =
            boost::dynamic_pointer_cast< EmpiricalAccelerationSettings >(
                accelerationSettings );
    if( empiricalSettings == NULL )
    {
        throw std::runtime_error( "Error, expected empirical acceleration settings when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration );
    }
    else
    {
        // Get pointer to gravity field of central body (for determining keplerian elememts)
        boost::shared_ptr< GravityFieldModel > gravityField = bodyExertingAcceleration->getGravityFieldModel( );

        if( gravityField == NULL )
        {
            throw std::runtime_error( "Error " + nameOfBodyExertingAcceleration + " does not have a gravity field " +
                                      "when making empirical acceleration on" + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            // Create acceleration model.
            accelerationModel = boost::make_shared< EmpiricalAcceleration >(
                        empiricalSettings->constantAcceleration_,
                        empiricalSettings->sineAcceleration_,
                        empiricalSettings->cosineAcceleration_,
                        boost::bind( &Body::getState, bodyUndergoingAcceleration ),
                        boost::bind( &GravityFieldModel::getGravitationalParameter, gravityField ),
                        boost::bind( &Body::getState, bodyExertingAcceleration ) );
        }
    }

    return accelerationModel;
}

//! Function to create a thrust acceleration model.
boost::shared_ptr< propulsion::ThrustAcceleration >
createThrustAcceleratioModel(
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyUndergoingThrust )
{
    // Check input consistency
    boost::shared_ptr< ThrustAccelerationSettings > thrustAccelerationSettings =
            boost::dynamic_pointer_cast< ThrustAccelerationSettings >( accelerationSettings );
    if( thrustAccelerationSettings == NULL )
    {
        throw std::runtime_error( "Error when creating thrust acceleration, input is inconsistent" );
    }

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > magnitudeUpdateSettings;
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > directionUpdateSettings;



    // Check if user-supplied interpolator for full thrust ius present.
    if( thrustAccelerationSettings->interpolatorInterface_ != NULL )
    {
        // Check input consisten
        if( thrustAccelerationSettings->thrustFrame_ == unspecified_thurst_frame )
        {
            throw std::runtime_error( "Error when creating thrust acceleration, input frame is inconsistent with interface" );
        }
        else if( thrustAccelerationSettings->thrustFrame_ != inertial_thurst_frame )
        {
            // Create rotation function from thrust-frame to propagation frame.
            if( thrustAccelerationSettings->thrustFrame_ == lvlh_thrust_frame )
            {
                boost::function< Eigen::Vector6d( ) > vehicleStateFunction =
                        boost::bind( &Body::getState, bodyMap.at( nameOfBodyUndergoingThrust ) );
                boost::function< Eigen::Vector6d( ) > centralBodyStateFunction;

                if( ephemerides::isFrameInertial( thrustAccelerationSettings->centralBody_ ) )
                {
                    centralBodyStateFunction =  boost::lambda::constant( Eigen::Vector6d::Zero( ) );
                }
                else
                {
                    if( bodyMap.count( thrustAccelerationSettings->centralBody_ ) == 0 )
                    {
                        throw std::runtime_error( "Error when creating thrust acceleration, input central body not found" );
                    }
                    centralBodyStateFunction =
                            boost::bind( &Body::getState, bodyMap.at( thrustAccelerationSettings->centralBody_ ) );
                }
                thrustAccelerationSettings->interpolatorInterface_->resetRotationFunction(
                            boost::bind( &reference_frames::getVelocityBasedLvlhToInertialRotationFromFunctions,
                                         vehicleStateFunction, centralBodyStateFunction, true ) );
            }
            else
            {
                throw std::runtime_error( "Error when creating thrust acceleration, input frame not recognized" );
            }
        }
    }

    // Create thrust direction model.
    boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance = createThrustGuidanceModel(
                thrustAccelerationSettings->thrustDirectionGuidanceSettings_, bodyMap, nameOfBodyUndergoingThrust,
                getBodyFixedThrustDirection( thrustAccelerationSettings->thrustMagnitudeSettings_, bodyMap,
                                             nameOfBodyUndergoingThrust ), magnitudeUpdateSettings );

    // Create thrust magnitude model
    boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitude = createThrustMagnitudeWrapper(
                thrustAccelerationSettings->thrustMagnitudeSettings_, bodyMap, nameOfBodyUndergoingThrust,
                directionUpdateSettings );

    // Add required updates of environemt models.
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > totalUpdateSettings;
    propagators::addEnvironmentUpdates( totalUpdateSettings, magnitudeUpdateSettings );
    propagators::addEnvironmentUpdates( totalUpdateSettings, directionUpdateSettings );

    // Set DependentOrientationCalculator for body if required.
    if( !( thrustAccelerationSettings->thrustDirectionGuidanceSettings_->thrustDirectionType_ ==
           thrust_direction_from_existing_body_orientation ) )
    {
        bodyMap.at( nameOfBodyUndergoingThrust )->setDependentOrientationCalculator( thrustDirectionGuidance );
    }

    // Create and return thrust acceleration object.
    boost::function< void( const double ) > updateFunction =
            boost::bind( &updateThrustMagnitudeAndDirection, thrustMagnitude, thrustDirectionGuidance, _1 );
    boost::function< void( const double ) > timeResetFunction =
            boost::bind( &resetThrustMagnitudeAndDirectionTime, thrustMagnitude, thrustDirectionGuidance, _1 );
    return boost::make_shared< propulsion::ThrustAcceleration >(
                boost::bind( &propulsion::ThrustMagnitudeWrapper::getCurrentThrustMagnitude, thrustMagnitude ),
                boost::bind( &propulsion::BodyFixedForceDirectionGuidance ::getCurrentForceDirectionInPropagationFrame, thrustDirectionGuidance ),
                boost::bind( &Body::getBodyMass, bodyMap.at( nameOfBodyUndergoingThrust ) ),
                boost::bind( &propulsion::ThrustMagnitudeWrapper::getCurrentMassRate, thrustMagnitude ),
                thrustAccelerationSettings->thrustMagnitudeSettings_->thrustOriginId_,
                updateFunction, timeResetFunction, totalUpdateSettings );
}

//! Function to create a direct tical acceleration model, according to approach of Lainey et al. (2007, 2009, ...)
boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > createDirectTidalDissipationAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  boost::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Check input consistency
    boost::shared_ptr< DirectTidalDissipationAccelerationSettings > tidalAccelerationSettings =
            boost::dynamic_pointer_cast< DirectTidalDissipationAccelerationSettings >( accelerationSettings );
    if( tidalAccelerationSettings == NULL )
    {
        throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, input is inconsistent" );
    }

    boost::function< double( ) > gravitationalParaterFunctionOfBodyExertingTide;
    boost::function< double( ) > gravitationalParaterFunctionOfBodyUndergoingTide;

    if( tidalAccelerationSettings->useTideRaisedOnPlanet_ )
    {
        if( bodyUndergoingAcceleration->getGravityFieldModel( ) == NULL )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, satellite " +
                                      nameOfBodyUndergoingAcceleration + " has no gravity field" );
        }
        else
        {
            gravitationalParaterFunctionOfBodyUndergoingTide = boost::bind(
                        &GravityFieldModel::getGravitationalParameter, bodyUndergoingAcceleration->getGravityFieldModel( ) );
        }
    }
    else
    {
        if( bodyExertingAcceleration->getGravityFieldModel( ) == NULL )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, satellite " +
                                      nameOfBodyExertingAcceleration + " has no gravity field" );
        }
        else
        {
            gravitationalParaterFunctionOfBodyExertingTide = boost::bind(
                        &GravityFieldModel::getGravitationalParameter, bodyExertingAcceleration->getGravityFieldModel( ) );
        }


        if( bodyUndergoingAcceleration->getGravityFieldModel( ) == NULL )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, satellite " +
                                      nameOfBodyUndergoingAcceleration + " has no gravity field" );
        }
        else
        {
            gravitationalParaterFunctionOfBodyUndergoingTide = boost::bind(
                        &GravityFieldModel::getGravitationalParameter, bodyUndergoingAcceleration->getGravityFieldModel( ) );
        }
    }

    double referenceRadius = TUDAT_NAN;
    if( tidalAccelerationSettings->useTideRaisedOnPlanet_ )
    {
        if( boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) ) == NULL )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, planet " +
                                      nameOfBodyExertingAcceleration + " has no s.h. gravity field" );
        }
        else
        {
            referenceRadius = boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        bodyExertingAcceleration->getGravityFieldModel( ) )->getReferenceRadius( );
        }
    }
    else
    {
        if( boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                    bodyUndergoingAcceleration->getGravityFieldModel( ) ) == NULL )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, planet " +
                                      nameOfBodyUndergoingAcceleration + " has no s.h. gravity field" );
        }
        else
        {
            referenceRadius = boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        bodyUndergoingAcceleration->getGravityFieldModel( ) )->getReferenceRadius( );
        }
    }
    

    if( tidalAccelerationSettings->useTideRaisedOnPlanet_ )
    {
        boost::function< Eigen::Vector3d( ) > planetAngularVelocityVectorFunction =
                boost::bind( &Body::getCurrentAngularVelocityVectorInGlobalFrame, bodyExertingAcceleration );


        return boost::make_shared< DirectTidalDissipationAcceleration >(
                    boost::bind( &Body::getState, bodyUndergoingAcceleration ),
                    boost::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParaterFunctionOfBodyUndergoingTide,
                    planetAngularVelocityVectorFunction,
                    tidalAccelerationSettings->k2LoveNumber_,
                    tidalAccelerationSettings->timeLag_,
                    referenceRadius,
                    tidalAccelerationSettings->includeDirectRadialComponent_);
    }
    else
    {
        return boost::make_shared< DirectTidalDissipationAcceleration >(
                    boost::bind( &Body::getState, bodyUndergoingAcceleration ),
                    boost::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParaterFunctionOfBodyExertingTide,
                    gravitationalParaterFunctionOfBodyUndergoingTide,
                    tidalAccelerationSettings->k2LoveNumber_,
                    tidalAccelerationSettings->timeLag_,
                    referenceRadius,
                    tidalAccelerationSettings->includeDirectRadialComponent_);
    }
}

//! Function to create acceleration model object.
boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > createAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfCentralBody,
        const NamedBodyMap& bodyMap )
{
    // Declare pointer to return object.
    boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;

    // Switch to call correct acceleration model type factory function.
    switch( accelerationSettings->accelerationType_ )
    {
    case central_gravity:
        accelerationModelPointer = createGravitationalAccelerationModel(
                    bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationSettings,
                    nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                    centralBody, nameOfCentralBody );
        break;
    case spherical_harmonic_gravity:
        accelerationModelPointer = createGravitationalAccelerationModel(
                    bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationSettings,
                    nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                    centralBody, nameOfCentralBody );
        break;
    case mutual_spherical_harmonic_gravity:
        accelerationModelPointer = createGravitationalAccelerationModel(
                    bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationSettings,
                    nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                    centralBody, nameOfCentralBody );
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
    case thrust_acceleration:
        accelerationModelPointer = createThrustAcceleratioModel(
                    accelerationSettings, bodyMap,
                    nameOfBodyUndergoingAcceleration );
        break;
    case relativistic_correction_acceleration:
        accelerationModelPointer = createRelativisticCorrectionAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings, bodyMap );
        break;
    case empirical_acceleration:
        accelerationModelPointer = createEmpiricalAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case direct_tidal_dissipation_acceleration:
        accelerationModelPointer = createDirectTidalDissipationAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    default:
        throw std::runtime_error(
                    std::string( "Error, acceleration model ") +
                    std::to_string( accelerationSettings->accelerationType_ ) +
                    " not recognized when making acceleration model of" +
                    nameOfBodyExertingAcceleration + " on " +
                    nameOfBodyUndergoingAcceleration );
        break;
    }
    return accelerationModelPointer;
}

//! Function to put SelectedAccelerationMap in correct order, to ensure correct model creation
SelectedAccelerationList orderSelectedAccelerationMap( const SelectedAccelerationMap& selectedAccelerationsPerBody )
{
    // Declare map of acceleration models acting on current body.
    SelectedAccelerationList orderedAccelerationsPerBody;

    // Iterate over all bodies which are undergoing acceleration
    for( SelectedAccelerationMap::const_iterator bodyIterator =
         selectedAccelerationsPerBody.begin( ); bodyIterator != selectedAccelerationsPerBody.end( );
         bodyIterator++ )
    {
        // Retrieve name of body undergoing acceleration.
        std::string bodyUndergoingAcceleration = bodyIterator->first;

        // Retrieve list of required acceleration model types and bodies exerting accelerationd on
        // current body.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >
                accelerationsForBody = bodyIterator->second;

        // Retrieve indices of all acceleration anf thrust models.
        std::vector< int > aerodynamicAccelerationIndices;
        std::vector< int > thrustAccelerationIndices;

        std::vector< std::pair< std::string, boost::shared_ptr< AccelerationSettings > > >
                currentBodyAccelerations;
        int counter = 0;
        // Iterate over all bodies exerting an acceleration
        for( std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >::
             iterator body2Iterator = accelerationsForBody.begin( );
             body2Iterator != accelerationsForBody.end( ); body2Iterator++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = body2Iterator->first;
            std::vector< boost::shared_ptr< AccelerationSettings > > accelerationList = body2Iterator->second;
            for( unsigned int i = 0; i < accelerationList.size( ); i++ )
            {
                if( accelerationList.at( i )->accelerationType_ == basic_astrodynamics::thrust_acceleration )
                {
                    thrustAccelerationIndices.push_back( counter );
                }
                else if( accelerationList.at( i )->accelerationType_ == basic_astrodynamics::aerodynamic )
                {
                    aerodynamicAccelerationIndices.push_back( counter );
                }
                std::pair< std::string, boost::shared_ptr< AccelerationSettings > >  currentAccelerationPair =
                        std::make_pair( bodyExertingAcceleration, accelerationList.at( i ) );
                currentBodyAccelerations.push_back( currentAccelerationPair );
                counter++;
            }
        }

        if( thrustAccelerationIndices.size( ) > 0 && aerodynamicAccelerationIndices.size( ) > 0 )
        {
            std::vector< int > indexList;
            for( unsigned int i = 0; i < aerodynamicAccelerationIndices.size( ); i++ )
            {
                indexList.push_back( aerodynamicAccelerationIndices.at( i ) );
            }
            for( unsigned int i = 0; i < thrustAccelerationIndices.size( ); i++ )
            {
                indexList.push_back( thrustAccelerationIndices.at( i ) );
            }

            std::vector< int > unorderedIndexList = indexList;
            std::sort( indexList.begin( ), indexList.end( ) );
            if( !( indexList == unorderedIndexList ) )
            {
                std::vector< std::pair< std::string, boost::shared_ptr< AccelerationSettings > > >
                        orderedAccelerationSettings = currentBodyAccelerations;

                int indexCounter = 0;
                for( unsigned int i = 0; i < aerodynamicAccelerationIndices.size( ); i++ )
                {
                    orderedAccelerationSettings[ indexList.at( indexCounter ) ]
                            = currentBodyAccelerations[ aerodynamicAccelerationIndices.at( i ) ];
                    indexCounter++;
                }

                for( unsigned int i = 0; i < thrustAccelerationIndices.size( ); i++ )
                {
                    orderedAccelerationSettings[ indexList.at( indexCounter ) ]
                            = currentBodyAccelerations[ thrustAccelerationIndices.at( i ) ];
                    indexCounter++;
                }

                currentBodyAccelerations = orderedAccelerationSettings;
            }
        }

        orderedAccelerationsPerBody[ bodyUndergoingAcceleration ] = currentBodyAccelerations;
    }

    return orderedAccelerationsPerBody;
}


//! Function to create a set of acceleration models from a map of bodies and acceleration model types.
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::map< std::string, std::string >& centralBodies )
{
    // Declare return map.
    basic_astrodynamics::AccelerationMap accelerationModelMap;

    // Put selectedAccelerationPerBody in correct order
    SelectedAccelerationList orderedAccelerationPerBody =
            orderSelectedAccelerationMap( selectedAccelerationPerBody );

    // Iterate over all bodies which are undergoing acceleration
    for( SelectedAccelerationList::const_iterator bodyIterator =
         orderedAccelerationPerBody.begin( ); bodyIterator != orderedAccelerationPerBody.end( );
         bodyIterator++ )
    {
        boost::shared_ptr< Body > currentCentralBody;

        // Retrieve name of body undergoing acceleration.
        std::string bodyUndergoingAcceleration = bodyIterator->first;

        // Retrieve name of current central body.
        std::string currentCentralBodyName = centralBodies.at( bodyUndergoingAcceleration );

        if( !ephemerides::isFrameInertial( currentCentralBodyName ) )
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
        basic_astrodynamics::SingleBodyAccelerationMap mapOfAccelerationsForBody;

        // Retrieve list of required acceleration model types and bodies exerting accelerationd on
        // current body.
        std::vector< std::pair< std::string, boost::shared_ptr< AccelerationSettings > > >
                accelerationsForBody = bodyIterator->second;

        std::vector< std::pair< std::string, boost::shared_ptr< AccelerationSettings > > > thrustAccelerationSettings;

        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > currentAcceleration;
        // Iterate over all bodies exerting an acceleration
        for( unsigned int i = 0; i < accelerationsForBody.size( ); i++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = accelerationsForBody.at( i ).first;

            // Check if body exerting acceleration is included in bodyMap
            if( bodyMap.count( bodyExertingAcceleration ) ==  0 )
            {
                throw std::runtime_error(
                            std::string( "Error when making acceleration models, requested forces ")
                            + "acting on body " + bodyUndergoingAcceleration  + " due to body " +
                            bodyExertingAcceleration +
                            ", but no such body found in map of bodies" );
            }

            if( !( accelerationsForBody.at( i ).second->accelerationType_ == basic_astrodynamics::thrust_acceleration ) )
            {
                currentAcceleration = createAccelerationModel( bodyMap.at( bodyUndergoingAcceleration ),
                                                               bodyMap.at( bodyExertingAcceleration ),
                                                               accelerationsForBody.at( i ).second,
                                                               bodyUndergoingAcceleration,
                                                               bodyExertingAcceleration,
                                                               currentCentralBody,
                                                               currentCentralBodyName,
                                                               bodyMap );


                // Create acceleration model.
                mapOfAccelerationsForBody[ bodyExertingAcceleration ].push_back(
                            currentAcceleration );
            }
            else
            {
                thrustAccelerationSettings.push_back( accelerationsForBody.at( i ) );
            }

        }

        for( unsigned int i = 0; i < thrustAccelerationSettings.size( ); i++ )
        {
            currentAcceleration = createAccelerationModel( bodyMap.at( bodyUndergoingAcceleration ),
                                                           bodyMap.at( thrustAccelerationSettings.at( i ).first ),
                                                           thrustAccelerationSettings.at( i ).second,
                                                           bodyUndergoingAcceleration,
                                                           thrustAccelerationSettings.at( i ).first,
                                                           currentCentralBody,
                                                           currentCentralBodyName,
                                                           bodyMap );


            // Create acceleration model.
            mapOfAccelerationsForBody[ thrustAccelerationSettings.at( i ).first  ].push_back(
                        currentAcceleration );
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
