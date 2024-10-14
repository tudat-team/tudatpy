/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>
#include <functional>
#include <memory>





#include "tudat/astro/aerodynamics/flightConditions.h"
#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/propulsion/thrustMagnitudeWrapper.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/astro/relativity/relativisticAccelerationCorrection.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/basics/utilities.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/simulation/environment_setup/createOccultationModel.h"
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"


namespace tudat
{

namespace simulation_setup
{

using namespace aerodynamics;
using namespace gravitation;
using namespace basic_astrodynamics;
using namespace electromagnetism;
using namespace ephemerides;


//! Function to create a direct (i.e. not third-body) gravitational acceleration (of any type)
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > createDirectGravitationalAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfCentralBody,
        const bool isCentralBody )
{

    // Check if sum of gravitational parameters (i.e. inertial force w.r.t. central body) should be used.
    bool sumGravitationalParameters = 0;
    if( ( nameOfCentralBody == nameOfBodyExertingAcceleration ) && bodyUndergoingAcceleration != nullptr )
    {
        sumGravitationalParameters = 1;
    }


    // Check type of acceleration model and create.
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel;
    switch( accelerationSettings->accelerationType_ )
    {
    case point_mass_gravity:
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
    case polyhedron_gravity:
        accelerationModel = createPolyhedronGravityAcceleration(
                bodyUndergoingAcceleration,
                bodyExertingAcceleration,
                nameOfBodyUndergoingAcceleration,
                nameOfBodyExertingAcceleration,
                sumGravitationalParameters);
        break;
    case ring_gravity:
        accelerationModel = createRingGravityAcceleration(
                bodyUndergoingAcceleration,
                bodyExertingAcceleration,
                nameOfBodyUndergoingAcceleration,
                nameOfBodyExertingAcceleration,
                sumGravitationalParameters);
        break;
    default:

        std::string errorMessage = "Error when making gravitional acceleration model, cannot parse type " +
                std::to_string( accelerationSettings->accelerationType_ );
        throw std::runtime_error( errorMessage );
    }
    return accelerationModel;
}

//! Function to create a third-body gravitational acceleration (of any type)
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > createThirdBodyGravitationalAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Check type of acceleration model and create.
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel;
    switch( accelerationSettings->accelerationType_ )
    {
    case point_mass_gravity:
        accelerationModel = std::make_shared< ThirdBodyCentralGravityAcceleration >(
                    std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                        createDirectGravitationalAcceleration(
                            bodyUndergoingAcceleration, bodyExertingAcceleration,
                            nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 0 ) ),
                    std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                        createDirectGravitationalAcceleration(
                            centralBody, bodyExertingAcceleration,
                            nameOfCentralBody, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 1 ) ), nameOfCentralBody );
        break;
    case spherical_harmonic_gravity:
        accelerationModel = std::make_shared< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                    std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            bodyUndergoingAcceleration, bodyExertingAcceleration,
                            nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 0 ) ),
                    std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            centralBody, bodyExertingAcceleration, nameOfCentralBody, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 1 ) ), nameOfCentralBody );
        break;
    case mutual_spherical_harmonic_gravity:
        accelerationModel = std::make_shared< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                    std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            bodyUndergoingAcceleration, bodyExertingAcceleration,
                            nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 0 ) ),
                    std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            centralBody, bodyExertingAcceleration, nameOfCentralBody, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 1 ) ), nameOfCentralBody );
        break;
    case polyhedron_gravity:
        accelerationModel = std::make_shared< ThirdBodyPolyhedronGravitationalAccelerationModel >(
                std::dynamic_pointer_cast< PolyhedronGravitationalAccelerationModel >(
                    createDirectGravitationalAcceleration(
                        bodyUndergoingAcceleration, bodyExertingAcceleration,
                        nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                        accelerationSettings, "", 0 ) ),
                std::dynamic_pointer_cast< PolyhedronGravitationalAccelerationModel >(
                    createDirectGravitationalAcceleration(
                        centralBody, bodyExertingAcceleration, nameOfCentralBody, nameOfBodyExertingAcceleration,
                        accelerationSettings, "", 1 ) ),
                nameOfCentralBody );
        break;
    case ring_gravity:
        accelerationModel = std::make_shared< ThirdBodyRingGravitationalAccelerationModel >(
                std::dynamic_pointer_cast< RingGravitationalAccelerationModel >(
                    createDirectGravitationalAcceleration(
                        bodyUndergoingAcceleration, bodyExertingAcceleration,
                        nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                        accelerationSettings, "", 0 ) ),
                std::dynamic_pointer_cast< RingGravitationalAccelerationModel >(
                    createDirectGravitationalAcceleration(
                        centralBody, bodyExertingAcceleration, nameOfCentralBody, nameOfBodyExertingAcceleration,
                        accelerationSettings, "", 1 ) ),
                nameOfCentralBody );
        break;
    default:

        std::string errorMessage = "Error when making third-body gravitional acceleration model, cannot parse type " +
                std::to_string( accelerationSettings->accelerationType_ );
        throw std::runtime_error( errorMessage );
    }
    return accelerationModel;
}

//! Function to create gravitational acceleration (of any type)
std::shared_ptr< AccelerationModel< Eigen::Vector3d > > createGravitationalAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfCentralBody )
{

    std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;
    if( accelerationSettings->accelerationType_ != point_mass_gravity &&
            accelerationSettings->accelerationType_ != spherical_harmonic_gravity &&
            accelerationSettings->accelerationType_ != mutual_spherical_harmonic_gravity &&
            accelerationSettings->accelerationType_ != polyhedron_gravity &&
            accelerationSettings->accelerationType_ != ring_gravity )
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
std::shared_ptr< CentralGravitationalAccelerationModel3d > createCentralGravityAcceleratioModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const bool useCentralBodyFixedFrame )
{
    // Declare pointer to return object.
    std::shared_ptr< CentralGravitationalAccelerationModel3d > accelerationModelPointer;

    // Check if body is endowed with a gravity field model (i.e. is capable of exerting
    // gravitation acceleration).
    if( bodyExertingAcceleration->getGravityFieldModel( ) == nullptr )
    {
        throw std::runtime_error(
                    std::string( "Error, gravity field model not set when making central ") +
                    " gravitational acceleration of " + nameOfBodyExertingAcceleration + " on " +
                    nameOfBodyUndergoingAcceleration );
    }
    else
    {
        std::function< double( ) > gravitationalParameterFunction;

        bool useMutualAttraction = useCentralBodyFixedFrame;
        if( bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr && useMutualAttraction )
        {
            useMutualAttraction = false;
        }

        // Set correct value for gravitational parameter.
        if( !useMutualAttraction )
        {
            gravitationalParameterFunction =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodyExertingAcceleration->getGravityFieldModel( ) );
        }
        else
        {
            std::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodyExertingAcceleration->getGravityFieldModel( ) );
            std::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodyUndergoingAcceleration->getGravityFieldModel( ) );
            gravitationalParameterFunction =
                    std::bind( &utilities::sumFunctionReturn< double >,
                               gravitationalParameterOfBodyExertingAcceleration,
                               gravitationalParameterOfBodyUndergoingAcceleration );
        }

        // Create acceleration object.
        std::function< void( Eigen::Vector3d& ) > bodyUndergoingAccelerationPositionFunction =
                std::bind( &Body::getPositionByReference, bodyUndergoingAcceleration, std::placeholders::_1 );
        std::function< void( Eigen::Vector3d& ) > bodyExertingAccelerationPositionFunction =
                std::bind( &Body::getPositionByReference, bodyExertingAcceleration, std::placeholders::_1 );

        accelerationModelPointer =
                std::make_shared< CentralGravitationalAccelerationModel3d >(
                    bodyUndergoingAccelerationPositionFunction,
                    gravitationalParameterFunction,
                    bodyExertingAccelerationPositionFunction,
                    useMutualAttraction );
    }


    return accelerationModelPointer;
}

//! Function to create spherical harmonic gravity acceleration model.
std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel >
createSphericalHarmonicsGravityAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame,
        const bool useDegreeZeroTerm )
{
    // Declare pointer to return object
    std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicsSettings =
            std::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >(
                accelerationSettings );
    if( sphericalHarmonicsSettings == nullptr )
    {
        throw std::runtime_error(
                    std::string( "Error, acceleration settings inconsistent ") +
                    " making sh gravitational acceleration of " + nameOfBodyExertingAcceleration +
                    " on " + nameOfBodyUndergoingAcceleration );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );

        std::shared_ptr< RotationalEphemeris> rotationalEphemeris =
                bodyExertingAcceleration->getRotationalEphemeris( );
        if( sphericalHarmonicsGravityField == nullptr )
        {
            throw std::runtime_error(
                        std::string( "Error, spherical harmonic gravity field model not set when ")
                        + " making sh gravitational acceleration of " +
                        nameOfBodyExertingAcceleration +
                        " on " + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            if( rotationalEphemeris == nullptr )
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

            std::function< double( ) > gravitationalParameterFunction;

            bool useMutualAttraction = useCentralBodyFixedFrame;
            if( bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr && useMutualAttraction )
            {
                useMutualAttraction = false;
            }

            // Check if mutual acceleration is to be used.
            if( !useMutualAttraction )
            {
                gravitationalParameterFunction =
                        std::bind( &SphericalHarmonicsGravityField::getGravitationalParameter,
                                   sphericalHarmonicsGravityField );
            }
            else
            {
                // Create function returning summed gravitational parameter of the two bodies.
                std::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                   sphericalHarmonicsGravityField );
                std::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                   bodyUndergoingAcceleration->getGravityFieldModel( ) );
                gravitationalParameterFunction =
                        std::bind( &utilities::sumFunctionReturn< double >,
                                   gravitationalParameterOfBodyExertingAcceleration,
                                   gravitationalParameterOfBodyUndergoingAcceleration );
            }

            std::function< Eigen::MatrixXd( ) > originalCosineCoefficientFunction =
                    std::bind( &SphericalHarmonicsGravityField::getCosineCoefficientsBlock,
                               sphericalHarmonicsGravityField,
                               sphericalHarmonicsSettings->maximumDegree_,
                               sphericalHarmonicsSettings->maximumOrder_ );

            std::function< Eigen::MatrixXd( ) > cosineCoefficientFunction;
            if( !useDegreeZeroTerm || sphericalHarmonicsSettings->removePointMass_ )
            {
                cosineCoefficientFunction =
                        std::bind( &setDegreeAndOrderCoefficientToZero, originalCosineCoefficientFunction );
            }
            else
            {
                cosineCoefficientFunction = originalCosineCoefficientFunction;
            }

            // Create acceleration object.
            accelerationModel =
                    std::make_shared< SphericalHarmonicsGravitationalAccelerationModel >(
                    std::bind( &Body::getPositionByReference, bodyUndergoingAcceleration, std::placeholders::_1 ),
                      gravitationalParameterFunction,
                      sphericalHarmonicsGravityField->getReferenceRadius( ),
                      cosineCoefficientFunction,
                      std::bind( &SphericalHarmonicsGravityField::getSineCoefficientsBlock,
                                 sphericalHarmonicsGravityField,
                                 sphericalHarmonicsSettings->maximumDegree_,
                                 sphericalHarmonicsSettings->maximumOrder_ ),
                    std::bind( &Body::getPositionByReference, bodyExertingAcceleration, std::placeholders::_1 ),
                      std::bind( &Body::getCurrentRotationToGlobalFrame,
                                 bodyExertingAcceleration ), useMutualAttraction );
        }
    }
    return accelerationModel;
}

//! Function to create mutual spherical harmonic gravity acceleration model.
std::shared_ptr< gravitation::MutualSphericalHarmonicsGravitationalAccelerationModel >
createMutualSphericalHarmonicsGravityAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame,
        const bool acceleratedBodyIsCentralBody )
{
    using namespace basic_astrodynamics;

    // Declare pointer to return object
    std::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicsSettings =
            std::dynamic_pointer_cast< MutualSphericalHarmonicAccelerationSettings >( accelerationSettings );
    if( mutualSphericalHarmonicsSettings == nullptr )
    {
        std::string errorMessage = "Error, expected mutual spherical harmonics acceleration settings when making acceleration model on " +
                nameOfBodyUndergoingAcceleration + "due to " + nameOfBodyExertingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyExertingAcceleration =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyUndergoingAcceleration->getGravityFieldModel( ) );

        if( sphericalHarmonicsGravityFieldOfBodyExertingAcceleration == nullptr )
        {

            std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );

        }
        else if( sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration == nullptr )
        {

            std::string errorMessage = "Error " + nameOfBodyUndergoingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::function< double( ) > gravitationalParameterFunction;


            // Create function returning summed gravitational parameter of the two bodies.
            if( useCentralBodyFixedFrame == false )
            {
                gravitationalParameterFunction =
                        std::bind( &SphericalHarmonicsGravityField::getGravitationalParameter,
                                   sphericalHarmonicsGravityFieldOfBodyExertingAcceleration );
            }
            else
            {
                // Create function returning summed gravitational parameter of the two bodies.
                std::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                   sphericalHarmonicsGravityFieldOfBodyExertingAcceleration );
                std::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                   sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration );
                gravitationalParameterFunction =
                        std::bind( &utilities::sumFunctionReturn< double >,
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

            accelerationModel = std::make_shared< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        std::bind( &Body::getPositionByReference, bodyUndergoingAcceleration, std::placeholders::_1 ),
                        std::bind( &Body::getPositionByReference, bodyExertingAcceleration, std::placeholders::_1 ),
                        gravitationalParameterFunction,
                        sphericalHarmonicsGravityFieldOfBodyExertingAcceleration->getReferenceRadius( ),
                        sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration->getReferenceRadius( ),
                        std::bind( &SphericalHarmonicsGravityField::getCosineCoefficientsBlock,
                                   sphericalHarmonicsGravityFieldOfBodyExertingAcceleration,
                                   mutualSphericalHarmonicsSettings->maximumDegreeOfBodyExertingAcceleration_,
                                   mutualSphericalHarmonicsSettings->maximumOrderOfBodyExertingAcceleration_ ),
                        std::bind( &SphericalHarmonicsGravityField::getSineCoefficientsBlock,
                                   sphericalHarmonicsGravityFieldOfBodyExertingAcceleration,
                                   mutualSphericalHarmonicsSettings->maximumDegreeOfBodyExertingAcceleration_,
                                   mutualSphericalHarmonicsSettings->maximumOrderOfBodyExertingAcceleration_ ),
                        std::bind( &SphericalHarmonicsGravityField::getCosineCoefficientsBlock,
                                   sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration,
                                   maximumDegreeOfUndergoingBody,
                                   maximumOrderOfUndergoingBody ),
                        std::bind( &SphericalHarmonicsGravityField::getSineCoefficientsBlock,
                                   sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration,
                                   maximumDegreeOfUndergoingBody,
                                   maximumOrderOfUndergoingBody ),
                        std::bind( &Body::getCurrentRotationToGlobalFrame,
                                   bodyExertingAcceleration ),
                        std::bind( &Body::getCurrentRotationToGlobalFrame,
                                   bodyUndergoingAcceleration ),
                        useCentralBodyFixedFrame );
        }
    }
    return accelerationModel;
}

std::shared_ptr< gravitation::PolyhedronGravitationalAccelerationModel >
createPolyhedronGravityAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const bool useCentralBodyFixedFrame)
{

    // Declare pointer to return object
    std::shared_ptr< PolyhedronGravitationalAccelerationModel > accelerationModel;

    // Get pointer to gravity field of central body and cast to required type.
    std::shared_ptr< PolyhedronGravityField > polyhedronGravityField =
            std::dynamic_pointer_cast< PolyhedronGravityField >(
                bodyExertingAcceleration->getGravityFieldModel( ) );

    std::shared_ptr< RotationalEphemeris> rotationalEphemeris = bodyExertingAcceleration->getRotationalEphemeris( );

    if( polyhedronGravityField == nullptr )
    {
        throw std::runtime_error(
                    std::string( "Error, polyhedron gravity field model not set when ")
                    + " making polyhedron gravitational acceleration of " +
                    nameOfBodyExertingAcceleration +
                    " on " + nameOfBodyUndergoingAcceleration );
    }
    else
    {
        if( rotationalEphemeris == nullptr )
        {
            throw std::runtime_error( "Warning when making polyhedron acceleration on body " +
                                      nameOfBodyUndergoingAcceleration + ", no rotation model found for " +
                                      nameOfBodyExertingAcceleration );
        }

        if( rotationalEphemeris->getTargetFrameOrientation( ) != polyhedronGravityField->getFixedReferenceFrame( ) )
        {
            throw std::runtime_error( "Warning when making polyhedron acceleration on body " +
                                      nameOfBodyUndergoingAcceleration + ", rotation model found for " +
                                      nameOfBodyExertingAcceleration + " is incompatible, frames are: " +
                                      rotationalEphemeris->getTargetFrameOrientation( ) + " and " +
                                      polyhedronGravityField->getFixedReferenceFrame( ) );
        }

        std::function< double( ) > gravitationalParameterFunction;

        // Check if mutual acceleration is to be used.
        if( useCentralBodyFixedFrame == false ||
                bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr )
        {
            gravitationalParameterFunction =
                    std::bind( &PolyhedronGravityField::getGravitationalParameter, polyhedronGravityField );
        }
        else
        {
            // Create function returning summed gravitational parameter of the two bodies.
            std::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               polyhedronGravityField );
            std::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodyUndergoingAcceleration->getGravityFieldModel( ) );
            gravitationalParameterFunction =
                    std::bind( &utilities::sumFunctionReturn< double >,
                               gravitationalParameterOfBodyExertingAcceleration,
                               gravitationalParameterOfBodyUndergoingAcceleration );
        }

        std::function< double( ) > volumeFunction =
                std::bind( &PolyhedronGravityField::getVolume, polyhedronGravityField );
        std::function< Eigen::MatrixXd( ) > verticesCoordinatesFunction =
                std::bind( &PolyhedronGravityField::getVerticesCoordinates, polyhedronGravityField );
        std::function< Eigen::MatrixXi( ) > verticesDefiningEachFacetFunction =
                std::bind( &PolyhedronGravityField::getVerticesDefiningEachFacet, polyhedronGravityField );
        std::function< Eigen::MatrixXi( ) > verticesDefiningEachEdgeFunction =
                std::bind( &PolyhedronGravityField::getVerticesDefiningEachEdge, polyhedronGravityField );
        std::function< std::vector< Eigen::MatrixXd >( ) > facetDyadsFunction =
                std::bind( &PolyhedronGravityField::getFacetDyads, polyhedronGravityField );
        std::function< std::vector< Eigen::MatrixXd >( ) > edgeDyadsFunction =
                std::bind( &PolyhedronGravityField::getEdgeDyads, polyhedronGravityField );

        // Create acceleration object.
        accelerationModel =
                std::make_shared< PolyhedronGravitationalAccelerationModel >(
                        std::bind( &Body::getPositionByReference, bodyUndergoingAcceleration, std::placeholders::_1 ),
                        gravitationalParameterFunction,
                        volumeFunction,
                        verticesCoordinatesFunction,
                        verticesDefiningEachFacetFunction,
                        verticesDefiningEachEdgeFunction,
                        facetDyadsFunction,
                        edgeDyadsFunction,
                        std::bind( &Body::getPositionByReference, bodyExertingAcceleration, std::placeholders::_1 ),
                        std::bind( &Body::getCurrentRotationToGlobalFrame, bodyExertingAcceleration ),
                        useCentralBodyFixedFrame );

    }
    return accelerationModel;
}

std::shared_ptr< gravitation::RingGravitationalAccelerationModel > createRingGravityAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const bool useCentralBodyFixedFrame)
{

    // Declare pointer to return object
    std::shared_ptr< RingGravitationalAccelerationModel > accelerationModel;

    // Get pointer to gravity field of central body and cast to required type.
    std::shared_ptr< RingGravityField > ringGravityField =
            std::dynamic_pointer_cast< RingGravityField >(
                bodyExertingAcceleration->getGravityFieldModel( ) );

    std::shared_ptr< RotationalEphemeris> rotationalEphemeris = bodyExertingAcceleration->getRotationalEphemeris( );

    if( ringGravityField == nullptr )
    {
        throw std::runtime_error(
                    std::string( "Error, ring gravity field model not set when ")
                    + " making ring gravitational acceleration of " +
                    nameOfBodyExertingAcceleration +
                    " on " + nameOfBodyUndergoingAcceleration );
    }
    else
    {
        if( rotationalEphemeris == nullptr )
        {
            throw std::runtime_error( "Warning when making ring acceleration on body " +
                                      nameOfBodyUndergoingAcceleration + ", no rotation model found for " +
                                      nameOfBodyExertingAcceleration );
        }

        if( rotationalEphemeris->getTargetFrameOrientation( ) != ringGravityField->getFixedReferenceFrame( ) )
        {
            throw std::runtime_error( "Warning when making ring acceleration on body " +
                                      nameOfBodyUndergoingAcceleration + ", rotation model found for " +
                                      nameOfBodyExertingAcceleration + " is incompatible, frames are: " +
                                      rotationalEphemeris->getTargetFrameOrientation( ) + " and " +
                                      ringGravityField->getFixedReferenceFrame( ) );
        }

        std::function< double( ) > gravitationalParameterFunction;

        // Check if mutual acceleration is to be used.
        if( useCentralBodyFixedFrame == false ||
                bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr )
        {
            gravitationalParameterFunction =
                    std::bind( &RingGravityField::getGravitationalParameter, ringGravityField );
        }
        else
        {
            // Create function returning summed gravitational parameter of the two bodies.
            std::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               ringGravityField );
            std::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodyUndergoingAcceleration->getGravityFieldModel( ) );
            gravitationalParameterFunction =
                    std::bind( &utilities::sumFunctionReturn< double >,
                               gravitationalParameterOfBodyExertingAcceleration,
                               gravitationalParameterOfBodyUndergoingAcceleration );
        }

        std::function< double( ) > ringRadiusFunction = std::bind( &RingGravityField::getRingRadius, ringGravityField );

        // Create acceleration object.
        accelerationModel =
                std::make_shared< RingGravitationalAccelerationModel >(
                        std::bind( &Body::getPositionByReference, bodyUndergoingAcceleration, std::placeholders::_1 ),
                        gravitationalParameterFunction,
                        ringRadiusFunction,
                        ringGravityField->getEllipticIntegralSFromDAndB( ),
                        std::bind( &Body::getPositionByReference, bodyExertingAcceleration, std::placeholders::_1 ),
                        std::bind( &Body::getCurrentRotationToGlobalFrame, bodyExertingAcceleration ),
                        useCentralBodyFixedFrame );

    }
    return accelerationModel;
}

//! Function to create a third body central gravity acceleration model.
std::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration >
createThirdBodyCentralGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody )
{
    // Declare pointer to return object.
    std::shared_ptr< ThirdBodyCentralGravityAcceleration > accelerationModelPointer;

    // Create acceleration object.
    accelerationModelPointer =  std::make_shared< ThirdBodyCentralGravityAcceleration >(
                std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createCentralGravityAcceleratioModel( bodyUndergoingAcceleration,
                                                          bodyExertingAcceleration,
                                                          nameOfBodyUndergoingAcceleration,
                                                          nameOfBodyExertingAcceleration, 0 ) ),
                std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createCentralGravityAcceleratioModel( centralBody, bodyExertingAcceleration,
                                                          nameOfCentralBody,
                                                          nameOfBodyExertingAcceleration, 0 ) ), nameOfCentralBody );

    return accelerationModelPointer;
}

//! Function to create a third body spheric harmonic gravity acceleration model.
std::shared_ptr< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >
createThirdBodySphericalHarmonicGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    using namespace basic_astrodynamics;

    // Declare pointer to return object
    std::shared_ptr< ThirdBodySphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicsSettings =
            std::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >( accelerationSettings );
    if( sphericalHarmonicsSettings == nullptr )
    {
        std::string errorMessage = "Error, expected spherical harmonics acceleration settings when making acceleration model on " +
                nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        if( sphericalHarmonicsGravityField == nullptr )
        {
            std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making third body spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else
        {

            accelerationModel =  std::make_shared< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                        std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                            createSphericalHarmonicsGravityAcceleration(
                                bodyUndergoingAcceleration, bodyExertingAcceleration, nameOfBodyUndergoingAcceleration,
                                nameOfBodyExertingAcceleration, sphericalHarmonicsSettings, 0 ) ),
                        std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                            createSphericalHarmonicsGravityAcceleration(
                                centralBody, bodyExertingAcceleration, nameOfCentralBody,
                                nameOfBodyExertingAcceleration, sphericalHarmonicsSettings, 0 ) ), nameOfCentralBody );
        }
    }
    return accelerationModel;
}

//! Function to create a third body mutual spheric harmonic gravity acceleration model.
std::shared_ptr< gravitation::ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >
createThirdBodyMutualSphericalHarmonicGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Declare pointer to return object
    std::shared_ptr< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicsSettings =
            std::dynamic_pointer_cast< MutualSphericalHarmonicAccelerationSettings >( accelerationSettings );
    if( mutualSphericalHarmonicsSettings == nullptr )
    {

        std::string errorMessage = "Error, expected mutual spherical harmonics acceleration settings when making acceleration model on " +
                nameOfBodyUndergoingAcceleration +
                " due to " + nameOfBodyExertingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyExertingAcceleration =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyUndergoingAcceleration->getGravityFieldModel( ) );
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfCentralBody =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    centralBody->getGravityFieldModel( ) );

        if( sphericalHarmonicsGravityFieldOfBodyExertingAcceleration == nullptr )
        {
            std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else if( sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration == nullptr )
        {
            std::string errorMessage = "Error " + nameOfBodyUndergoingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else if( sphericalHarmonicsGravityFieldOfCentralBody == nullptr )
        {
            std::string errorMessage = "Error " + nameOfCentralBody + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< MutualSphericalHarmonicAccelerationSettings > accelerationSettingsForCentralBodyAcceleration =
                    std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                        mutualSphericalHarmonicsSettings->maximumDegreeOfBodyExertingAcceleration_,
                        mutualSphericalHarmonicsSettings->maximumOrderOfBodyExertingAcceleration_,
                        mutualSphericalHarmonicsSettings->maximumDegreeOfCentralBody_,
                        mutualSphericalHarmonicsSettings->maximumOrderOfCentralBody_ );
            accelerationModel =  std::make_shared< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                        std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                            createMutualSphericalHarmonicsGravityAcceleration(
                                bodyUndergoingAcceleration, bodyExertingAcceleration, nameOfBodyUndergoingAcceleration,
                                nameOfBodyExertingAcceleration, mutualSphericalHarmonicsSettings, 0, 0 ) ),
                        std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                            createMutualSphericalHarmonicsGravityAcceleration(
                                centralBody, bodyExertingAcceleration, nameOfCentralBody,
                                nameOfBodyExertingAcceleration, accelerationSettingsForCentralBodyAcceleration, 0, 1 ) ),
                        nameOfCentralBody );
        }
    }
    return accelerationModel;
}

std::shared_ptr< gravitation::ThirdBodyPolyhedronGravitationalAccelerationModel >
createThirdBodyPolyhedronGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody )
{
    using namespace basic_astrodynamics;

    // Declare pointer to return object
    std::shared_ptr< ThirdBodyPolyhedronGravitationalAccelerationModel > accelerationModel;

    // Get pointer to gravity field of central body and cast to required type.
    std::shared_ptr< PolyhedronGravityField > polyhedronGravityField =
            std::dynamic_pointer_cast< PolyhedronGravityField >( bodyExertingAcceleration->getGravityFieldModel( ) );
    if( polyhedronGravityField == nullptr )
    {
        std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a polyhedron gravity field " +
                "when making third body polyhedron gravity acceleration on " +
                nameOfBodyUndergoingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {

        accelerationModel =  std::make_shared< ThirdBodyPolyhedronGravitationalAccelerationModel >(
            std::dynamic_pointer_cast< PolyhedronGravitationalAccelerationModel >(
                createPolyhedronGravityAcceleration(
                    bodyUndergoingAcceleration, bodyExertingAcceleration, nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration, 0 ) ),
            std::dynamic_pointer_cast< PolyhedronGravitationalAccelerationModel >(
                createPolyhedronGravityAcceleration(
                    centralBody, bodyExertingAcceleration, nameOfCentralBody,
                    nameOfBodyExertingAcceleration, 0 ) ),
            nameOfCentralBody );
    }

    return accelerationModel;
}

std::shared_ptr< gravitation::ThirdBodyRingGravitationalAccelerationModel > createThirdBodyRingGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody )
{
    using namespace basic_astrodynamics;

    // Declare pointer to return object
    std::shared_ptr< ThirdBodyRingGravitationalAccelerationModel > accelerationModel;

    // Get pointer to gravity field of central body and cast to required type.
    std::shared_ptr< RingGravityField > ringGravityField =
            std::dynamic_pointer_cast< RingGravityField >( bodyExertingAcceleration->getGravityFieldModel( ) );
    if( ringGravityField == nullptr )
    {
        std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a ring gravity field " +
                "when making third body ring gravity acceleration on " +
                nameOfBodyUndergoingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        accelerationModel =  std::make_shared< ThirdBodyRingGravitationalAccelerationModel >(
            std::dynamic_pointer_cast< RingGravitationalAccelerationModel >(
                createPolyhedronGravityAcceleration(
                    bodyUndergoingAcceleration, bodyExertingAcceleration, nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration, 0 ) ),
            std::dynamic_pointer_cast< RingGravitationalAccelerationModel >(
                createPolyhedronGravityAcceleration(
                    centralBody, bodyExertingAcceleration, nameOfCentralBody,
                    nameOfBodyExertingAcceleration, 0 ) ),
            nameOfCentralBody );
    }

    return accelerationModel;
}

//! Function to create an aerodynamic acceleration model.
std::shared_ptr< aerodynamics::AerodynamicAcceleration > createAerodynamicAcceleratioModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    // Check existence of required environment models
    if( bodyUndergoingAcceleration->getAerodynamicCoefficientInterface( ) == nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, body " +
                                  nameOfBodyUndergoingAcceleration +
                                  "has no aerodynamic coefficients." );
    }

    if( bodyExertingAcceleration->getAtmosphereModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, central body " +
                                  nameOfBodyExertingAcceleration + " has no atmosphere model.");
    }

    if( bodyExertingAcceleration->getShapeModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, central body " +
                                  nameOfBodyExertingAcceleration + " has no shape model." );
    }

    // Retrieve flight conditions; create object if not yet extant.
    std::shared_ptr< AtmosphericFlightConditions > bodyFlightConditions =
            std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodyUndergoingAcceleration->getFlightConditions( ) );

    if( bodyFlightConditions == nullptr && bodyUndergoingAcceleration->getFlightConditions( ) == nullptr )
    {

        bodyFlightConditions = createAtmosphericFlightConditions( bodyUndergoingAcceleration,
                                                                  bodyExertingAcceleration,
                                                                  nameOfBodyUndergoingAcceleration,
                                                                  nameOfBodyExertingAcceleration );
        bodyUndergoingAcceleration->setFlightConditions( bodyFlightConditions );
    }
    else if( bodyFlightConditions == nullptr && bodyUndergoingAcceleration->getFlightConditions( ) != nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, found flight conditions that are not atmospheric." );
    }

    // Retrieve frame in which aerodynamic coefficients are defined.
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > aerodynamicCoefficients =
            bodyUndergoingAcceleration->getAerodynamicCoefficientInterface( );
    reference_frames::AerodynamicsReferenceFrames accelerationFrame = aerodynamics::getCompleteFrameForCoefficients(
            aerodynamicCoefficients->getForceCoefficientsFrame( ) );

    // Create function to transform from frame of aerodynamic coefficienrs to that of propagation.
    std::function< void( Eigen::Vector3d&, const Eigen::Vector3d& ) > toPropagationFrameTransformation;
    toPropagationFrameTransformation =
            reference_frames::getAerodynamicForceTransformationReferenceFunction(
                bodyFlightConditions->getAerodynamicAngleCalculator( ),
                accelerationFrame,
                std::bind( &Body::getCurrentRotationToGlobalFrameReference, bodyExertingAcceleration ),
                reference_frames::inertial_frame );

    std::function< Eigen::Vector3d&( ) > coefficientFunction =
            std::bind( &AerodynamicCoefficientInterface::getCurrentForceCoefficientsReference,
                       aerodynamicCoefficients );
    std::function< void( Eigen::Vector3d& ) > coefficientInPropagationFrameFunction =
            std::bind( &reference_frames::transformVectorFunctionFromVectorReferenceFunctions,
                       std::placeholders::_1,
                       coefficientFunction, toPropagationFrameTransformation );

    // Create acceleration model.
    return std::make_shared< AerodynamicAcceleration >(
                coefficientInPropagationFrameFunction,
                std::bind( &AtmosphericFlightConditions::getCurrentDensity, bodyFlightConditions ),
                std::bind( &AtmosphericFlightConditions::getCurrentAirspeed, bodyFlightConditions ),
                std::bind( &Body::getBodyMass, bodyUndergoingAcceleration ),
                std::bind( &AerodynamicCoefficientInterface::getReferenceArea,
                           aerodynamicCoefficients ),
                aerodynamics::areCoefficientsInNegativeDirection(
                        aerodynamicCoefficients->getForceCoefficientsFrame( ) ) );
}

std::shared_ptr< RadiationPressureAcceleration >
createRadiationPressureAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const SystemOfBodies& bodies,
        const std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    RadiationPressureTargetModelType targetModelType = undefined_target;
    if( std::dynamic_pointer_cast< RadiationPressureAccelerationSettings >( accelerationSettings  ) != nullptr )
    {
        targetModelType = std::dynamic_pointer_cast< RadiationPressureAccelerationSettings >( accelerationSettings  )->targetModelType_;
    }
    if( targetModelType == multi_type_target )
    {
        throw std::runtime_error( "Error when creating radiation pressure acceleration, cannot select multi-type target" );
    }

    // Create references to bodies in radiation pressure terms
    const auto& sourceName = nameOfBodyExertingAcceleration;
    const auto& targetName = nameOfBodyUndergoingAcceleration;
    const auto& source = bodyExertingAcceleration;
    const auto& target = bodyUndergoingAcceleration;

    if( source->getRadiationSourceModel() == nullptr )
    {
        throw std::runtime_error( "Error when making radiation pressure acceleration, body " +
                                  sourceName +
                                  " has no radiation source model." );
    }

    std::shared_ptr<electromagnetism::RadiationPressureTargetModel> targetModel;
    if( target->getRadiationPressureTargetModels().size( ) == 0 )
    {
        throw std::runtime_error( "Error when making radiation pressure acceleration, body " +
                                  targetName +
                                  " has no radiation pressure target models." );
    }
    else
    {
        targetModel = getRadiationPressureTargetModelOfType( target, targetModelType, " when making radiation pressure acceleration due to " + sourceName + " " );
    }

    // Cast source and target models for type checks
    auto isotropicPointRadiationSourceModel =
            std::dynamic_pointer_cast<electromagnetism::IsotropicPointRadiationSourceModel>(
                    source->getRadiationSourceModel());
    auto paneledRadiationSourceModel =
            std::dynamic_pointer_cast<electromagnetism::PaneledRadiationSourceModel>(
                    source->getRadiationSourceModel());
    auto cannonballRadiationPressureTargetModel =
            std::dynamic_pointer_cast<electromagnetism::CannonballRadiationPressureTargetModel>(
                    targetModel);

    // Get target rotation function
    std::function<Eigen::Quaterniond()> targetRotationFromLocalToGlobalFrameFunction;
    if (cannonballRadiationPressureTargetModel != nullptr) {
        // Cannonball target is rotation-invariant and can use identity rotation
        // Enables use of cannonball target without target body rotation model
        // TODO-DOMINIK maybe move to RP acceleration subclass so that cannonball does not need target rotation, then move into if below
        targetRotationFromLocalToGlobalFrameFunction = [] () { return Eigen::Quaterniond::Identity(); };
    }
    else
    {
        targetRotationFromLocalToGlobalFrameFunction = [target] { return target->getCurrentRotationToGlobalFrame(); };
    }

    // Find occulting bodies corresponding to this source
    auto targetOccultingBodiesMap = targetModel->getSourceToTargetOccultingBodies();
    std::vector<std::string> sourceToTargetOccultingBodies = {};
    if (targetOccultingBodiesMap.count(sourceName) > 0)
    {
        sourceToTargetOccultingBodies = targetOccultingBodiesMap.at(sourceName);
    }
    else if (targetOccultingBodiesMap.count("") > 0)
    {
        // Use the same occulting bodies for all sources of this target
        sourceToTargetOccultingBodies = targetOccultingBodiesMap.at("");
    }

    // Check if occulting bodies are not source or target
    for (auto& occultingBodyName : sourceToTargetOccultingBodies)
    {
        if (occultingBodyName == sourceName)
        {
            throw std::runtime_error( "Error when making radiation pressure acceleration, source body cannot "
                                      "act as occulting body.");
        }
        if (occultingBodyName == targetName)
        {
            throw std::runtime_error( "Error when making radiation pressure acceleration, target body cannot "
                                      "act as occulting body.");
        }
    }
    auto sourceToTargetOccultationModel = createOccultationModel(sourceToTargetOccultingBodies, bodies);
    
    // Create acceleration model
    if (isotropicPointRadiationSourceModel != nullptr)
    {
        return std::make_shared<IsotropicPointSourceRadiationPressureAcceleration>(
                isotropicPointRadiationSourceModel,
                source->getShapeModel(),
                [source] { return source->getPosition(); },
                targetModel,
                [target] { return target->getPosition(); },
                targetRotationFromLocalToGlobalFrameFunction,
                [target] { return target->getBodyMass(); },
                sourceToTargetOccultationModel);
    }
    else if (paneledRadiationSourceModel != nullptr)
    {
        return std::make_shared<PaneledSourceRadiationPressureAcceleration>(
                paneledRadiationSourceModel,
                [source] { return source->getPosition(); },
                [source] { return source->getCurrentRotationToGlobalFrame(); },
                targetModel,
                [target] { return target->getPosition(); },
                targetRotationFromLocalToGlobalFrameFunction,
                [target] { return target->getBodyMass(); },
                sourceToTargetOccultationModel);
    }
    else
    {
        throw std::runtime_error( "Error when making radiation pressure acceleration, radiation source model "
                                  "for body " + sourceName +
                                  " is not supported." );
    }
}

//! Function to create a cannonball radiation pressure acceleration model.
// RP-OLD
std::shared_ptr< AccelerationModel3d >
createCannonballRadiationPressureAcceleratioModel(
    const std::shared_ptr< Body > bodyUndergoingAcceleration,
    const std::shared_ptr< Body > bodyExertingAcceleration,
    const std::string& nameOfBodyUndergoingAcceleration,
    const std::string& nameOfBodyExertingAcceleration,
    const SystemOfBodies& bodies )
{
    std::cerr<<"Warning, you are using the deprecated (as of tudatpy v0.8) version of the cannonball radiation pressure model"<<
               ", the interface you are using will be dropped from v0.9 onwards. To learn how to convert your code to the new interfaces"<<
               ", and be able to use the powerful new radiation pressure framework, see "<<
               "https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/radiation_pressure_acceleration.html#backwards-compatibility"<<std::endl;
    // Retrieve radiation pressure interface
    if( bodyUndergoingAcceleration->getRadiationPressureInterfaces( ).count(
        nameOfBodyExertingAcceleration ) == 0 )
    {
        throw std::runtime_error(
            "Error when making radiation pressure, no radiation pressure interface found  in " +
            nameOfBodyUndergoingAcceleration +
            " for body " + nameOfBodyExertingAcceleration );
    }
    std::shared_ptr< RadiationPressureInterface > radiationPressureInterface =
        bodyUndergoingAcceleration->getRadiationPressureInterfaces( ).at(
            nameOfBodyExertingAcceleration );

    std::map<std::string, std::vector<std::string> > occultingBodies;
    occultingBodies[ nameOfBodyExertingAcceleration ] = radiationPressureInterface->getOccultingBodies( );
    std::shared_ptr< CannonballRadiationPressureTargetModel > cannonBallTargetModel =
        std::make_shared< CannonballRadiationPressureTargetModel >(
            radiationPressureInterface->getArea( ),
            radiationPressureInterface->getRadiationPressureCoefficient( ),
            occultingBodies );

    std::vector< std::shared_ptr<electromagnetism::RadiationPressureTargetModel> > targetModels = bodyUndergoingAcceleration->getRadiationPressureTargetModels( );
    bool addModel = false;
    if(  targetModels.size( ) == 0 )
    {
        addModel = true;
    }
    else if( targetModels.size( ) == 1  )
    {
        std::shared_ptr< CannonballRadiationPressureTargetModel > existingCannonBallTargetModel =
            std::dynamic_pointer_cast< CannonballRadiationPressureTargetModel >( targetModels.at( 0 ) );
        if( existingCannonBallTargetModel == nullptr )
        {
            throw std::runtime_error( "Error when create deprecated cannonball radiation pressure model; existing target model found of non-cannonball type." );
        }
        else if( !( cannonBallTargetModel->getArea( ) == existingCannonBallTargetModel->getArea( ) &&
            cannonBallTargetModel->getCoefficient( ) == existingCannonBallTargetModel->getCoefficient( ) &&
            cannonBallTargetModel->getSourceToTargetOccultingBodies( ) == existingCannonBallTargetModel->getSourceToTargetOccultingBodies( ) ) )
        {
            throw std::runtime_error( "Error when create deprecated cannonball radiation pressure model; existing target model found of cannonball type with inconsistent settings." );
        }
    }
    else
    {
        throw std::runtime_error( "Error when create deprecated cannonball radiation pressure model; deprecated interface does not permit multiple existing target models." );
    }

    if( addModel )
    {
        bodyUndergoingAcceleration->addRadiationPressureTargetModel( cannonBallTargetModel );
    }

    return createRadiationPressureAccelerationModel(
        bodyUndergoingAcceleration, bodyExertingAcceleration, nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
        bodies );
}

//! Function to create Yarkovsky acceleration model, based on the simplified model of https://doi.org/10.1038/s43247-021-00337-x.
std::shared_ptr< electromagnetism::YarkovskyAcceleration > createYarkovskyAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  std::shared_ptr< AccelerationSettings > accelerationSettings )
{

    // Declare pointer to return object
    std::shared_ptr< electromagnetism::YarkovskyAcceleration > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< simulation_setup::YarkovskyAccelerationSettings > yarkovskySettings =
            std::dynamic_pointer_cast< simulation_setup::YarkovskyAccelerationSettings >(
                accelerationSettings );

    if( yarkovskySettings == nullptr )
    {
        throw std::runtime_error( "Error, expected Yarkovsky acceleration settings when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration );
    }
    else
    {
        accelerationModel = std::make_shared< electromagnetism::YarkovskyAcceleration >(
                    yarkovskySettings->yarkovskyParameter_,
                    std::bind( &Body::getState, bodyUndergoingAcceleration ),
                    std::bind( &Body::getState, bodyExertingAcceleration ) );
        // }
    }

    return accelerationModel;
}


std::shared_ptr< basic_astrodynamics::CustomAccelerationModel > createCustomAccelerationModel(
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration )
{
    std::shared_ptr< CustomAccelerationSettings > customAccelerationSettings =
            std::dynamic_pointer_cast< CustomAccelerationSettings >(
                accelerationSettings );
    if( customAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error, expected custom acceleration settings when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration  );
    }
    return std::make_shared< CustomAccelerationModel >( customAccelerationSettings->accelerationFunction_ );

}

//! Function to create an orbiter relativistic correction acceleration model
std::shared_ptr< relativity::RelativisticAccelerationCorrection > createRelativisticCorrectionAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const SystemOfBodies& bodies )
{
    using namespace relativity;

    // Declare pointer to return object
    std::shared_ptr< RelativisticAccelerationCorrection > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< RelativisticAccelerationCorrectionSettings > relativisticAccelerationSettings =
            std::dynamic_pointer_cast< RelativisticAccelerationCorrectionSettings >(
                accelerationSettings );
    if( relativisticAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error, expected relativistic acceleration settings when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration );
    }
    else
    {

        // Retrieve function pointers for properties of bodies exerting/undergoing acceleration.
        std::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingAcceleration =
                std::bind( &Body::getState, bodyExertingAcceleration );
        std::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingAcceleration =
                std::bind( &Body::getState, bodyUndergoingAcceleration );

        std::function< double( ) > centralBodyGravitationalParameterFunction;
        std::shared_ptr< GravityFieldModel > gravityField = bodyExertingAcceleration->getGravityFieldModel( );
        if( gravityField == nullptr )
        {
            throw std::runtime_error( "Error " + nameOfBodyExertingAcceleration + " does not have a gravity field " +
                                      "when making relativistic acceleration on" + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            centralBodyGravitationalParameterFunction =
                    std::bind( &GravityFieldModel::getGravitationalParameter, bodyExertingAcceleration->getGravityFieldModel( ) );
        }

        // Create acceleration model if only schwarzschild term is to be used.
        if( relativisticAccelerationSettings->calculateLenseThirringCorrection_ == false &&
                relativisticAccelerationSettings->calculateDeSitterCorrection_ == false )
        {
            std::function< double( ) > ppnGammaFunction = std::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet );
            std::function< double( ) > ppnBetaFunction = std::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet );

            // Create acceleration model.
            accelerationModel = std::make_shared< RelativisticAccelerationCorrection >
                    ( stateFunctionOfBodyUndergoingAcceleration,
                      stateFunctionOfBodyExertingAcceleration,
                      centralBodyGravitationalParameterFunction,
                      ppnGammaFunction, ppnBetaFunction );

        }
        else
        {

            // Retrieve parameters of primary body if de Sitter term is to be used.
            std::function< Eigen::Vector6d( ) > stateFunctionOfPrimaryBody;
            std::function< double( ) > primaryBodyGravitationalParameterFunction;
            if( relativisticAccelerationSettings->calculateDeSitterCorrection_ == true )
            {
                if(  bodies.count( relativisticAccelerationSettings->primaryBody_ ) == 0 )
                {
                    throw std::runtime_error( "Error, no primary body " + relativisticAccelerationSettings->primaryBody_ +
                                              " found when making de Sitter acceleration correction" );
                }
                stateFunctionOfPrimaryBody =
                        std::bind( &Body::getState, bodies.at( relativisticAccelerationSettings->primaryBody_ ) );

                if(  bodies.at( relativisticAccelerationSettings->primaryBody_ )->getGravityFieldModel( ) == nullptr )
                {
                    throw std::runtime_error( "Error, primary body " + relativisticAccelerationSettings->primaryBody_ +
                                              " has no gravity field when making de Sitter acceleration correction" );
                }

                primaryBodyGravitationalParameterFunction =
                        std::bind( &GravityFieldModel::getGravitationalParameter,
                                   bodies.at( relativisticAccelerationSettings->primaryBody_ )->getGravityFieldModel( ) );


            }

            // Retrieve angular momentum vector if Lense-Thirring
            std::function< Eigen::Vector3d( ) > angularMomentumFunction;
            if( relativisticAccelerationSettings->calculateLenseThirringCorrection_ == true  )
            {
                angularMomentumFunction = [ = ]( ){ return
                            relativisticAccelerationSettings->centralBodyAngularMomentum_; };
            }

            if( relativisticAccelerationSettings->calculateDeSitterCorrection_ == true )
            {
                // Create acceleration model with Lense-Thirring and de Sitter terms.
                accelerationModel = std::make_shared< RelativisticAccelerationCorrection >
                        ( stateFunctionOfBodyUndergoingAcceleration,
                          stateFunctionOfBodyExertingAcceleration,
                          stateFunctionOfPrimaryBody,
                          centralBodyGravitationalParameterFunction,
                          primaryBodyGravitationalParameterFunction,
                          relativisticAccelerationSettings->primaryBody_,
                          angularMomentumFunction,
                          std::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet ),
                          std::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet ),
                          relativisticAccelerationSettings->calculateSchwarzschildCorrection_ );
            }
            else
            {
                // Create acceleration model with Lense-Thirring and term.
                accelerationModel = std::make_shared< RelativisticAccelerationCorrection >
                        ( stateFunctionOfBodyUndergoingAcceleration,
                          stateFunctionOfBodyExertingAcceleration,
                          centralBodyGravitationalParameterFunction,
                          angularMomentumFunction,
                          std::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet ),
                          std::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet ),
                          relativisticAccelerationSettings->calculateSchwarzschildCorrection_ );
            }
        }
    }
    return accelerationModel;
}


//! Function to create empirical acceleration model.
std::shared_ptr< EmpiricalAcceleration > createEmpiricalAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Declare pointer to return object
    std::shared_ptr< EmpiricalAcceleration > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< EmpiricalAccelerationSettings > empiricalSettings =
            std::dynamic_pointer_cast< EmpiricalAccelerationSettings >(
                accelerationSettings );
    if( empiricalSettings == nullptr )
    {
        throw std::runtime_error( "Error, expected empirical acceleration settings when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration );
    }
    else
    {
        // Get pointer to gravity field of central body (for determining keplerian elememts)
        std::shared_ptr< GravityFieldModel > gravityField = bodyExertingAcceleration->getGravityFieldModel( );

        if( gravityField == nullptr )
        {
            throw std::runtime_error( "Error " + nameOfBodyExertingAcceleration + " does not have a gravity field " +
                                      "when making empirical acceleration on" + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            // Create acceleration model.
            accelerationModel = std::make_shared< EmpiricalAcceleration >(
                        empiricalSettings->constantAcceleration_,
                        empiricalSettings->sineAcceleration_,
                        empiricalSettings->cosineAcceleration_,
                        std::bind( &Body::getState, bodyUndergoingAcceleration ),
                        std::bind( &GravityFieldModel::getGravitationalParameter, gravityField ),
                        std::bind( &Body::getState, bodyExertingAcceleration ) );
        }
    }

    return accelerationModel;
}

//! Function to create a thrust acceleration model.
std::shared_ptr< propulsion::ThrustAcceleration >
createThrustAcceleratioModel(
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const SystemOfBodies& bodies,
        const std::string& nameOfBodyUndergoingThrust )
{
    // Check input consistency
    std::shared_ptr< ThrustAccelerationSettings > thrustAccelerationSettings =
            std::dynamic_pointer_cast< ThrustAccelerationSettings >( accelerationSettings );
    if( thrustAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating thrust acceleration, input is inconsistent" );
    }

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > magnitudeUpdateSettings;
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > directionUpdateSettings;


    // Create thrust magnitude model
    std::shared_ptr< propulsion::ThrustDirectionCalculator > thrustDirectionWrapper;
    if( std::dynamic_pointer_cast< ephemerides::DirectionBasedRotationalEphemeris >(
                bodies.at( nameOfBodyUndergoingThrust )->getRotationalEphemeris( ) ) != nullptr )
    {
        thrustDirectionWrapper = std::make_shared< propulsion::DirectThrustDirectionCalculator >(
                    std::dynamic_pointer_cast< ephemerides::DirectionBasedRotationalEphemeris >(
                    bodies.at( nameOfBodyUndergoingThrust )->getRotationalEphemeris( ) ) );
    }
    else
    {
        thrustDirectionWrapper = std::make_shared< propulsion::OrientationBasedThrustDirectionCalculator >(
                    std::bind( &Body::getCurrentRotationToGlobalFrame, bodies.at( nameOfBodyUndergoingThrust ) ) );
        directionUpdateSettings[ propagators::body_rotational_state_update ].push_back( nameOfBodyUndergoingThrust );
    }

    std::vector< std::shared_ptr< system_models::EngineModel > > engineModelsForAcceleration;
    std::shared_ptr< system_models::VehicleSystems > vehicleSystems = bodies.at( nameOfBodyUndergoingThrust )->getVehicleSystems( );
    if( vehicleSystems == nullptr )
    {
        throw std::runtime_error( "Error when creating thrust acceleration model for " + nameOfBodyUndergoingThrust +
                                  ", body has no systems" );
    }

    std::map< std::string, std::shared_ptr< system_models::EngineModel > > engineModels = vehicleSystems->getEngineModels( );

    if( thrustAccelerationSettings->useAllEngines_ == true )
    {
        engineModelsForAcceleration = utilities::createVectorFromMapValues( engineModels );
    }
    else
    {
        for( unsigned int i = 0; i < thrustAccelerationSettings->engineIds_.size( ); i++ )
        {
            if( engineModels.count( thrustAccelerationSettings->engineIds_.at( i ) ) == 0 )
            {
                throw std::runtime_error( "Error when retrieving engine " + thrustAccelerationSettings->engineIds_.at( i ) +
                                          " for thrust acceleration. Engine with this name not found. " );
            }
            else
            {
                engineModelsForAcceleration.push_back(
                            engineModels.at( thrustAccelerationSettings->engineIds_.at( i ) ) );
            }
        }
    }

    // Add required updates of environemt models.
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > totalUpdateSettings;
    propagators::addEnvironmentUpdates( totalUpdateSettings, magnitudeUpdateSettings );
    propagators::addEnvironmentUpdates( totalUpdateSettings, directionUpdateSettings );

    return std::make_shared< propulsion::ThrustAcceleration >(
                engineModelsForAcceleration,
                thrustDirectionWrapper,
                std::bind( &Body::getBodyMass, bodies.at( nameOfBodyUndergoingThrust ) ),
                totalUpdateSettings );
}

//! Function to create a direct tical acceleration model, according to approach of Lainey et al. (2007, 2009, ...)
std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > createDirectTidalDissipationAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Check input consistency
    std::shared_ptr< DirectTidalDissipationAccelerationSettings > tidalAccelerationSettings =
            std::dynamic_pointer_cast< DirectTidalDissipationAccelerationSettings >( accelerationSettings );
    if( tidalAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, input is inconsistent" );
    }

    std::function< double( ) > gravitationalParameterFunctionOfBodyExertingTide;
    std::function< double( ) > gravitationalParameterFunctionOfBodyUndergoingTide;

    if( tidalAccelerationSettings->useTideRaisedOnPlanet_ )
    {
        if( bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, satellite " +
                                      nameOfBodyUndergoingAcceleration + " has no gravity field" );
        }
        else
        {
            gravitationalParameterFunctionOfBodyUndergoingTide = std::bind(
                        &GravityFieldModel::getGravitationalParameter, bodyUndergoingAcceleration->getGravityFieldModel( ) );
        }
    }
    else
    {
        if( bodyExertingAcceleration->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, satellite " +
                                      nameOfBodyExertingAcceleration + " has no gravity field" );
        }
        else
        {
            gravitationalParameterFunctionOfBodyExertingTide = std::bind(
                        &GravityFieldModel::getGravitationalParameter, bodyExertingAcceleration->getGravityFieldModel( ) );
        }


        if( bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, satellite " +
                                      nameOfBodyUndergoingAcceleration + " has no gravity field" );
        }
        else
        {
            gravitationalParameterFunctionOfBodyUndergoingTide = std::bind(
                        &GravityFieldModel::getGravitationalParameter, bodyUndergoingAcceleration->getGravityFieldModel( ) );
        }
    }

    double referenceRadius = TUDAT_NAN;
    if( tidalAccelerationSettings->useTideRaisedOnPlanet_ )
    {
        if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, planet " +
                                      nameOfBodyExertingAcceleration + " has no s.h. gravity field" );
        }
        else
        {
            referenceRadius = std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        bodyExertingAcceleration->getGravityFieldModel( ) )->getReferenceRadius( );
        }
    }
    else
    {
        if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                    bodyUndergoingAcceleration->getGravityFieldModel( ) ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, planet " +
                                      nameOfBodyUndergoingAcceleration + " has no s.h. gravity field" );
        }
        else
        {
            referenceRadius = std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        bodyUndergoingAcceleration->getGravityFieldModel( ) )->getReferenceRadius( );
        }
    }
    

    if( tidalAccelerationSettings->useTideRaisedOnPlanet_ )
    {
        std::function< Eigen::Vector3d( ) > planetAngularVelocityVectorFunction =
                std::bind( &Body::getCurrentAngularVelocityVectorInGlobalFrame, bodyExertingAcceleration );

        // Create direct tidal model from tidal time lag directly
        if ( std::isnan( tidalAccelerationSettings->inverseTidalQualityFactor_ ) && std::isnan( tidalAccelerationSettings->tidalPeriod_ ) )
        {
            return std::make_shared< DirectTidalDissipationAcceleration >(
                    std::bind( &Body::getState, bodyUndergoingAcceleration ), std::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParameterFunctionOfBodyUndergoingTide, planetAngularVelocityVectorFunction,
                    tidalAccelerationSettings->k2LoveNumber_, tidalAccelerationSettings->timeLag_,
                    referenceRadius, tidalAccelerationSettings->includeDirectRadialComponent_);
        }
        // Create direct tidal model from Q and T
        else
        {
            return std::make_shared< DirectTidalDissipationAcceleration >(
                    std::bind( &Body::getState, bodyUndergoingAcceleration ), std::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParameterFunctionOfBodyUndergoingTide, planetAngularVelocityVectorFunction,
                    tidalAccelerationSettings->k2LoveNumber_, tidalAccelerationSettings->inverseTidalQualityFactor_,
                    tidalAccelerationSettings->tidalPeriod_, referenceRadius, tidalAccelerationSettings->includeDirectRadialComponent_);
        }
    }

    else if( !tidalAccelerationSettings->explicitLibraionalTideOnSatellite_ )
    {
        // Create direct tidal model from tidal time lag directly
        if ( std::isnan( tidalAccelerationSettings->inverseTidalQualityFactor_ ) && std::isnan( tidalAccelerationSettings->tidalPeriod_ ) )
        {
            return std::make_shared< DirectTidalDissipationAcceleration >(
                    std::bind( &Body::getState, bodyUndergoingAcceleration ), std::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParameterFunctionOfBodyExertingTide, gravitationalParameterFunctionOfBodyUndergoingTide,
                    tidalAccelerationSettings->k2LoveNumber_, tidalAccelerationSettings->timeLag_,
                    referenceRadius, tidalAccelerationSettings->includeDirectRadialComponent_);
        }
        // Create direct tidal model from Q and T
        else
        {
            return std::make_shared< DirectTidalDissipationAcceleration >(
                    std::bind( &Body::getState, bodyUndergoingAcceleration ), std::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParameterFunctionOfBodyExertingTide, gravitationalParameterFunctionOfBodyUndergoingTide,
                    tidalAccelerationSettings->k2LoveNumber_, tidalAccelerationSettings->inverseTidalQualityFactor_, tidalAccelerationSettings->tidalPeriod_,
                    referenceRadius, tidalAccelerationSettings->includeDirectRadialComponent_);
        }
    }

    else
    {
        std::function< Eigen::Vector3d( ) > moonAngularVelocityVectorFunction =
                std::bind( &Body::getCurrentAngularVelocityVectorInGlobalFrame, bodyExertingAcceleration );

        // Create direct tidal model from tidal time lag directly
        if ( std::isnan( tidalAccelerationSettings->inverseTidalQualityFactor_ ) && std::isnan( tidalAccelerationSettings->tidalPeriod_ ) )
        {
            return std::make_shared< DirectTidalDissipationAcceleration >(
                    std::bind( &Body::getState, bodyUndergoingAcceleration ), std::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParameterFunctionOfBodyExertingTide, gravitationalParameterFunctionOfBodyUndergoingTide,
                    moonAngularVelocityVectorFunction, tidalAccelerationSettings->k2LoveNumber_,
                    tidalAccelerationSettings->timeLag_, referenceRadius, tidalAccelerationSettings->includeDirectRadialComponent_);
        }
        // Create direct tidal model from Q and T
        else
        {
            return std::make_shared< DirectTidalDissipationAcceleration >(
                    std::bind( &Body::getState, bodyUndergoingAcceleration ), std::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParameterFunctionOfBodyExertingTide, gravitationalParameterFunctionOfBodyUndergoingTide,
                    moonAngularVelocityVectorFunction, tidalAccelerationSettings->k2LoveNumber_, tidalAccelerationSettings->inverseTidalQualityFactor_,
                    tidalAccelerationSettings->tidalPeriod_, referenceRadius, tidalAccelerationSettings->includeDirectRadialComponent_);
        }
    }
}

//! Function to create a momentum wheel desaturation acceleration model.
std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > createMomentumWheelDesaturationAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Check input consistency
    std::shared_ptr< MomentumWheelDesaturationAccelerationSettings > desaturationAccelerationSettings =
            std::dynamic_pointer_cast< MomentumWheelDesaturationAccelerationSettings >( accelerationSettings );
    if( desaturationAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating momentum wheel desaturation acceleration, input is inconsistent" );
    }

    if( nameOfBodyUndergoingAcceleration != nameOfBodyExertingAcceleration )
    {
        throw std::runtime_error( "Error when creating momentum wheel desaturation acceleration, exerting and undergoing bodies are not the same" );
    }

    // Return desaturation acceleration model.
    return std::make_shared< propulsion::MomentumWheelDesaturationThrustAcceleration >(
                desaturationAccelerationSettings->thrustMidTimes_,
                desaturationAccelerationSettings->deltaVValues_,
                desaturationAccelerationSettings->totalManeuverTime_,
                desaturationAccelerationSettings->maneuverRiseTime_ );
}

//! Function to create acceleration model object.
std::shared_ptr< AccelerationModel< Eigen::Vector3d > > createAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfCentralBody,
        const SystemOfBodies& bodies )
{
    // Declare pointer to return object.
    std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;

    // Switch to call correct acceleration model type factory function.
    switch( accelerationSettings->accelerationType_ )
    {
    case point_mass_gravity:
    case spherical_harmonic_gravity:
    case mutual_spherical_harmonic_gravity:
    case polyhedron_gravity:
    case ring_gravity:
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
    case radiation_pressure:
        accelerationModelPointer = createRadiationPressureAccelerationModel(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    bodies,
                    accelerationSettings );
        break;
    case cannon_ball_radiation_pressure:
        accelerationModelPointer = createCannonballRadiationPressureAcceleratioModel(
            bodyUndergoingAcceleration,
            bodyExertingAcceleration,
            nameOfBodyUndergoingAcceleration,
            nameOfBodyExertingAcceleration,
            bodies);
        break;
    case thrust_acceleration:
        accelerationModelPointer = createThrustAcceleratioModel(
                    accelerationSettings, bodies,
                    nameOfBodyUndergoingAcceleration );
        break;
    case relativistic_correction_acceleration:
        accelerationModelPointer = createRelativisticCorrectionAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings, bodies );
        break;
    case empirical_acceleration:
        accelerationModelPointer = createEmpiricalAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case direct_tidal_dissipation_in_central_body_acceleration:
        accelerationModelPointer = createDirectTidalDissipationAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case direct_tidal_dissipation_in_orbiting_body_acceleration:
        accelerationModelPointer = createDirectTidalDissipationAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case momentum_wheel_desaturation_acceleration:
        accelerationModelPointer = createMomentumWheelDesaturationAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case yarkovsky_acceleration:
        accelerationModelPointer = createYarkovskyAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case custom_acceleration:
        accelerationModelPointer = createCustomAccelerationModel(
                    accelerationSettings,
                    nameOfBodyUndergoingAcceleration );
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


void addEihAccelerations(
    const SystemOfBodies& bodies,
    const std::map< std::string, std::vector< std::string > > orderedEihBodies,
    const std::map< std::string, std::string >& centralBodies,
    basic_astrodynamics::AccelerationMap& accelerationMap )
{
    std::vector< std::string > eihExertingBodies;
    std::vector< std::string > eihUndergoingBodies;
    std::vector< std::vector< std::string > > eihBodyListToCheck;

    for( auto it : centralBodies )
    {
        if( it.second != "SSB" )
        {
            throw std::runtime_error( "Error when creating EIH accelerations, implementation currently requires central body equal to SSB" );
        }
    }
    if( bodies.getFrameOrigin( ) != "SSB" )
    {
        throw std::runtime_error( "Error when creating EIH accelerations, implementation currently requires global frame origin to be equal to SSB" );
    }

    // Create list of bodies exerting acceleration, and bodies involved in the acceleration of each body
    for( auto it : orderedEihBodies )
    {
        eihUndergoingBodies.push_back( it.first );
        eihBodyListToCheck.push_back( it.second );
        eihBodyListToCheck.at( eihBodyListToCheck.size( ) - 1 ).push_back( it.first );
    }

    // Create sorted list of bodies involved in first body's accelerations, and check against duplicates
    std::vector< std::string > firstBodyInvolvedBodies = eihBodyListToCheck.at( 0 );
    std::sort( firstBodyInvolvedBodies.begin( ), firstBodyInvolvedBodies.end( ) );
    bool duplicateEntriesExist = std::adjacent_find(firstBodyInvolvedBodies.begin(), firstBodyInvolvedBodies.end()) != firstBodyInvolvedBodies.end();
    if( duplicateEntriesExist )
    {
        throw std::runtime_error( "Error when creating EIH equations, acceleration settings for " +
                                  eihUndergoingBodies.at( 0 ) + " duplicate entries found for involved bodies" );
    }

    // Check if all bodies require the same accelerating bodies
    for( unsigned int i = 1; i < eihBodyListToCheck.size( ); i++ )
    {
        if( eihBodyListToCheck.at( 0 ).size( ) != eihBodyListToCheck.at( i ).size( ) )
        {
            throw std::runtime_error( "Error when creating EIH equations, acceleration settings for " +
                                      eihUndergoingBodies.at( 0 ) + " and " +
                                      eihUndergoingBodies.at( i ) + " are incompatible (different number of bodies exerting EIH acceleration found)" );
        }
        std::vector< std::string > currentBodyInvolvedBodies = eihBodyListToCheck.at( i );
        std::sort( currentBodyInvolvedBodies.begin( ), currentBodyInvolvedBodies.end( ) );

        for( unsigned int j = 0; j < firstBodyInvolvedBodies.size( ); j++ )
        {
            if( firstBodyInvolvedBodies.at( j ) != currentBodyInvolvedBodies.at( j ) )
            {
                throw std::runtime_error( "Error when creating EIH equations, acceleration settings for " +
                                          eihExertingBodies.at( 0 ) + " and " +
                                          eihExertingBodies.at( i ) + " are incompatible (ordered list of involved bodies are different)" );
            }
        }
    }


    eihExertingBodies = eihUndergoingBodies;
    for( unsigned int j = 0; j < firstBodyInvolvedBodies.size( ); j++ )
    {
        if( std::find( eihExertingBodies.begin( ),
                       eihExertingBodies.end( ),
                       firstBodyInvolvedBodies.at( j ) ) == eihExertingBodies.end( ) )
        {
            eihExertingBodies.push_back( firstBodyInvolvedBodies.at( j ) );
        }
    }


    std::vector< std::function< double( ) > > gravitationalParameterFunction;
    std::vector< std::function< Eigen::Matrix< double, 6, 1 >( ) > > bodyStateFunctions;
    for( unsigned int i = 0; i < eihExertingBodies.size( ); i++ )
    {
        gravitationalParameterFunction.push_back( std::bind( &Body::getGravitationalParameter, bodies.at( eihExertingBodies.at( i ) ) ) );
        bodyStateFunctions.push_back( std::bind( &Body::getState, bodies.at( eihExertingBodies.at( i ) ) ) );

    }

    std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations =
        std::make_shared< relativity::EinsteinInfeldHoffmannEquations >(
            eihUndergoingBodies, eihExertingBodies,
            gravitationalParameterFunction,
            bodyStateFunctions,
            std::bind( &relativity::PPNParameterSet::getParameterGamma, relativity::ppnParameterSet ),
            std::bind( &relativity::PPNParameterSet::getParameterBeta, relativity::ppnParameterSet ) );

    for( unsigned int i = 0; i < eihUndergoingBodies.size( ); i++ )
    {
        std::vector< std::string > currentEihExertingAccelerations = eihExertingBodies;
        std::vector< std::string >::iterator entryToRemove = std::find(
            currentEihExertingAccelerations.begin(), currentEihExertingAccelerations.end(), eihUndergoingBodies.at( i ) );
        if (entryToRemove != currentEihExertingAccelerations.end())
        {
            currentEihExertingAccelerations.erase( entryToRemove );
        }

        accelerationMap[ eihUndergoingBodies.at( i ) ][ "" ].push_back( std::make_shared< relativity::EinsteinInfeldHoffmannAcceleration >(
            eihEquations, eihUndergoingBodies.at( i ), eihExertingBodies ) );
    }
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
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > >
                accelerationsForBody = bodyIterator->second;

        // Retrieve indices of all acceleration anf thrust models.
        std::vector< int > aerodynamicAccelerationIndices;
        std::vector< int > thrustAccelerationIndices;

        std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > >
                currentBodyAccelerations;
        int counter = 0;
        // Iterate over all bodies exerting an acceleration
        for( std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > >::
             iterator body2Iterator = accelerationsForBody.begin( );
             body2Iterator != accelerationsForBody.end( ); body2Iterator++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = body2Iterator->first;
            std::vector< std::shared_ptr< AccelerationSettings > > accelerationList = body2Iterator->second;
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
                std::pair< std::string, std::shared_ptr< AccelerationSettings > >  currentAccelerationPair =
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
                std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > >
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

inline basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
    const SystemOfBodies& bodies,
    const SelectedAccelerationMap& selectedAccelerationPerBody,
    const std::map< std::string, std::string >& centralBodies )
{
    // Declare return map.
    basic_astrodynamics::AccelerationMap accelerationModelMap;
    std::map< std::string, std::vector< std::string > > orderedEihBodies;

    // Put selectedAccelerationPerBody in correct order
    SelectedAccelerationList orderedAccelerationPerBody =
        orderSelectedAccelerationMap( selectedAccelerationPerBody );

    // Iterate over all bodies which are undergoing acceleration
    for( SelectedAccelerationList::const_iterator bodyIterator =
        orderedAccelerationPerBody.begin( ); bodyIterator != orderedAccelerationPerBody.end( );
         bodyIterator++ )
    {
        std::shared_ptr< Body > currentCentralBody;

        // Retrieve name of body undergoing acceleration.
        std::string bodyUndergoingAcceleration = bodyIterator->first;

        // Retrieve name of current central body.
        std::string currentCentralBodyName = centralBodies.at( bodyUndergoingAcceleration );

        if( !ephemerides::isFrameInertial( currentCentralBodyName ) )
        {
            if( bodies.count( currentCentralBodyName ) == 0 )
            {
                throw std::runtime_error(
                    std::string( "Error, could not find non-inertial central body ") +
                    currentCentralBodyName + " of " + bodyUndergoingAcceleration +
                    " when making acceleration model." );
            }
            else
            {
                currentCentralBody = bodies.at( currentCentralBodyName );
            }
        }

        // Check if body undergoing acceleration is included in bodies
        if( bodies.count( bodyUndergoingAcceleration ) ==  0 )
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
        std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > >
            accelerationsForBody = bodyIterator->second;

        std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > > thrustAccelerationSettings;

        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > currentAcceleration;

        // Iterate over all bodies exerting an acceleration
        for( unsigned int i = 0; i < accelerationsForBody.size( ); i++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = accelerationsForBody.at( i ).first;

            // Check if body exerting acceleration is included in bodies
            if( bodies.count( bodyExertingAcceleration ) ==  0 )
            {
                throw std::runtime_error(
                    std::string( "Error when making acceleration models, requested forces ")
                    + "acting on body " + bodyUndergoingAcceleration  + " due to body " +
                    bodyExertingAcceleration +
                    ", but no such body found in map of bodies" );
            }

            if( accelerationsForBody.at( i ).second->accelerationType_ == basic_astrodynamics::thrust_acceleration )
            {
                thrustAccelerationSettings.push_back( accelerationsForBody.at( i ) );
            }
            else if( accelerationsForBody.at( i ).second->accelerationType_ == basic_astrodynamics::einstein_infeld_hoffmann_acceleration )
            {
                if( orderedEihBodies.count( bodyUndergoingAcceleration ) > 0 )
                {
                    if( std::find( orderedEihBodies.at( bodyUndergoingAcceleration ).begin( ),
                                   orderedEihBodies.at( bodyUndergoingAcceleration ).end( ),
                                   bodyExertingAcceleration ) !=
                        orderedEihBodies.at( bodyUndergoingAcceleration ).end( ) )
                    {
                        throw std::runtime_error( "Error when parsing EIH acceleration settings, found combinatin of bodies " +
                                                  bodyUndergoingAcceleration + ", " + bodyExertingAcceleration  + " multiple times." );
                    }
                }

                orderedEihBodies[ bodyUndergoingAcceleration ].push_back( bodyExertingAcceleration );
            }
            else
            {
                currentAcceleration = createAccelerationModel( bodies.at( bodyUndergoingAcceleration ),
                                                               bodies.at( bodyExertingAcceleration ),
                                                               accelerationsForBody.at( i ).second,
                                                               bodyUndergoingAcceleration,
                                                               bodyExertingAcceleration,
                                                               currentCentralBody,
                                                               currentCentralBodyName,
                                                               bodies );


                // Create acceleration model.
                mapOfAccelerationsForBody[ bodyExertingAcceleration ].push_back(
                    currentAcceleration );
            }
        }

        // Create thrust accelerations last
        for( unsigned int i = 0; i < thrustAccelerationSettings.size( ); i++ )
        {
            currentAcceleration = createAccelerationModel( bodies.at( bodyUndergoingAcceleration ),
                                                           bodies.at( thrustAccelerationSettings.at( i ).first ),
                                                           thrustAccelerationSettings.at( i ).second,
                                                           bodyUndergoingAcceleration,
                                                           thrustAccelerationSettings.at( i ).first,
                                                           currentCentralBody,
                                                           currentCentralBodyName,
                                                           bodies );


            // Create acceleration model.
            mapOfAccelerationsForBody[ thrustAccelerationSettings.at( i ).first  ].push_back(
                currentAcceleration );
        }
        // Put acceleration models on current body in return map.
        accelerationModelMap[ bodyUndergoingAcceleration ] = mapOfAccelerationsForBody;
    }

    if( orderedEihBodies.size( ) > 0 )
    {
        addEihAccelerations(
            bodies, orderedEihBodies, centralBodies, accelerationModelMap );
    }

    return accelerationModelMap;
}

//! Function to create acceleration models from a map of bodies and acceleration model types.
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const SystemOfBodies& bodies,
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

    return createAccelerationModelsMap( bodies, selectedAccelerationPerBody, centralBodyMap );
}

} // namespace simulation_setup

} // namespace tudat
