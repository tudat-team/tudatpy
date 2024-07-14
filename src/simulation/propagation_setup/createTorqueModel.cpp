/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create an aerodynamic torque model.
std::shared_ptr< basic_astrodynamics::InertialTorqueModel > createInertialTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::string& nameOfBodyUndergoingTorque )
{
    std::function< Eigen::Vector3d( ) > angularVelocityFunction =
            std::bind( &Body::getCurrentAngularVelocityVectorInLocalFrame, bodyUndergoingTorque );
    std::function< Eigen::Matrix3d( ) > inertiaTensorFunction =
            std::bind( &Body::getBodyInertiaTensor, bodyUndergoingTorque );
    return std::make_shared< basic_astrodynamics::InertialTorqueModel >(
                angularVelocityFunction, inertiaTensorFunction );
}

//! Function to create an aerodynamic torque model.
std::shared_ptr< aerodynamics::AerodynamicTorque > createAerodynamicTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque )
{
    // Check existence of required environment models
    if( bodyUndergoingTorque->getAerodynamicCoefficientInterface( ) == nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic torque, body " +
                                  nameOfBodyUndergoingTorque +
                                  "has no aerodynamic coefficients." );
    }

    if( bodyExertingTorque->getAtmosphereModel( ) == nullptr )
    {
        throw std::runtime_error(  "Error when making aerodynamic torque, central body " +
                                   nameOfBodyExertingTorque + " has no atmosphere model.");
    }

    if( bodyExertingTorque->getShapeModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic torque, central body " +
                                  nameOfBodyExertingTorque + " has no shape model." );
    }

    // Retrieve flight conditions; create object if not yet extant.
    std::shared_ptr< aerodynamics::AtmosphericFlightConditions > bodyFlightConditions =
            std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                bodyUndergoingTorque->getFlightConditions( ) );

    if( bodyFlightConditions == nullptr && bodyUndergoingTorque->getFlightConditions( ) == nullptr )
    {
        bodyFlightConditions = createAtmosphericFlightConditions( bodyUndergoingTorque,
                                                       bodyExertingTorque,
                                                       nameOfBodyUndergoingTorque,
                                                       nameOfBodyExertingTorque ) ;
        bodyUndergoingTorque->setFlightConditions(
                    bodyFlightConditions );
    }
    else if( bodyFlightConditions == nullptr && bodyUndergoingTorque->getFlightConditions( ) != nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic torque, found flight conditions that are not atmospheric." );
    }

    // Retrieve frame in which aerodynamic coefficients are defined.
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > aerodynamicCoefficients =
            bodyUndergoingTorque->getAerodynamicCoefficientInterface( );
    reference_frames::AerodynamicsReferenceFrames torqueFrame = aerodynamics::getCompleteFrameForCoefficients(
            aerodynamicCoefficients->getMomentCoefficientsFrame( ) );

    // Create function to transform from frame of aerodynamic coefficienrs to that of propagation.
    std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > toPropagationFrameTransformation;
    toPropagationFrameTransformation =
            reference_frames::getAerodynamicForceTransformationFunction(
                bodyFlightConditions->getAerodynamicAngleCalculator( ),
                torqueFrame,
                std::bind( &Body::getCurrentRotationToGlobalFrame, bodyExertingTorque ),
                reference_frames::body_frame );


    std::function< Eigen::Vector3d( ) > coefficientFunction =
            std::bind( &aerodynamics::AerodynamicCoefficientInterface::getCurrentMomentCoefficients,
                         aerodynamicCoefficients );
    std::function< Eigen::Vector3d( ) > coefficientInPropagationFrameFunction =
            std::bind( &reference_frames::transformVectorFunctionFromVectorFunctions,
                         coefficientFunction, toPropagationFrameTransformation );

    // Create torque model.
    return std::make_shared< aerodynamics::AerodynamicTorque >(
                coefficientInPropagationFrameFunction,
                std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentDensity, bodyFlightConditions ),
                std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeed, bodyFlightConditions ),
                std::bind( &aerodynamics::AerodynamicCoefficientInterface::getReferenceArea, aerodynamicCoefficients ),
                std::bind( &aerodynamics::AerodynamicCoefficientInterface::getReferenceLengths, aerodynamicCoefficients ),
                aerodynamics::areCoefficientsInNegativeDirection(
                        aerodynamicCoefficients->getMomentCoefficientsFrame( ) ) );
}

//! Function to create a second-degree gravitational torque.
std::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > createSecondDegreeGravitationalTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque )
{
    // Retrieve state functions
    std::function< Eigen::Vector3d( ) > positionOfBodySubjectToTorqueFunction =
            std::bind( &simulation_setup::Body::getPosition, bodyUndergoingTorque );
    std::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction =
            std::bind( &simulation_setup::Body::getPosition, bodyExertingTorque );

    // Check model availability
    std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel = bodyExertingTorque->getGravityFieldModel( );
    std::function< double( ) > gravitationalParameterOfAttractingBodyFunction;
    if( gravityFieldModel ==  nullptr )
    {
        throw std::runtime_error( "Error when making second degree gravitational torque, " + nameOfBodyExertingTorque +
                                  " does not possess a gravity field" );
    }
    else
    {
        gravitationalParameterOfAttractingBodyFunction = std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                                      gravityFieldModel );
    }

    // Retrieve environment parameters
    std::function< Eigen::Matrix3d( ) > inertiaTensorOfRotatingBodyFunction  =
            std::bind( &simulation_setup::Body::getBodyInertiaTensor, bodyUndergoingTorque );
    std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameFunction =
            std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodyUndergoingTorque );

    return std::make_shared<gravitation::SecondDegreeGravitationalTorqueModel >(
                positionOfBodySubjectToTorqueFunction, gravitationalParameterOfAttractingBodyFunction,
                inertiaTensorOfRotatingBodyFunction, positionOfBodyExertingTorqueFunction, rotationToBodyFixedFrameFunction );

}


//! Function to create a spherical harmonic gravitational torque
std::shared_ptr< gravitation::SphericalHarmonicGravitationalTorqueModel > createSphericalHarmonicGravitationalTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::shared_ptr< TorqueSettings > torqueSettings,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque )
{
    std::shared_ptr< SphericalHarmonicTorqueSettings > sphericalHarmonicTorqueSettings =
            std::dynamic_pointer_cast< SphericalHarmonicTorqueSettings >( torqueSettings );

    if( sphericalHarmonicTorqueSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating spherical harmonic torque, input is inconsistent" );
    }
    std::shared_ptr< AccelerationSettings > sphericalHarmonicAccelerationSettings =
            std::make_shared< SphericalHarmonicAccelerationSettings >(
                sphericalHarmonicTorqueSettings->maximumDegree_,
                sphericalHarmonicTorqueSettings->maximumOrder_ );
    std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                 createSphericalHarmonicsGravityAcceleration(
                    bodyExertingTorque, bodyUndergoingTorque, nameOfBodyExertingTorque, nameOfBodyUndergoingTorque,
                    sphericalHarmonicAccelerationSettings, false, false ) );

    return std::make_shared< gravitation::SphericalHarmonicGravitationalTorqueModel >(
                sphericalHarmonicAcceleration,
                std::bind( &Body::getCurrentRotationToLocalFrame, bodyUndergoingTorque ),
                std::bind( &Body::getBodyMass, bodyExertingTorque ) );
}


//! Function to create a spherical harmonic gravitational torque
std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureTorque > createRadiationPressureTorqueModel(
    const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
    const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
    const std::shared_ptr< TorqueSettings > torqueSettings,
    const std::string& nameOfBodyUndergoingTorque,
    const std::string& nameOfBodyExertingTorque,
    const SystemOfBodies& bodies )
{

    std::shared_ptr< AccelerationSettings > sphericalHarmonicAccelerationSettings =
        std::make_shared< AccelerationSettings >( basic_astrodynamics::radiation_pressure );

    std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration > radiationPressureAcceleration =
        std::dynamic_pointer_cast< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration >(
            createRadiationPressureAccelerationModel(
                bodyExertingTorque, bodyUndergoingTorque, nameOfBodyExertingTorque, nameOfBodyUndergoingTorque,
                bodies ) );

    if( radiationPressureAcceleration == nullptr )
    {
        throw std::runtime_error( "Error when creating radiation pressure torque, only isotropic source is supported" );
    }

    std::function< Eigen::Vector3d( ) > centerOfMassFunction = nullptr;
    std::shared_ptr< RigidBodyProperties > massProperties = bodyUndergoingTorque->getMassProperties( );
    if( massProperties == nullptr )
    {
        throw std::runtime_error( "Error when making radiation pressure torque model, no mass properties are defined, which are required for center of mass definition" );
    }
    else
    {
        centerOfMassFunction = std::bind( &RigidBodyProperties::getCurrentCenterOfMass, massProperties );
    }

    return std::make_shared< electromagnetism::IsotropicPointSourceRadiationPressureTorque >(
        radiationPressureAcceleration, centerOfMassFunction );
}

std::shared_ptr< basic_astrodynamics::CustomTorqueModel > createCustomTorqueModel(
        const std::shared_ptr< TorqueSettings > torqueSettings,
        const std::string& nameOfBodyUndergoingTorque )
{
    std::shared_ptr< CustomTorqueSettings > customTorqueSettings =
            std::dynamic_pointer_cast< CustomTorqueSettings >(
                torqueSettings );
    if( customTorqueSettings == nullptr )
    {
        throw std::runtime_error( "Error, expected custom torque settings when making torque model on " +
                                  nameOfBodyUndergoingTorque  );
    }
    return std::make_shared< basic_astrodynamics::CustomTorqueModel >( customTorqueSettings->torqueFunction_ );

}

//! Function to create torque model object.
std::shared_ptr< basic_astrodynamics::TorqueModel > createTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::shared_ptr< TorqueSettings > torqueSettings,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque,
        const SystemOfBodies& bodies )
{
    std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel;

    switch( torqueSettings->torqueType_ )
    {
    case basic_astrodynamics::inertial_torque:
    {
        torqueModel = createInertialTorqueModel(
                    bodyUndergoingTorque, nameOfBodyUndergoingTorque );
        break;
    }
    case basic_astrodynamics::second_order_gravitational_torque:
    {
        torqueModel = createSecondDegreeGravitationalTorqueModel(
                    bodyUndergoingTorque, bodyExertingTorque, nameOfBodyUndergoingTorque, nameOfBodyExertingTorque );
        break;
    }
    case basic_astrodynamics::aerodynamic_torque:
    {
        torqueModel = createAerodynamicTorqueModel(
                    bodyUndergoingTorque, bodyExertingTorque, nameOfBodyUndergoingTorque, nameOfBodyExertingTorque );
        break;
    }
    case basic_astrodynamics::spherical_harmonic_gravitational_torque:
    {
        torqueModel = createSphericalHarmonicGravitationalTorqueModel(
                    bodyUndergoingTorque, bodyExertingTorque, torqueSettings, nameOfBodyUndergoingTorque, nameOfBodyExertingTorque );
        break;
    }
    case basic_astrodynamics::radiation_pressure_torque:
    {
        torqueModel = createRadiationPressureTorqueModel(
            bodyUndergoingTorque, bodyExertingTorque, torqueSettings, nameOfBodyUndergoingTorque, nameOfBodyExertingTorque, bodies );
        break;
    }
    case basic_astrodynamics::custom_torque:
    {
        torqueModel = createCustomTorqueModel(
                    torqueSettings, nameOfBodyUndergoingTorque );
        break;
    }
    default:
        throw std::runtime_error(
                    "Error, did not recognize type " + std::to_string( torqueSettings->torqueType_ ) +
                    " when making torque model" );
    }

    return torqueModel;
}


//! Function to create torque models from a map of bodies and torque model settings.
basic_astrodynamics::TorqueModelMap createTorqueModelsMap(
        const SystemOfBodies& bodies,
        SelectedTorqueMap selectedTorquePerBody,
        const std::vector< std::string >& propagatedBodies )
{
    for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
    {
        if( selectedTorquePerBody.count( propagatedBodies.at( i ) ) == 0 )
        {
            selectedTorquePerBody[ propagatedBodies.at( i ) ] =
                    std::map< std::string, std::vector< std::shared_ptr< TorqueSettings > > >( );
        }
    }

    basic_astrodynamics::TorqueModelMap torqueModelMap;

    for( SelectedTorqueMap::const_iterator acceleratedBodyIterator = selectedTorquePerBody.begin( );
         acceleratedBodyIterator != selectedTorquePerBody.end( ); acceleratedBodyIterator++ )
    {
        if( bodies.count( acceleratedBodyIterator->first ) == 0 )
        {
            throw std::runtime_error(
                        "Error, could not find body " + acceleratedBodyIterator->first + " when making torque model map." );
        }
        else
        {
            torqueModelMap[ acceleratedBodyIterator->first ][ acceleratedBodyIterator->first ].push_back(
                        createTorqueModel(
                        bodies.at( acceleratedBodyIterator->first ), bodies.at( acceleratedBodyIterator->first ),
                        std::make_shared< TorqueSettings >( basic_astrodynamics::inertial_torque ),
                        acceleratedBodyIterator->first, acceleratedBodyIterator->first, bodies ) );

            for( std::map< std::string, std::vector< std::shared_ptr< TorqueSettings > > >::const_iterator
                 acceleratingBodyIterator = acceleratedBodyIterator->second.begin( );
                 acceleratingBodyIterator != acceleratedBodyIterator->second.end( ); acceleratingBodyIterator++ )
            {
                if( bodies.count( acceleratingBodyIterator->first ) == 0 )
                {
                    throw std::runtime_error(
                                "Error, could not find body " + acceleratingBodyIterator->first + " when making torque model map." );
                }

                for( unsigned int i = 0; i < acceleratingBodyIterator->second.size( ); i++ )
                {
                    torqueModelMap[ acceleratedBodyIterator->first ][ acceleratingBodyIterator->first ].push_back(
                                createTorqueModel(
                                bodies.at( acceleratedBodyIterator->first ), bodies.at( acceleratingBodyIterator->first ),
                                acceleratingBodyIterator->second.at( i ),
                                acceleratedBodyIterator->first, acceleratingBodyIterator->first, bodies ) );
                }
            }
        }
    }

    return torqueModelMap;
}


}  // namespace simulation_setup

}  // namespace tudat

