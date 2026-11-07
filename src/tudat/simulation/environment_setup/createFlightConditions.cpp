/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>
#include <string>
#include <iostream>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/simulation/environment_setup/environmentModelClosure.h"

namespace tudat
{

namespace simulation_setup
{

std::shared_ptr< reference_frames::AerodynamicAngleCalculator > createAerodynamicAngleCalculatorFromBodyStates(
        const std::shared_ptr< Body > bodyWithFlightConditions,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyExertingAcceleration )
{
    // Create function to rotate state from intertial to body-fixed frame.
    std::function< Eigen::Quaterniond( ) > rotationToFrameFunction = std::bind( &Body::getCurrentRotationToLocalFrame, centralBody );
    std::function< Eigen::Matrix3d( ) > rotationMatrixToFrameDerivativeFunction =
            std::bind( &Body::getCurrentRotationMatrixDerivativeToLocalFrame, centralBody );

    std::function< Eigen::Matrix< double, 6, 1 >( ) > bodyStateFunction = std::bind( &Body::getState, bodyWithFlightConditions );
    std::function< Eigen::Matrix< double, 6, 1 >( ) > centralBodyStateFunction = std::bind( &Body::getState, centralBody );

    std::function< Eigen::Matrix< double, 6, 1 >( ) > relativeBodyFixedStateFunction =
            std::bind( &ephemerides::transformRelativeStateToFrame< double >,
                       bodyStateFunction,
                       centralBodyStateFunction,
                       rotationToFrameFunction,
                       rotationMatrixToFrameDerivativeFunction );

    // Create aerodynamic angles calculator and set in flight conditions.
    return std::make_shared< reference_frames::AerodynamicAngleCalculator >(
            relativeBodyFixedStateFunction,
            std::bind( &simulation_setup::Body::getCurrentRotationToGlobalFrame, centralBody ),
            nameOfBodyExertingAcceleration,
            1 );
}

std::shared_ptr< reference_frames::AerodynamicAngleCalculator > createAerodynamicAngleCalculator(
        const std::shared_ptr< Body > bodyWithFlightConditions,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    if( std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >(
                bodyWithFlightConditions->getRotationalEphemeris( ) ) != nullptr )
    {
        return std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >(
                       bodyWithFlightConditions->getRotationalEphemeris( ) )
                ->getAerodynamicAngleCalculator( );
    }
    else
    {
        return createAerodynamicAngleCalculatorFromBodyStates( bodyWithFlightConditions, centralBody, nameOfBodyExertingAcceleration );
    }
}

bool applyFlightConditionsClosure( const std::shared_ptr< aerodynamics::FlightConditions > flightConditions,
                                   const std::shared_ptr< Body > bodyWithFlightConditions,
                                   const std::shared_ptr< Body > centralBody,
                                   const std::string& nameOfBodyUndergoingAcceleration,
                                   const std::string& nameOfBodyExertingAcceleration )
{
    // Check whether all required environment models are set.
    if( centralBody->getShapeModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making flight conditions, central body " + nameOfBodyExertingAcceleration +
                                  " has no shape model." );
    }
    if( centralBody->getRotationalEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when making flight conditions, central body " + nameOfBodyExertingAcceleration +
                                  " has no rotation model." );
    }

    // Retrieve aerodynamic angle calculator from aerodynamic angle based rotation model of body (if any), which must be closed
    std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > aerodynamicAngleRotationModel =
            std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >(
                    bodyWithFlightConditions->getRotationalEphemeris( ) );
    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator;
    if( aerodynamicAngleRotationModel != nullptr )
    {
        aerodynamicAngleCalculator = aerodynamicAngleRotationModel->getAerodynamicAngleCalculator( );
        if( aerodynamicAngleCalculator == nullptr )
        {
            return false;
        }
    }
    else
    {
        aerodynamicAngleCalculator =
                createAerodynamicAngleCalculatorFromBodyStates( bodyWithFlightConditions, centralBody, nameOfBodyExertingAcceleration );
    }

    std::shared_ptr< aerodynamics::AtmosphericFlightConditions > atmosphericFlightConditions =
            std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >( flightConditions );
    if( atmosphericFlightConditions != nullptr )
    {
        if( centralBody->getAtmosphereModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when making flight conditions, body " + nameOfBodyExertingAcceleration +
                                      " has no atmosphere model." );
        }

        // Add wind model (always exists)
        aerodynamicAngleCalculator->setWindModel( centralBody->getAtmosphereModel( )->getWindModel( ), centralBody->getShapeModel( ) );

        // Set rotation matrix derivative function for proper handling of non-rotating atmosphere
        aerodynamicAngleCalculator->setRotationMatrixDerivativeFunction(
                std::bind( &Body::getCurrentRotationMatrixDerivativeToLocalFrame, centralBody ) );

        // Set atmospheric rotation flag from wind model
        aerodynamicAngleCalculator->setIncludeAtmosphericRotation(
                centralBody->getAtmosphereModel( )->getWindModel( )->getIncludeCorotation( ) );

        std::function< double( const std::string& ) > controlSurfaceDeflectionFunction;
        if( bodyWithFlightConditions->getVehicleSystems( ) != nullptr )
        {
            controlSurfaceDeflectionFunction = std::bind( &system_models::VehicleSystems::getCurrentControlSurfaceDeflection,
                                                          bodyWithFlightConditions->getVehicleSystems( ),
                                                          std::placeholders::_1 );
        }
        atmosphericFlightConditions->applyEnvironmentModelClosure( centralBody->getAtmosphereModel( ),
                                                                   centralBody->getShapeModel( ),
                                                                   bodyWithFlightConditions->getAerodynamicCoefficientInterface( ),
                                                                   aerodynamicAngleCalculator,
                                                                   controlSurfaceDeflectionFunction );
    }
    else
    {
        flightConditions->applyEnvironmentModelClosure( centralBody->getShapeModel( ), aerodynamicAngleCalculator );
    }

    // Set aerodynamic angles w.r.t. body-fixed frame from rotation model of body (set by aerodynamic angle based rotation model
    // itself, if relevant)
    if( bodyWithFlightConditions->getRotationalEphemeris( ) != nullptr && aerodynamicAngleRotationModel == nullptr )
    {
        aerodynamicAngleCalculator->setBodyFixedAngleInterface(
                std::make_shared< reference_frames::FromGenericEphemerisAerodynamicAngleInterface >(
                        bodyWithFlightConditions->getRotationalEphemeris( ) ) );
    }
    return true;
}

bool applyEnvironmentModelClosure( const std::shared_ptr< aerodynamics::FlightConditions > flightConditions,
                                   const std::string& bodyName,
                                   const SystemOfBodies& bodies )
{
    if( bodies.count( flightConditions->getCentralBody( ) ) == 0 )
    {
        throw std::runtime_error( "Error when applying closure of flight conditions of body " + bodyName + ", central body " +
                                  flightConditions->getCentralBody( ) + " not found" );
    }
    return applyFlightConditionsClosure( flightConditions,
                                         bodies.at( bodyName ),
                                         bodies.at( flightConditions->getCentralBody( ) ),
                                         bodyName,
                                         flightConditions->getCentralBody( ) );
}

//! Function to create an atmospheric flight conditions object
std::shared_ptr< aerodynamics::AtmosphericFlightConditions > createAtmosphericFlightConditions(
        const std::shared_ptr< Body > bodyWithFlightConditions,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions =
            std::make_shared< aerodynamics::AtmosphericFlightConditions >( nameOfBodyExertingAcceleration );
    if( !applyFlightConditionsClosure( flightConditions,
                                       bodyWithFlightConditions,
                                       centralBody,
                                       nameOfBodyUndergoingAcceleration,
                                       nameOfBodyExertingAcceleration ) )
    {
        throw std::runtime_error( "Error when making flight conditions for " + nameOfBodyUndergoingAcceleration +
                                  ", aerodynamic angle based rotation model of the body is not closed" );
    }
    return flightConditions;
}

//! Function to create a flight conditions object
std::shared_ptr< aerodynamics::FlightConditions > createFlightConditions( const std::shared_ptr< Body > bodyWithFlightConditions,
                                                                          const std::shared_ptr< Body > centralBody,
                                                                          const std::string& nameOfBodyUndergoingAcceleration,
                                                                          const std::string& nameOfBodyExertingAcceleration )
{
    std::shared_ptr< aerodynamics::FlightConditions > flightConditions =
            std::make_shared< aerodynamics::FlightConditions >( nameOfBodyExertingAcceleration );
    if( !applyFlightConditionsClosure( flightConditions,
                                       bodyWithFlightConditions,
                                       centralBody,
                                       nameOfBodyUndergoingAcceleration,
                                       nameOfBodyExertingAcceleration ) )
    {
        throw std::runtime_error( "Error when making flight conditions for " + nameOfBodyUndergoingAcceleration +
                                  ", aerodynamic angle based rotation model of the body is not closed" );
    }
    return flightConditions;
}

std::shared_ptr< aerodynamics::FlightConditions > createFlightConditionsWithoutClosure( const SystemOfBodies& bodies,
                                                                                        const std::string& bodyName,
                                                                                        const std::string& centralBodyName )
{
    if( bodies.count( centralBodyName ) == 0 )
    {
        throw std::runtime_error( "Error when creating flight conditions for body " + bodyName + ", central body " + centralBodyName +
                                  " not found" );
    }
    if( bodies.at( centralBodyName )->getAtmosphereModel( ) != nullptr )
    {
        return std::make_shared< aerodynamics::AtmosphericFlightConditions >( centralBodyName );
    }
    else
    {
        return std::make_shared< aerodynamics::FlightConditions >( centralBodyName );
    }
}

void addFlightConditions( const SystemOfBodies& bodies, const std::string bodyName, const std::string centralBodyName )
{
    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when creating flight conditions for body " + bodyName + ", body no found" );
    }

    if( bodies.count( centralBodyName ) == 0 )
    {
        throw std::runtime_error( "Error when creating flight conditions for body " + bodyName + ", central body " + centralBodyName +
                                  " not found" );
    }

    std::shared_ptr< aerodynamics::FlightConditions > currentBodyFlightConditions = bodies.at( bodyName )->getFlightConditions( );
    if( currentBodyFlightConditions != nullptr )
    {
        if( currentBodyFlightConditions->getCentralBody( ) != centralBodyName )
        {
            throw std::runtime_error( "Error, flight conditions for " + bodyName + "w.r.t." + centralBodyName +
                                      " cannot be created, flight conditions w.r.t." + currentBodyFlightConditions->getCentralBody( ) +
                                      "already created." );
        }
        else if( bodies.at( centralBodyName )->getAtmosphereModel( ) != nullptr &&
                 std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >( currentBodyFlightConditions ) == nullptr )
        {
            throw std::runtime_error( "Error when creating flight conditions for body " + bodyName + "w.r.t." + centralBodyName +
                                      ", flight conditions already exist, but they are for an atmosphere-free body" );
        }
    }
    else
    {
        // Create flight conditions without closure, and apply the closure (as well as that of models waiting for these flight
        // conditions)
        bodies.at( bodyName )->setFlightConditions( createFlightConditionsWithoutClosure( bodies, bodyName, centralBodyName ) );
        bodies.at( bodyName )->setEnvironmentModelRequiresClosure( flight_conditions_closure );
    }

    closeSystemOfBodies( bodies );
    if( bodies.at( bodyName )->getEnvironmentModelClosureStates( ).count( flight_conditions_closure ) > 0 &&
        !bodies.at( bodyName )->getEnvironmentModelClosureStates( ).at( flight_conditions_closure ) )
    {
        throw std::runtime_error( "Error when creating flight conditions for body " + bodyName + " w.r.t. " + centralBodyName +
                                  ", closure of flight conditions could not be applied" );
    }
}

void addAtmosphericFlightConditions( const SystemOfBodies& bodies, const std::string bodyName, const std::string centralBodyName )
{
    addFlightConditions( bodies, bodyName, centralBodyName );
    if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >( bodies.at( bodyName )->getFlightConditions( ) ) == nullptr )
    {
        if( bodies.at( centralBodyName )->getAtmosphereModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when adding atmospheric flight conditions for " + bodyName + " w.r.t. " + centralBodyName +
                                      ", conditions could not be created, central body has no atmosphere" );
        }
        else
        {
            throw std::runtime_error( "Error when adding atmospheric flight conditions for " + bodyName + " w.r.t. " + centralBodyName +
                                      ", conditions could not be created" );
        }
    }
}

}  // namespace simulation_setup

}  // namespace tudat
