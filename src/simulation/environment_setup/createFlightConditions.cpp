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

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind/bind.hpp>

using namespace boost::placeholders;


#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/basics/deprecationWarnings.h"

namespace tudat
{

namespace simulation_setup
{

std::shared_ptr< reference_frames::AerodynamicAngleCalculator >  createAerodynamicAngleCalculator(
        const std::shared_ptr< Body > bodyWithFlightConditions,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    if( std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >(
                bodyWithFlightConditions->getRotationalEphemeris( ) ) != nullptr )
    {
        return std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >(
                    bodyWithFlightConditions->getRotationalEphemeris( ) )->getAerodynamicAngleCalculator( );
    }
    else
    {
        // Create function to rotate state from intertial to body-fixed frame.
        std::function< Eigen::Quaterniond( ) > rotationToFrameFunction =
                std::bind( &Body::getCurrentRotationToLocalFrame, centralBody );
        std::function< Eigen::Matrix3d( ) > rotationMatrixToFrameDerivativeFunction =
                std::bind( &Body::getCurrentRotationMatrixDerivativeToLocalFrame, centralBody );

        std::function< Eigen::Matrix< double, 6, 1 >( ) > bodyStateFunction = std::bind( &Body::getState, bodyWithFlightConditions );
        std::function< Eigen::Matrix< double, 6, 1 >( ) > centralBodyStateFunction = std::bind( &Body::getState, centralBody );

        std::function< Eigen::Matrix< double, 6, 1 >( ) > relativeBodyFixedStateFunction =
                std::bind( &ephemerides::transformRelativeStateToFrame< double >,
                           bodyStateFunction, centralBodyStateFunction,
                           rotationToFrameFunction,
                           rotationMatrixToFrameDerivativeFunction );

        // Create aerodynamic angles calculator and set in flight conditions.
        return std::make_shared< reference_frames::AerodynamicAngleCalculator >(
                    relativeBodyFixedStateFunction,
                    std::bind( &simulation_setup::Body::getCurrentRotationToGlobalFrame, centralBody ),
                    nameOfBodyExertingAcceleration, 1 );
    }

}

//! Function to create an atmospheric flight conditions object
std::shared_ptr< aerodynamics::AtmosphericFlightConditions > createAtmosphericFlightConditions(
        const std::shared_ptr< Body > bodyWithFlightConditions,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    // Check whether all required environment models are set.
    if( centralBody->getAtmosphereModel( ) == nullptr )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, body " + nameOfBodyExertingAcceleration +
                    " has no atmosphere model." );
    }

    if( centralBody->getShapeModel( ) == nullptr )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, body " + nameOfBodyExertingAcceleration +
                    " has no shape model." );
    }

    if( centralBody->getRotationalEphemeris( ) == nullptr )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, body " + nameOfBodyExertingAcceleration +
                    " has no rotation model." );
    }

    if( bodyWithFlightConditions->getAerodynamicCoefficientInterface( ) == nullptr )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, body " + nameOfBodyUndergoingAcceleration +
                    " has no aerodynamic coefficients." );
    }


    // Create aerodynamic angles calculator and set in flight conditions.
    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
            createAerodynamicAngleCalculator(
                bodyWithFlightConditions, centralBody, nameOfBodyUndergoingAcceleration,
                nameOfBodyExertingAcceleration );

    // Add wind model if present
    if( centralBody->getAtmosphereModel( )->getWindModel( ) != nullptr )
    {
        if( centralBody->getShapeModel( ) == nullptr )
        {
            std::cerr << "Warnning, body " << nameOfBodyExertingAcceleration << " has wind model, but no shape model, cannot compute wind as function of altitude " << std::endl;
        }
        else
        {
           aerodynamicAngleCalculator->setWindModel(
                        centralBody->getAtmosphereModel( )->getWindModel( ), centralBody->getShapeModel( ) );
        }
    }


    // Create flight conditions.
    std::function< double( const std::string& )> controlSurfaceDeflectionFunction;
    if( bodyWithFlightConditions->getVehicleSystems( ) != nullptr )
    {
        controlSurfaceDeflectionFunction = std::bind(
                    &system_models::VehicleSystems::getCurrentControlSurfaceDeflection,
                    bodyWithFlightConditions->getVehicleSystems( ), std::placeholders::_1 );
    }
    std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions =
            std::make_shared< aerodynamics::AtmosphericFlightConditions >(
                centralBody->getAtmosphereModel( ), centralBody->getShapeModel( ),
                bodyWithFlightConditions->getAerodynamicCoefficientInterface( ), aerodynamicAngleCalculator,
                controlSurfaceDeflectionFunction );

    return flightConditions;


}

//! Function to create a flight conditions object
std::shared_ptr< aerodynamics::FlightConditions >  createFlightConditions(
        const std::shared_ptr< Body > bodyWithFlightConditions,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    // Check whether all required environment models are set.
    if( centralBody->getShapeModel( ) == nullptr )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, central body " + nameOfBodyExertingAcceleration +
                    " has no shape model." );
    }

    if( centralBody->getRotationalEphemeris( ) == nullptr )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, central body " + nameOfBodyExertingAcceleration +
                    " has no rotation model." );
    }

    // Create aerodynamic angles calculator and set in flight conditions.
    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
            createAerodynamicAngleCalculator(
                bodyWithFlightConditions, centralBody, nameOfBodyUndergoingAcceleration,
                nameOfBodyExertingAcceleration );

    return std::make_shared< aerodynamics::FlightConditions >(
                centralBody->getShapeModel( ), aerodynamicAngleCalculator );

}


void addFlightConditions(
        const SystemOfBodies& bodies,
        const std::string bodyName,
        const std::string centralBodyName )
{
    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when creating flight conditions for body " + bodyName + ", body no found" );
    }

    if( bodies.count( centralBodyName ) == 0 )
    {
        throw std::runtime_error( "Error when creating flight conditions for body " + bodyName + ", central body " + centralBodyName + " not found" );
    }

    std::shared_ptr< Body > body = bodies.at( bodyName );
    std::shared_ptr< Body > centralBody = bodies.at( centralBodyName );
    if( centralBody->getAtmosphereModel( ) == nullptr )
    {
        body->setFlightConditions(
                createFlightConditions( body, centralBody, bodyName, centralBodyName ) );
    }
    else
    {
        body->setFlightConditions(
                createAtmosphericFlightConditions( body, centralBody, bodyName, centralBodyName ) );
    }
}

void addAtmosphericFlightConditions(
        const SystemOfBodies& bodies,
        const std::string bodyName,
        const std::string centralBodyName )
{
    addFlightConditions( bodies, bodyName, centralBodyName );
    if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                bodies.at( bodyName )->getFlightConditions( ) ) == nullptr )
    {
        if( bodies.at( centralBodyName )->getAtmosphereModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when adding atmospheric flight conditions for " + bodyName + " w.r.t. " +
                                      centralBodyName + ", conditions could not be created, central body has no atmosphere" );
        }
        else
        {
            throw std::runtime_error( "Error when adding atmospheric flight conditions for " + bodyName + " w.r.t. " +
                                      centralBodyName + ", conditions could not be created" );
        }

    }
}


void setGuidanceAnglesFunctions(
        const std::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const std::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator,
        const bool silenceWarnings  )
{
    utilities::printDeprecationError(
                "tudatpy.numerical_simulation.environment_setup.set_aerodynamic_guidance",
                "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#aerodynamic-guidance" );
}

void setGuidanceAnglesFunctions(
        const std::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const std::shared_ptr< simulation_setup::Body > bodyWithAngles,
        const bool silenceWarnings )
{
    utilities::printDeprecationError(
                "tudatpy.numerical_simulation.environment_setup.set_aerodynamic_guidance",
                "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#aerodynamic-guidance" );
}

void setAerodynamicOrientationFunctions(
        const std::shared_ptr< simulation_setup::Body > body,
        const std::function< double( ) > angleOfAttackFunction,
        const std::function< double( ) > angleOfSideslipFunction,
        const std::function< double( ) > bankAngleFunction,
        const std::function<void(const double)> updateFunction )
{
    utilities::printDeprecationError(
                "tudatpy.numerical_simulation.environment_setup.set_aerodynamic_orientation_functions",
                "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#aerodynamic-guidance" );
}

void setConstantAerodynamicOrientation(
        const std::shared_ptr< simulation_setup::Body > body,
        const double angleOfAttack,
        const double sideslipAngle,
        const double bankAngle,
        const bool silenceWarnings )
{
    utilities::printDeprecationError(
                "tudatpy.numerical_simulation.environment_setup.set_constant_aerodynamic_orientation",
                "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#aerodynamic-guidance" );
}

} // namespace simulation_setup

} // namespace tudat
