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

namespace tudat
{

namespace simulation_setup
{

//! Function to create an atmospheric flight conditions object
std::shared_ptr< aerodynamics::AtmosphericFlightConditions > createAtmosphericFlightConditions(
        const std::shared_ptr< Body > bodyWithFlightConditions,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::function< double( ) > angleOfAttackFunction,
        const std::function< double( ) > angleOfSideslipFunction,
        const std::function< double( ) > bankAngleFunction,
        const std::function< void( const double ) > angleUpdateFunction )
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
    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
            std::make_shared< reference_frames::AerodynamicAngleCalculator >(
                relativeBodyFixedStateFunction,
                std::bind( &simulation_setup::Body::getCurrentRotationToGlobalFrame, centralBody ),
                nameOfBodyExertingAcceleration, 1,
                angleOfAttackFunction, angleOfSideslipFunction, bankAngleFunction, angleUpdateFunction );

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
                    "Error when making flight conditions, body " + nameOfBodyExertingAcceleration +
                    " has no shape model." );
    }

    if( centralBody->getRotationalEphemeris( ) == nullptr )
    {
        throw std::runtime_error(
                    "Error when making flight conditions, body " + nameOfBodyExertingAcceleration +
                    " has no rotation model." );
    }

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
    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
            std::make_shared< reference_frames::AerodynamicAngleCalculator >(
                relativeBodyFixedStateFunction,
                std::bind( &simulation_setup::Body::getCurrentRotationToGlobalFrame, centralBody ),
                nameOfBodyExertingAcceleration, 1 );

    return std::make_shared< aerodynamics::FlightConditions >(
                centralBody->getShapeModel( ), aerodynamicAngleCalculator );

}


//! Function to set the angle of attack to trimmed conditions.
std::shared_ptr< aerodynamics::TrimOrientationCalculator > setTrimmedConditions(
        const std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions )
{
    // Create trim object.
    std::shared_ptr< aerodynamics::TrimOrientationCalculator > trimOrientation =
            std::make_shared< aerodynamics::TrimOrientationCalculator >(
                flightConditions->getAerodynamicCoefficientInterface( ) );

    // Create angle-of-attack function from trim object.
    std::function< std::vector< double >( ) > untrimmedIndependentVariablesFunction =
            std::bind( &aerodynamics::AtmosphericFlightConditions::getAerodynamicCoefficientIndependentVariables,
                         flightConditions );
    std::function< std::map< std::string, std::vector< double > >( ) > untrimmedControlSurfaceIndependentVariableFunction =
            std::bind( &aerodynamics::AtmosphericFlightConditions::getControlSurfaceAerodynamicCoefficientIndependentVariables,
                         flightConditions );

    flightConditions->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                std::bind( &aerodynamics::TrimOrientationCalculator::findTrimAngleOfAttackFromFunction, trimOrientation,
                             untrimmedIndependentVariablesFunction, untrimmedControlSurfaceIndependentVariableFunction ) );

    return trimOrientation;
}

//! Function to set the angle of attack to trimmed conditions.
std::shared_ptr< aerodynamics::TrimOrientationCalculator > setTrimmedConditions(
        const std::shared_ptr< Body > bodyWithFlightConditions )
{
    if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                bodyWithFlightConditions->getFlightConditions( ) ) == nullptr )
    {
        throw std::runtime_error( "Error, body does not have FlightConditions when setting trim conditions." );
    }

    return setTrimmedConditions( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                     bodyWithFlightConditions->getFlightConditions( ) ));
}

//! Function that must be called to link the AerodynamicGuidance object to the simulation
void setGuidanceAnglesFunctions(
        const std::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const std::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator,
        const bool silenceWarnings  )
{
    angleCalculator->setOrientationAngleFunctions(
                std::bind( &aerodynamics::AerodynamicGuidance::getCurrentAngleOfAttack, aerodynamicGuidance ),
                std::bind( &aerodynamics::AerodynamicGuidance::getCurrentAngleOfSideslip, aerodynamicGuidance ),
                std::bind( &aerodynamics::AerodynamicGuidance::getCurrentBankAngle, aerodynamicGuidance ),
                std::bind( &aerodynamics::AerodynamicGuidance::updateGuidance, aerodynamicGuidance, std::placeholders::_1 ),
                silenceWarnings );
}

//! Function that must be called to link the AerodynamicGuidance object to the simulation
void setGuidanceAnglesFunctions(
        const std::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const std::shared_ptr< simulation_setup::Body > bodyWithAngles,
        const bool silenceWarnings )
{
    std::shared_ptr< reference_frames::DependentOrientationCalculator >  orientationCalculator =
            bodyWithAngles->getDependentOrientationCalculator( );
    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator =
            std::dynamic_pointer_cast< reference_frames::AerodynamicAngleCalculator >( orientationCalculator );

    if( angleCalculator == nullptr )
    {
        throw std::runtime_error( "Error, body does not have AerodynamicAngleCalculator when setting aerodynamic guidance" );
    }
    else
    {
        setGuidanceAnglesFunctions( aerodynamicGuidance, angleCalculator, silenceWarnings );
    }
}

void setAerodynamicOrientationFunctions(
        const std::shared_ptr< simulation_setup::Body > body,
        const std::function< double( ) > angleOfAttackFunction,
        const std::function< double( ) > angleOfSideslipFunction,
        const std::function< double( ) > bankAngleFunction,
        const std::function< void( const double ) > angleUpdateFunction )
{
    if( body->getFlightConditions( ) == nullptr )
    {
        throw std::runtime_error( "Error when setting aerodynamic angle functions, body " + body->getBodyName( ) + " has no FlightConditions" );
    }

    std::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions =
            body->getFlightConditions( );
    vehicleFlightConditions->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                angleOfAttackFunction, angleOfSideslipFunction, bankAngleFunction, angleUpdateFunction );
}

void setConstantAerodynamicOrientation(
        const std::shared_ptr< simulation_setup::Body > body,
        const double angleOfAttack,
        const double sideslipAngle,
        const double bankAngle,
        const bool silenceWarnings )
{
    if( body->getFlightConditions( ) == nullptr )
    {
        throw std::runtime_error( "Error when setting constant aerodynamic angles, body " + body->getBodyName( ) + " has no FlightConditions" );
    }

    std::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions =
            body->getFlightConditions( );
    vehicleFlightConditions->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                angleOfAttack, sideslipAngle, bankAngle, silenceWarnings );
}

} // namespace simulation_setup

} // namespace tudat
