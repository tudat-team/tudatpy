/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"

namespace tudat
{

namespace simulation_setup
{

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
                    "Error when making flight conditions, body " + nameOfBodyUndergoingAcceleration +
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

    // Add wind model if present
    if( centralBody->getAtmosphereModel( )->getWindModel( ) != NULL )
    {
        if( centralBody->getShapeModel( ) == NULL )
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
