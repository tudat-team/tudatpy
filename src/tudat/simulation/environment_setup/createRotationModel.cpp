/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/astro/ephemerides/iauRotationModel.h"
#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"

#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"
#include "tudat/astro/ephemerides/customRotationalEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

std::shared_ptr< ephemerides::RotationalEphemeris > createGcrsToItrsRotationModel(
        const std::shared_ptr< GcrsToItrsRotationModelSettings >& gcrsToItrsRotationSettings );

std::shared_ptr< RotationModelSettings > simpleRotationModelFromSpiceSettings( const std::string& originalFrame,
                                                                               const std::string& targetFrame,
                                                                               const std::string& targetFrameSpice,
                                                                               const double initialTime )
{
    return std::make_shared< SimpleRotationModelSettings >(
            originalFrame,
            targetFrame,
            spice_interface::computeRotationQuaternionBetweenFrames( originalFrame, targetFrameSpice, initialTime ),
            initialTime,
            spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame( originalFrame, targetFrameSpice, initialTime ).norm( ) );
}

std::function< Eigen::Matrix3d( const double ) > getRotationFunctionFromSatelliteBasedFrame(
        const ephemerides::SatelliteBasedFrames frameId,
        const SystemOfBodies& bodies,
        const std::string& body,
        const std::string& centralBody )
{
    std::function< Eigen::Matrix3d( const double ) > rotationToInertialFrameFunction;
    switch( frameId )
    {
        case ephemerides::inertial_satellite_based_frame: {
            rotationToInertialFrameFunction = []( const double ) { return Eigen::Matrix3d::Identity( ); };
            break;
        }
        case ephemerides::tnw_satellite_based_frame: {
            if( bodies.count( body ) == 0 )
            {
                throw std::runtime_error( "Error when getting rotation function from TNW, input body " + body + " not found" );
            }

            if( ephemerides::isFrameInertial( centralBody ) )
            {
                rotationToInertialFrameFunction = [ = ]( const double ) {
                    return reference_frames::getTnwToInertialRotation( bodies.at( body )->getState( ), true );
                };
            }
            else
            {
                if( bodies.count( centralBody ) == 0 )
                {
                    throw std::runtime_error( "Error when getting rotation function from TNW, input central body " + centralBody +
                                              " not found" );
                }
                rotationToInertialFrameFunction = [ = ]( const double ) {
                    return reference_frames::getTnwToInertialRotation(
                            bodies.at( body )->getState( ) - bodies.at( centralBody )->getState( ), true );
                };
            }
            break;
        }
        case ephemerides::rsw_satellite_based_frame: {
            if( bodies.count( body ) == 0 )
            {
                throw std::runtime_error( "Error when getting rotation function from TNW, input body " + body + " not found" );
            }

            if( ephemerides::isFrameInertial( centralBody ) )
            {
                rotationToInertialFrameFunction = [ = ]( const double ) {
                    return reference_frames::getRswSatelliteCenteredToInertialFrameRotationMatrix( bodies.at( body )->getState( ) );
                };
            }
            else
            {
                if( bodies.count( centralBody ) == 0 )
                {
                    throw std::runtime_error( "Error when getting rotation function from TNW, input central body " + centralBody +
                                              " not found" );
                }
                rotationToInertialFrameFunction = [ = ]( const double ) {
                    return reference_frames::getRswSatelliteCenteredToInertialFrameRotationMatrix( bodies.at( body )->getState( ) -
                                                                                                   bodies.at( centralBody )->getState( ) );
                };
            }
            break;
        }
        default:
            throw std::runtime_error( "Error when getting rotation function from satellite fixedd frame, input type " +
                                      std::to_string( frameId ) + " not found" );
    }
    return rotationToInertialFrameFunction;
}

//! Function to retrieve a state from one of two functions
Eigen::Vector6d getStateFromSelectedStateFunction( const double currentTime,
                                                   const bool useFirstFunction,
                                                   const std::function< Eigen::Vector6d( const double ) > stateFunction1,
                                                   const std::function< Eigen::Vector6d( const double ) > stateFunction2 )
{
    return ( useFirstFunction ) ? ( stateFunction1( currentTime ) ) : ( stateFunction2( currentTime ) );
}

//! Function to create a state function for a body, valid both during propagation, and outside propagation
std::function< Eigen::Vector6d( const double, bool ) > createRelativeStateFunction( const SystemOfBodies& bodies,
                                                                                    const std::string orbitingBody,
                                                                                    const std::string centralBody )
{
    // Retrieve state functions for relevant bodies (obtained from current state of body objects)
    std::function< Eigen::Vector6d( const double ) > bodyInertialStateFunction = std::bind( &Body::getState, bodies.at( orbitingBody ) );
    std::function< Eigen::Vector6d( const double ) > centralBodyInertialStateFunction =
            std::bind( &Body::getState, bodies.at( centralBody ) );

    // Define relative state function from body object
    std::function< Eigen::Vector6d( const double ) > fromBodyStateFunction = std::bind(
            &ephemerides::getDifferenceBetweenStates, bodyInertialStateFunction, centralBodyInertialStateFunction, std::placeholders::_1 );

    // Define state function from ephemeris
    std::function< Eigen::Vector6d( const double ) > fromEphemerisStateFunction;

    if( bodies.at( orbitingBody )->getEphemeris( )->getReferenceFrameOrigin( ) == centralBody )
    {
        fromEphemerisStateFunction =
                std::bind( &ephemerides::Ephemeris::getCartesianState, bodies.at( orbitingBody )->getEphemeris( ), std::placeholders::_1 );
    }
    else
    {
        std::function< Eigen::Vector6d( const double ) > ephemerisInertialStateFunction =
                std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, bodies.at( orbitingBody ), std::placeholders::_1 );
        std::function< Eigen::Vector6d( const double ) > ephemerisCentralBodyInertialStateFunction =
                std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, bodies.at( centralBody ), std::placeholders::_1 );
        fromEphemerisStateFunction = std::bind( &ephemerides::getDifferenceBetweenStates,
                                                ephemerisInertialStateFunction,
                                                ephemerisCentralBodyInertialStateFunction,
                                                std::placeholders::_1 );
    }

    return std::bind( &getStateFromSelectedStateFunction,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      fromBodyStateFunction,
                      fromEphemerisStateFunction );
}

//! Function to set the angle of attack to trimmed conditions.
std::shared_ptr< aerodynamics::TrimOrientationCalculator > createTrimCalculator( const std::shared_ptr< Body > bodyWithFlightConditions )
{
    if( bodyWithFlightConditions->getAerodynamicCoefficientInterface( ) == nullptr )
    {
        throw std::runtime_error( "Error, body does not have AerodynamicCoefficientInterface when creating trim conditions." );
    }

    // Create trim object.
    return std::make_shared< aerodynamics::TrimOrientationCalculator >( bodyWithFlightConditions->getAerodynamicCoefficientInterface( ) );
}

//! Function to set the angle of attack to trimmed conditions.
void linkTrimmedConditions( const std::shared_ptr< aerodynamics::TrimOrientationCalculator > trimCalculator,
                            const std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions,
                            const std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > rotationModel,
                            const std::function< Eigen::Vector2d( const double ) > sideslipAndBankAngleFunction = nullptr )
{
    // Create angle-of-attack function from trim object.
    std::function< std::vector< double >( ) > untrimmedIndependentVariablesFunction =
            std::bind( &aerodynamics::AtmosphericFlightConditions::getAerodynamicCoefficientIndependentVariables, flightConditions );
    std::function< std::map< std::string, std::vector< double > >( ) > untrimmedControlSurfaceIndependentVariableFunction = std::bind(
            &aerodynamics::AtmosphericFlightConditions::getControlSurfaceAerodynamicCoefficientIndependentVariables, flightConditions );

    std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction = [ = ]( const double currentTime ) {
        Eigen::Vector2d sideslipBankAngles = Eigen::Vector2d::Zero( );
        if( sideslipAndBankAngleFunction != nullptr )
        {
            sideslipBankAngles = sideslipAndBankAngleFunction( currentTime );
            flightConditions->getAerodynamicAngleCalculator( )->setSideslipAndBankAngles( sideslipBankAngles( 0 ),
                                                                                          sideslipBankAngles( 1 ) );
        }

        double angleOfAttack = trimCalculator->findTrimAngleOfAttackFromFunction( untrimmedIndependentVariablesFunction,
                                                                                  untrimmedControlSurfaceIndependentVariableFunction );
        return ( Eigen::Vector3d( ) << angleOfAttack, sideslipBankAngles( 0 ), sideslipBankAngles( 1 ) ).finished( );
    };
    rotationModel->setAerodynamicAngleFunction( aerodynamicAngleFunction );
}

std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > createAerodynamicAngleBasedRotationModel(
        const std::string& body,
        const std::string& centralBody,
        const SystemOfBodies& bodies,
        const std::string& originalFrame,
        const std::string& targetFrame )
{
    if( bodies.count( body ) == 0 )
    {
        throw std::runtime_error( "Error when creating aerodynamic angle based rotation model, body " + body + " not found" );
    }

    std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > rotationModel =
            std::make_shared< ephemerides::AerodynamicAngleRotationalEphemeris >( centralBody, originalFrame, targetFrame );
    applyEnvironmentModelClosure( rotationModel, body, bodies );
    return rotationModel;
}

std::shared_ptr< ephemerides::InertialBodyFixedDirectionCalculator > createInertialDirectionCalculator(
        const std::shared_ptr< InertialDirectionSettings > directionSettings,
        const std::string& body )
{
    return createInertialDirectionCalculator( directionSettings, body, SystemOfBodies( ) );
}

std::shared_ptr< ephemerides::InertialBodyFixedDirectionCalculator > createInertialDirectionCalculatorWithoutClosure(
        const std::shared_ptr< InertialDirectionSettings > directionSettings,
        const std::string& body )
{
    std::shared_ptr< ephemerides::InertialBodyFixedDirectionCalculator > directionCalculator;
    switch( directionSettings->inertialDirectionType_ )
    {
        case custom_inertial_direction: {
            std::shared_ptr< CustomInertialDirectionSettings > customDirectionSettings =
                    std::dynamic_pointer_cast< CustomInertialDirectionSettings >( directionSettings );
            if( customDirectionSettings == nullptr )
            {
                throw std::runtime_error(
                        "Error when making direction calculator for direction-based rotation model, expected type "
                        "CustomInertialDirectionSettings" );
            }
            directionCalculator = std::make_shared< ephemerides::CustomBodyFixedDirectionCalculator >(
                    customDirectionSettings->inertialBodyAxisDirectionFunction_ );
            break;
        }
        case state_based_inertial_direction: {
            std::shared_ptr< StateBasedInertialDirectionSettings > stateBasedDirectionSettings =
                    std::dynamic_pointer_cast< StateBasedInertialDirectionSettings >( directionSettings );
            if( stateBasedDirectionSettings == nullptr )
            {
                throw std::runtime_error(
                        "Error when making direction calculator for direction-based rotation model, expected type "
                        "StateBasedInertialDirectionSettings" );
            }
            directionCalculator = std::make_shared< ephemerides::StateBasedBodyFixedDirectionCalculator >(
                    stateBasedDirectionSettings->centralBody_,
                    stateBasedDirectionSettings->isColinearWithVelocity_,
                    stateBasedDirectionSettings->directionIsOppositeToVector_,
                    nullptr );
            break;
        }
        case bilinear_tangent_inertial_direction: {
            throw std::runtime_error(
                    "Error when making direction calculator for direction-based rotation model, bilinear_tangent_inertial_direction not "
                    "yet implemented" );
            break;
        }
        default: {
            throw std::runtime_error(
                    "Error when making direction calculator for direction-based rotation model, type not yet implemented" );
            break;
        }
    }
    directionCalculator->setDirectionFrame( directionSettings->directionFrame_ );
    return directionCalculator;
}

void applyEnvironmentModelClosure( const std::shared_ptr< ephemerides::InertialBodyFixedDirectionCalculator > directionCalculator,
                                   const std::string& body,
                                   const SystemOfBodies& bodies )
{
    // Retrieve rotation from frame in which direction is defined (if not inertial)
    std::pair< ephemerides::SatelliteBasedFrames, std::string > directionFrame = directionCalculator->getDirectionFrame( );
    std::function< Eigen::Matrix3d( const double ) > rotationFunction = nullptr;
    if( directionFrame.first != ephemerides::inertial_satellite_based_frame )
    {
        rotationFunction = getRotationFunctionFromSatelliteBasedFrame( directionFrame.first, bodies, body, directionFrame.second );
    }

    std::shared_ptr< ephemerides::StateBasedBodyFixedDirectionCalculator > stateBasedDirectionCalculator =
            std::dynamic_pointer_cast< ephemerides::StateBasedBodyFixedDirectionCalculator >( directionCalculator );
    if( stateBasedDirectionCalculator != nullptr )
    {
        // Retrieve state function of body for which thrust is to be computed.
        std::function< Eigen::Vector6d( ) > bodyStateFunction = std::bind( &Body::getState, bodies.at( body ) );
        std::function< Eigen::Vector6d( ) > centralBodyStateFunction;

        // Retrieve state function of central body (or set to zero if inertial)
        std::string centralBody = stateBasedDirectionCalculator->getCentralBody( );
        if( centralBody != "SSB" )
        {
            centralBodyStateFunction = std::bind( &Body::getState, bodies.at( centralBody ) );
        }
        else if( bodies.getFrameOrigin( ) == "SSB" )
        {
            centralBodyStateFunction = []( ) { return Eigen::Vector6d::Zero( ); };
        }
        else
        {
            throw std::runtime_error(
                    "Error when getting state-direction-based rotation model, requested state w.r.t. SSB, but SSB is not the global "
                    "origin" );
        }

        // Define relative state function
        stateBasedDirectionCalculator->applyEnvironmentModelClosure(
                std::bind( &ephemerides::getRelativeState, std::placeholders::_1, bodyStateFunction, centralBodyStateFunction ),
                rotationFunction );
    }
    else
    {
        directionCalculator->applyEnvironmentModelClosure( rotationFunction );
    }
}

std::shared_ptr< ephemerides::InertialBodyFixedDirectionCalculator > createInertialDirectionCalculator(
        const std::shared_ptr< InertialDirectionSettings > directionSettings,
        const std::string& body,
        const SystemOfBodies& bodies )
{
    std::shared_ptr< ephemerides::InertialBodyFixedDirectionCalculator > directionCalculator =
            createInertialDirectionCalculatorWithoutClosure( directionSettings, body );
    applyEnvironmentModelClosure( directionCalculator, body, bodies );
    return directionCalculator;
}

std::shared_ptr< ephemerides::DirectionBasedRotationalEphemeris > createStateDirectionBasedRotationModel(
        const std::string& body,
        const std::string& centralBody,
        const SystemOfBodies& bodies,
        const Eigen::Vector3d& associatedBodyFixedDirection,
        const std::string& originalFrame,
        const std::string& targetFrame,
        const bool isColinearWithVelocity,
        const bool directionIsOppositeToVector,
        const std::function< double( const double ) > freeRotationAngleFunction )
{
    std::shared_ptr< ephemerides::InertialBodyFixedDirectionCalculator > directionCalculator = createInertialDirectionCalculator(
            std::make_shared< StateBasedInertialDirectionSettings >( centralBody, isColinearWithVelocity, directionIsOppositeToVector ),
            body,
            bodies );
    return std::make_shared< ephemerides::DirectionBasedRotationalEphemeris >(
            directionCalculator, associatedBodyFixedDirection, originalFrame, targetFrame, freeRotationAngleFunction );
}

std::shared_ptr< ephemerides::RotationalEphemeris > createTrimmedAerodynamicAngleBasedRotationModel(
        const std::string& body,
        const std::string& centralBody,
        const SystemOfBodies& bodies,
        const std::string& originalFrame,
        const std::string& targetFrame,
        const std::function< Eigen::Vector2d( const double ) > sideslipAndBankAngleFunction )
{
    std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > rotationModel =
            std::make_shared< ephemerides::AerodynamicAngleRotationalEphemeris >( centralBody, originalFrame, targetFrame );
    rotationModel->setIsPitchTrimmed( true, sideslipAndBankAngleFunction );

    std::shared_ptr< Body > trimmedBody = bodies.at( body );
    trimmedBody->setRotationalEphemeris( rotationModel );

    // Close aerodynamic angle calculator, then flight conditions (using this calculator), then trim conditions
    applyEnvironmentModelClosure( rotationModel, body, bodies );
    if( trimmedBody->getFlightConditions( ) == nullptr )
    {
        addFlightConditions( bodies, body, centralBody );
    }
    else if( trimmedBody->getFlightConditions( )->getAerodynamicAngleCalculator( ) == nullptr )
    {
        applyEnvironmentModelClosure( trimmedBody->getFlightConditions( ), body, bodies );
    }

    if( !applyEnvironmentModelClosure( rotationModel, body, bodies ) )
    {
        throw std::runtime_error( "Error, body does not have FlightConditions when setting trim conditions." );
    }
    return rotationModel;
}

bool applyEnvironmentModelClosure( const std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel,
                                   const std::string& body,
                                   const SystemOfBodies& bodies )
{
    using namespace ephemerides;

    std::shared_ptr< SynchronousRotationalEphemeris > synchronousRotationModel =
            std::dynamic_pointer_cast< SynchronousRotationalEphemeris >( rotationModel );
    if( synchronousRotationModel != nullptr )
    {
        synchronousRotationModel->applyEnvironmentModelClosure(
                createRelativeStateFunction( bodies, body, synchronousRotationModel->getCentralBodyName( ) ) );
        return true;
    }

    std::shared_ptr< AerodynamicAngleRotationalEphemeris > aerodynamicAngleRotationModel =
            std::dynamic_pointer_cast< AerodynamicAngleRotationalEphemeris >( rotationModel );
    if( aerodynamicAngleRotationModel != nullptr )
    {
        std::string centralBody = aerodynamicAngleRotationModel->getCentralBodyName( );
        if( aerodynamicAngleRotationModel->getAerodynamicAngleCalculator( ) == nullptr )
        {
            // Use aerodynamic angle calculator of existing flight conditions, or create new one
            std::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator;
            if( bodies.at( body )->getFlightConditions( ) != nullptr &&
                bodies.at( body )->getFlightConditions( )->getAerodynamicAngleCalculator( ) != nullptr )
            {
                angleCalculator = bodies.at( body )->getFlightConditions( )->getAerodynamicAngleCalculator( );
            }
            else
            {
                if( bodies.count( centralBody ) == 0 )
                {
                    throw std::runtime_error( "Error when creating aerodynamic angle based rotation model, central body " + centralBody +
                                              " not found" );
                }
                angleCalculator =
                        createAerodynamicAngleCalculatorFromBodyStates( bodies.at( body ), bodies.at( centralBody ), centralBody );
            }
            aerodynamicAngleRotationModel->applyEnvironmentModelClosure( angleCalculator );
            angleCalculator->setBodyFixedAngleInterface(
                    std::make_shared< reference_frames::FromAeroEphemerisAerodynamicAngleInterface >( aerodynamicAngleRotationModel ) );
        }

        // Trimmed angle of attack requires atmospheric flight conditions to be closed
        if( aerodynamicAngleRotationModel->getIsPitchTrimmed( ) )
        {
            std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions =
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >( bodies.at( body )->getFlightConditions( ) );
            if( flightConditions == nullptr || flightConditions->getAerodynamicAngleCalculator( ) == nullptr )
            {
                return false;
            }
            linkTrimmedConditions( createTrimCalculator( bodies.at( body ) ),
                                   flightConditions,
                                   aerodynamicAngleRotationModel,
                                   aerodynamicAngleRotationModel->getTrimSideslipAndBankAngleFunction( ) );
        }
        return true;
    }

    std::shared_ptr< DirectionBasedRotationalEphemeris > directionBasedRotationModel =
            std::dynamic_pointer_cast< DirectionBasedRotationalEphemeris >( rotationModel );
    if( directionBasedRotationModel != nullptr )
    {
        applyEnvironmentModelClosure( directionBasedRotationModel->getInertialBodyAxisDirectionCalculator( ), body, bodies );
        return true;
    }

    // Other rotation models require no closure
    return true;
}

//! Function to create a rotation model.
std::shared_ptr< ephemerides::RotationalEphemeris > createRotationModel(
        const std::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body )
{
    return createRotationModel( rotationModelSettings, body, SystemOfBodies( ) );
}

//! Function to create a rotation model (with closure applied).
std::shared_ptr< ephemerides::RotationalEphemeris > createRotationModel(
        const std::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body,
        const SystemOfBodies& bodies )
{
    std::shared_ptr< PitchTrimRotationSettings > trimRotationSettings =
            std::dynamic_pointer_cast< PitchTrimRotationSettings >( rotationModelSettings );
    if( trimRotationSettings != nullptr )
    {
        if( rotationModelSettings->getOriginalFrame( ) == "" )
        {
            rotationModelSettings->resetOriginalFrame( bodies.getFrameOrientation( ) );
        }
        return createTrimmedAerodynamicAngleBasedRotationModel( body,
                                                                trimRotationSettings->centralBody_,
                                                                bodies,
                                                                trimRotationSettings->getOriginalFrame( ),
                                                                trimRotationSettings->getTargetFrame( ),
                                                                trimRotationSettings->sideslipAndBankAngleFunction_ );
    }

    std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel =
            createRotationModelWithoutClosure( rotationModelSettings, body, bodies );
    if( !applyEnvironmentModelClosure( rotationModel, body, bodies ) )
    {
        throw std::runtime_error( "Error when creating rotation model of body " + body + ", closure could not be applied" );
    }
    return rotationModel;
}

//! Function to create a rotation model, without closure (links to other environment models)
std::shared_ptr< ephemerides::RotationalEphemeris > createRotationModelWithoutClosure(
        const std::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body,
        const SystemOfBodies& bodies )
{
    if( rotationModelSettings->getOriginalFrame( ) == "" )
    {
        rotationModelSettings->resetOriginalFrame( bodies.getFrameOrientation( ) );
    }

    using namespace tudat::ephemerides;

    // Declare return object.
    std::shared_ptr< RotationalEphemeris > rotationalEphemeris;

    // Check which type of rotation model is to be created.
    switch( rotationModelSettings->getRotationType( ) )
    {
        case simple_rotation_model: {
            // Check whether settings for simple rotation model are consistent with its type.
            std::shared_ptr< SimpleRotationModelSettings > simpleRotationSettings =
                    std::dynamic_pointer_cast< SimpleRotationModelSettings >( rotationModelSettings );
            if( simpleRotationSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected simple rotation model settings for " + body );
            }
            else
            {
                // Create and initialize simple rotation model.
                rotationalEphemeris = std::make_shared< SimpleRotationalEphemeris >( simpleRotationSettings->getInitialOrientation( ),
                                                                                     simpleRotationSettings->getRotationRate( ),
                                                                                     simpleRotationSettings->getInitialTime( ),
                                                                                     simpleRotationSettings->getOriginalFrame( ),
                                                                                     simpleRotationSettings->getTargetFrame( ) );
            }
            break;
        }
            // #if TUDAT_BUILD_WITH_SOFA_INTERFACE
        case gcrs_to_itrs_rotation_model: {
            // Check whether settings for simple rotation model are consistent with its type.
            std::shared_ptr< GcrsToItrsRotationModelSettings > gcrsToItrsRotationSettings =
                    std::dynamic_pointer_cast< GcrsToItrsRotationModelSettings >( rotationModelSettings );
            if( gcrsToItrsRotationSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected GCRS to ITRS rotation model settings for " + body );
            }
            else
            {
                rotationalEphemeris = createGcrsToItrsRotationModel( gcrsToItrsRotationSettings );
                break;
            }
        }
            // #endif

        case spice_rotation_model: {
            std::shared_ptr< SpiceRotationModelSettings > spiceRotationModelSettings =
                    std::dynamic_pointer_cast< SpiceRotationModelSettings >( rotationModelSettings );
            std::string spiceFrameName = "";
            if( spiceRotationModelSettings != nullptr )
            {
                spiceFrameName = spiceRotationModelSettings->getSpiceFrameName( );
            }
            // Create rotational ephemeris directly from Spice.
            rotationalEphemeris = std::make_shared< SpiceRotationalEphemeris >(
                    rotationModelSettings->getOriginalFrame( ), rotationModelSettings->getTargetFrame( ), spiceFrameName );
            break;
        }
        case planetary_rotation_model: {
            std::shared_ptr< PlanetaryRotationModelSettings > planetaryRotationModelSettings =
                    std::dynamic_pointer_cast< PlanetaryRotationModelSettings >( rotationModelSettings );
            if( planetaryRotationModelSettings == nullptr )
            {
                std::cerr << "Error, expected planetary rotation model settings for " << body << std::endl;
            }
            else
            {
                Eigen::VectorXd bodyKeplerElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                        spice_interface::getBodyCartesianStateAtEpoch(
                                body, planetaryRotationModelSettings->getCentralBody( ), "ECLIPJ2000", "NONE", 0.0 ),
                        spice_interface::getBodyGravitationalParameter( planetaryRotationModelSettings->getCentralBody( ) ) );

                double meanAnomalyAtJ2000 = orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
                        orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly( bodyKeplerElements( 5 ),
                                                                                           bodyKeplerElements( 1 ) ),
                        bodyKeplerElements( 1 ) );

                double meanMotion = orbital_element_conversions::convertSemiMajorAxisToEllipticalMeanMotion(
                        bodyKeplerElements( 0 ),
                        spice_interface::getBodyGravitationalParameter( planetaryRotationModelSettings->getCentralBody( ) ) +
                                spice_interface::getBodyGravitationalParameter( body ) );

                std::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationAnglesCalculator =
                        std::make_shared< PlanetaryOrientationAngleCalculator >(
                                planetaryRotationModelSettings->getAnglePsiAtEpoch( ),
                                planetaryRotationModelSettings->getAnglePsiRateAtEpoch( ),
                                planetaryRotationModelSettings->getAngleIAtEpoch( ),
                                planetaryRotationModelSettings->getAngleIRateAtEpoch( ),
                                planetaryRotationModelSettings->getAnglePhiAtEpoch( ),
                                planetaryRotationModelSettings->getAnglePhiRateAtEpoch( ),
                                planetaryRotationModelSettings->getCoreFactor( ),
                                planetaryRotationModelSettings->getFreeCoreNutationRate( ),
                                meanMotion,
                                meanAnomalyAtJ2000,
                                planetaryRotationModelSettings->getOriginalFrame( ),
                                planetaryRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                                planetaryRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                                planetaryRotationModelSettings->getTimeDependentPhaseCorrectionFunctions( ),
                                planetaryRotationModelSettings->getRotationRateCorrections( ),
                                planetaryRotationModelSettings->getxPolarMotionCoefficients( ),
                                planetaryRotationModelSettings->getyPolarMotionCoefficients( ) );

                rotationalEphemeris = std::make_shared< PlanetaryRotationModel >( planetaryRotationModelSettings->getAngleN( ),
                                                                                  planetaryRotationModelSettings->getAngleJ( ),
                                                                                  planetaryOrientationAnglesCalculator,
                                                                                  planetaryRotationModelSettings->getOriginalFrame( ),
                                                                                  planetaryRotationModelSettings->getTargetFrame( ) );
            }

            break;
        }
        case iau_rotation_model: {
            std::shared_ptr< IauRotationModelSettings > iauRotationModelSettings =
                    std::dynamic_pointer_cast< IauRotationModelSettings >( rotationModelSettings );
            if( iauRotationModelSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected IAU rotation model settings for " + body );
            }
            else
            {
                rotationalEphemeris = std::make_shared< IauRotationModel >( iauRotationModelSettings->getOriginalFrame( ),
                                                                            iauRotationModelSettings->getTargetFrame( ),
                                                                            iauRotationModelSettings->nominalMeridian_,
                                                                            iauRotationModelSettings->nominalPole_,
                                                                            iauRotationModelSettings->rotationRate_,
                                                                            iauRotationModelSettings->polePrecession_,
                                                                            iauRotationModelSettings->meridianPeriodicTerms_,
                                                                            iauRotationModelSettings->polePeriodicTerms_,
                                                                            iauRotationModelSettings->referenceEpochJ2000_,
                                                                            iauRotationModelSettings->anglesBaseFrame_ );
            }

            break;
        }
        case synchronous_rotation_model: {
            std::shared_ptr< SynchronousRotationModelSettings > synchronousRotationSettings =
                    std::dynamic_pointer_cast< SynchronousRotationModelSettings >( rotationModelSettings );
            if( synchronousRotationSettings == NULL )
            {
                throw std::runtime_error( "Error, expected synchronous rotation model settings for " + body );
            }
            else
            {
                if( bodies.count( synchronousRotationSettings->getCentralBodyName( ) ) == 0 )
                {
                    throw std::runtime_error( "Error, when making synchronous rotation model for " + body + ", central body " +
                                              synchronousRotationSettings->getCentralBodyName( ) + " not found." );
                }
                if( bodies.at( body )->getEphemeris( ) == nullptr )
                {
                    throw std::runtime_error( "Error, when making synchronous rotation model for " + body +
                                              ", central body, no ephemeris found." );
                }
                if( bodies.at( body )->getEphemeris( )->getReferenceFrameOrigin( ) == synchronousRotationSettings->getCentralBodyName( ) )
                {
                    if( bodies.at( body )->getEphemeris( )->getReferenceFrameOrientation( ) !=
                        synchronousRotationSettings->getOriginalFrame( ) )
                    {
                        throw std::runtime_error( "Error, ephemeris of body " + body + " is in " +
                                                  bodies.at( body )->getEphemeris( )->getReferenceFrameOrientation( ) +
                                                  " frame when making synchronous rotation model, expected " +
                                                  synchronousRotationSettings->getOriginalFrame( ) + " frame." );
                    }
                }
                std::shared_ptr< SynchronousRotationalEphemeris > synchronousRotationalEphemeris =
                        std::make_shared< SynchronousRotationalEphemeris >( nullptr,
                                                                            synchronousRotationSettings->getCentralBodyName( ),
                                                                            synchronousRotationSettings->getOriginalFrame( ),
                                                                            synchronousRotationSettings->getTargetFrame( ) );

                rotationalEphemeris = synchronousRotationalEphemeris;
            }
            break;
        }
        case tabulated_rotation_model: {
            // Check whether settings for simple rotation model are consistent with its type.
            std::shared_ptr< TabulatedRotationSettings > tabulatedRotationSettings =
                    std::dynamic_pointer_cast< TabulatedRotationSettings >( rotationModelSettings );
            if( tabulatedRotationSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected tabulated rotation model settings for " + body );
            }
            else
            {
                // Create and initialize simple rotation model.
                rotationalEphemeris = std::make_shared< TabulatedRotationalEphemeris< double, double > >(
                        interpolators::createOneDimensionalInterpolator( tabulatedRotationSettings->getBodyStateHistory( ),
                                                                         tabulatedRotationSettings->getInterpolatorSettings( ) ),
                        tabulatedRotationSettings->getOriginalFrame( ),
                        tabulatedRotationSettings->getTargetFrame( ) );
            }
            break;
        }
        case aerodynamic_angle_based_rotation_model: {
            // Check whether settings for simple rotation model are consistent with its type.
            std::shared_ptr< AerodynamicAngleRotationSettings > aerodynamicAngleRotationSettings =
                    std::dynamic_pointer_cast< AerodynamicAngleRotationSettings >( rotationModelSettings );
            if( aerodynamicAngleRotationSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected aerdynamic angle based rotation model settings for " + body );
            }
            else
            {
                rotationalEphemeris = std::make_shared< AerodynamicAngleRotationalEphemeris >(
                        aerodynamicAngleRotationSettings->centralBody_,
                        aerodynamicAngleRotationSettings->getOriginalFrame( ),
                        aerodynamicAngleRotationSettings->getTargetFrame( ),
                        aerodynamicAngleRotationSettings->aerodynamicAngleFunction_ );
            }
            break;
        }
        case pitch_trim_rotation_model: {
            // Check whether settings for simple rotation model are consistent with its type.
            std::shared_ptr< PitchTrimRotationSettings > trimRotationSettings =
                    std::dynamic_pointer_cast< PitchTrimRotationSettings >( rotationModelSettings );
            if( trimRotationSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected pitch trim rotation model settings for " + body );
            }
            else
            {
                std::shared_ptr< AerodynamicAngleRotationalEphemeris > trimmedRotationalEphemeris =
                        std::make_shared< AerodynamicAngleRotationalEphemeris >( trimRotationSettings->centralBody_,
                                                                                 trimRotationSettings->getOriginalFrame( ),
                                                                                 trimRotationSettings->getTargetFrame( ) );
                trimmedRotationalEphemeris->setIsPitchTrimmed( true, trimRotationSettings->sideslipAndBankAngleFunction_ );
                rotationalEphemeris = trimmedRotationalEphemeris;

                // Trimmed conditions require flight conditions (closed together with the rotation model)
                if( bodies.count( body ) > 0 && bodies.at( body )->getFlightConditions( ) == nullptr )
                {
                    bodies.at( body )->setFlightConditions(
                            createFlightConditionsWithoutClosure( bodies, body, trimRotationSettings->centralBody_ ) );
                }
            }
            break;
        }
        case body_fixed_direction_based_rotation_model: {
            // Check whether settings for simple rotation model are consistent with its type.
            std::shared_ptr< BodyFixedDirectionBasedRotationSettings > bodyFixedDirectionBasedRotationSettings =
                    std::dynamic_pointer_cast< BodyFixedDirectionBasedRotationSettings >( rotationModelSettings );
            if( bodyFixedDirectionBasedRotationSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected body-fixed direction rotation model settings for " + body );
            }
            else
            {
                std::shared_ptr< InertialBodyFixedDirectionCalculator > directionCalculator =
                        createInertialDirectionCalculatorWithoutClosure(
                                bodyFixedDirectionBasedRotationSettings->inertialDirectionSettings_, body );

                // Create and initialize simple rotation model.
                rotationalEphemeris = std::make_shared< DirectionBasedRotationalEphemeris >(
                        directionCalculator,
                        Eigen::Vector3d::UnitX( ),
                        bodyFixedDirectionBasedRotationSettings->getOriginalFrame( ),
                        bodyFixedDirectionBasedRotationSettings->getTargetFrame( ),
                        bodyFixedDirectionBasedRotationSettings->freeRotationAngleFunction_ );
            }
            break;
        }
        case orbital_state_based_rotation_model: {
            // Check whether settings for simple rotation model are consistent with its type.
            std::shared_ptr< OrbitalStateBasedRotationSettings > orbitalStateBasedRotationSettings =
                    std::dynamic_pointer_cast< OrbitalStateBasedRotationSettings >( rotationModelSettings );
            if( orbitalStateBasedRotationSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected orbital state based rotation model settings for " + body );
            }
            else
            {
                rotationalEphemeris = std::make_shared< DirectionBasedRotationalEphemeris >(
                        createInertialDirectionCalculatorWithoutClosure(
                                std::make_shared< StateBasedInertialDirectionSettings >(
                                        orbitalStateBasedRotationSettings->centralBody_,
                                        orbitalStateBasedRotationSettings->isColinearWithVelocity_,
                                        orbitalStateBasedRotationSettings->directionIsOppositeToVector_ ),
                                body ),
                        Eigen::Vector3d::UnitX( ),
                        orbitalStateBasedRotationSettings->getOriginalFrame( ),
                        orbitalStateBasedRotationSettings->getTargetFrame( ),
                        orbitalStateBasedRotationSettings->freeRotationAngleFunction_ );
            }
            break;
        }
        case custom_rotation_model: {
            // Check whether settings for simple rotation model are consistent with its type.
            std::shared_ptr< CustomRotationModelSettings > customRotationSettings =
                    std::dynamic_pointer_cast< CustomRotationModelSettings >( rotationModelSettings );
            if( customRotationSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected custom rotation model settings for " + body );
            }
            else
            {
                std::function< Eigen::Quaterniond( const double ) > customOrientationFunction = [ = ]( const double time ) {
                    return Eigen::Quaterniond( customRotationSettings->customOrientationFunction_( time ) );
                };

                rotationalEphemeris = std::make_shared< CustomRotationalEphemeris >( customOrientationFunction,
                                                                                     customRotationSettings->getOriginalFrame( ),
                                                                                     customRotationSettings->getTargetFrame( ),
                                                                                     customRotationSettings->finiteDifferenceTimeStep_ );
            }
            break;
        }

        default:
            throw std::runtime_error( "Error, did not recognize rotation model settings type " +
                                      std::to_string( rotationModelSettings->getRotationType( ) ) );
    }

    return rotationalEphemeris;
}

}  // namespace simulation_setup

}  // namespace tudat
