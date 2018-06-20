/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEBODIES_H
#define TUDAT_CREATEBODIES_H

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodyShapeModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRadiationPressureInterface.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"

namespace tudat
{

namespace simulation_setup
{

//! Struct holding settings for a body to be created.
/*!
 *  Struct holding settings for a body to be created. From the settings, a Body object is
 *  created by the createBodies function. Default values can be generated from the function in
 *  defaultBodies.h.
 */
struct BodySettings
{
    //! Constant mass.
    double constantMass = TUDAT_NAN;

    //! Settings for the atmosphere model that the body is to contain.
    std::shared_ptr< AtmosphereSettings > atmosphereSettings;

    //! Settings for the ephemeris model that the body is to contain.
    std::shared_ptr< EphemerisSettings > ephemerisSettings;

    //! Settings for the gravity field model that the body is to contain.
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings;

    //! Settings for the rotation model that the body is to contain.
    std::shared_ptr< RotationModelSettings > rotationModelSettings;

    //! Settings for the shape model that the body is to contain.
    std::shared_ptr< BodyShapeSettings > shapeModelSettings;

    //! Settings for the radiations pressure interfaces that the body is to contain (source body as key).
    std::map< std::string,
    std::shared_ptr< RadiationPressureInterfaceSettings > > radiationPressureSettings;

    //! Settings for the aerodynamic coefficients that the body is to contain.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;

    //! Settings for variations of the gravity field of the body.
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;

    std::vector< std::shared_ptr< GroundStationSettings > > groundStationSettings;

};

//! Function that determines the order in which bodies are to be created
/*!
 * Function that determines the order in which bodies are to be created, to ensure that any dependency between body models are
 * correctly handled. Currently, no dependencies exist that force any particular creation order.
 * \param bodySettings Map of body settings (with map key the body name)
 * \return List of pairs: name and body settings of that body
 */
std::vector< std::pair< std::string, std::shared_ptr< BodySettings > > > determineBodyCreationOrder(
        const std::map< std::string, std::shared_ptr< BodySettings > >& bodySettings );

//! Function to create a map of bodies objects.
/*!
 *  Function to create a msap of body objects based on model-specific settings for the bodies,
 *  containing settings for each relevant environment model.
 *  \param bodySettings List of settings for the bodies that are to be created, defined as a map of
 *  pointers to an object of class BodySettings
 *  \return List of bodies created according to settings in bodySettings.
 */
NamedBodyMap createBodies(
        const std::map< std::string, std::shared_ptr< BodySettings > >& bodySettings );


//! Function to define the global origin and orientation of the reference frame
/*!
 * Function to define the global origin and orientation of the reference frame that is to be used in
 * the simulations.  This function checks the origin and orientation of the Ephemeris and
 * RotationalEphemeris, and checks whether their origin/orientation is the same as that
 * globalFrameOrigin and globalFrameOrientation provided as input.  In particular, this function
 * sets the ephemerisFrameToBaseFrameFunction_ anf ephemerisFrameToBaseFrameLongFunction_ variables
 * of the Body objects, which provide a time-dependent translation of the global origin to the
 * body's ephemeris origin. In case of an inconsistency in the current and requried frames, this
 * function throws an error.
 * \param bodyMap List of body objects that constitute the environment.
 * \param globalFrameOrigin Global reference frame origin.
 * \param globalFrameOrientation Global referencef frame orientation.
 */
template< typename StateScalarType = double, typename TimeType = double >
void setGlobalFrameBodyEphemerides( const NamedBodyMap& bodyMap,
                                    const std::string& globalFrameOrigin,
                                    const std::string& globalFrameOrientation )
{
    using namespace tudat::simulation_setup;
    std::string ephemerisFrameOrigin;
    std::string ephemerisFrameOrientation;
    std::string rotationModelFrame;

    std::vector< std::string > globalFrameOriginChain;

    // Get chain of ephemeris frame origins of global frame origin (if it is not SSB
    if( globalFrameOrigin != "SSB" )
    {
        if( bodyMap.count( globalFrameOrigin ) == 0 )
        {
            throw std::runtime_error(
                        "Error, body non-barycentric global frame origin selected, but this body " + globalFrameOrigin +
                        " is not found." );
        }
        else
        {
            std::string currentOrigin = globalFrameOrigin;
            while( currentOrigin != "SSB" )
            {
                std::shared_ptr< ephemerides::Ephemeris > currentEphemeris =
                        bodyMap.at( currentOrigin )->getEphemeris( );
                if( currentEphemeris == nullptr )
                {
                    throw std::runtime_error(
                                "Error, body non-barycentric global frame origin selected, but body " + currentOrigin +
                                " in chain has no ephemeris." );
                }
                else
                {
                    currentOrigin = currentEphemeris->getReferenceFrameOrigin( );
                    ephemerisFrameOrientation = currentEphemeris->getReferenceFrameOrientation( );
                    if( ephemerisFrameOrientation != globalFrameOrientation )
                    {
                        throw std::runtime_error(
                                    "Error, ephemeris orientation of body " + currentOrigin
                                    + " is not the same as global orientation " + ephemerisFrameOrientation
                                    + ", " + globalFrameOrientation );
                    }
                }

                if( std::find( globalFrameOriginChain.begin( ), globalFrameOriginChain.end( ), currentOrigin ) !=
                        globalFrameOriginChain.end( ) )
                {
                    throw std::runtime_error(
                                "Error, body non-barycentric global frame origin selected, but body " + currentOrigin +
                                " already found in origin chain." );
                }
                else
                {
                    globalFrameOriginChain.push_back( currentOrigin );
                }
            }
        }
    }

    // Iterate over all bodies
    for( NamedBodyMap::const_iterator bodyIterator = bodyMap.begin( );
         bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        // Check id body contains an ephemeris
        if( bodyIterator->second->getEphemeris( ) != nullptr )
        {
            // Retrieve ephemeris origin
            ephemerisFrameOrigin = bodyIterator->second->getEphemeris( )->getReferenceFrameOrigin( );

            // Check if ephemeris origin differs from global origin.
            if( ephemerisFrameOrigin != globalFrameOrigin )
            {
                // Make correction to SSB if it is global frame origin
                if( globalFrameOrigin == "SSB" )
                {
                    // Check if correction can be made
                    if( bodyMap.count( ephemerisFrameOrigin ) == 0 )
                    {
                        throw std::runtime_error(
                                    "Error, body " + bodyIterator->first + " has ephemeris in frame " +
                                    ephemerisFrameOrigin + ", but no conversion to frame " + globalFrameOrigin +
                                    " can be made" );
                    }
                    else
                    {
                        std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction =
                                std::bind( &Body::getStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                                             bodyMap.at( ephemerisFrameOrigin ), std::placeholders::_1 );
                        std::shared_ptr< BaseStateInterface > baseStateInterface =
                                std::make_shared< BaseStateInterfaceImplementation< TimeType, StateScalarType > >(
                                    ephemerisFrameOrigin, stateFunction );
                        bodyIterator->second->setEphemerisFrameToBaseFrame( baseStateInterface );
                    }
                }
                // Make correction to global frame origin (if not SSB)
                else
                {
                    // Set barycentric state function of global frame origin
                    if( globalFrameOrigin == bodyIterator->first )
                    {
                        std::shared_ptr< ephemerides::ReferenceFrameManager > frameManager =
                                propagators::createFrameManager( bodyMap );
                        frameManager->getEphemeris( globalFrameOrigin, "SSB" );

                        std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction =
                                std::bind( &ephemerides::Ephemeris::getTemplatedStateFromEphemeris< StateScalarType, TimeType >,
                                             frameManager->getEphemeris( globalFrameOrigin, "SSB" ), std::placeholders::_1 );

                        std::shared_ptr< BaseStateInterface > baseStateInterface =
                                std::make_shared< BaseStateInterfaceImplementation< TimeType, StateScalarType > >(
                                    globalFrameOrigin, stateFunction, true );
                        bodyIterator->second->setEphemerisFrameToBaseFrame( baseStateInterface );

                    }
                    // Set correction function if ephemeris origin is SSB
                    else if( ephemerisFrameOrigin == "SSB" )
                    {

                        std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction =
                               std::bind( &Body::getGlobalFrameOriginBarycentricStateFromEphemeris< StateScalarType, TimeType >,
                                             bodyMap.at( globalFrameOrigin ), std::placeholders::_1 );
                        std::shared_ptr< BaseStateInterface > baseStateInterface =
                                std::make_shared< BaseStateInterfaceImplementation< TimeType, StateScalarType > >(
                                    globalFrameOrigin, stateFunction, true );
                        bodyIterator->second->setEphemerisFrameToBaseFrame( baseStateInterface );
                    }
                    else
                    {
                        // Check if correction can be made
                        if( bodyMap.count( ephemerisFrameOrigin ) == 0 )
                        {
                            throw std::runtime_error(
                                        "Error, body " + bodyIterator->first + " has ephemeris in frame " +
                                        ephemerisFrameOrigin + ", but no conversion to frame " + globalFrameOrigin +
                                        " can be made" );
                        }
                        else
                        {
                            // Set correction function from ephemeris origin to global frame origin
                            std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction =
                                    std::bind( &Body::getStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                                                 bodyMap.at( ephemerisFrameOrigin ), std::placeholders::_1 );
                            std::shared_ptr< BaseStateInterface > baseStateInterface =
                                    std::make_shared< BaseStateInterfaceImplementation< TimeType, StateScalarType > >(
                                        ephemerisFrameOrigin, stateFunction, false );
                            bodyIterator->second->setEphemerisFrameToBaseFrame( baseStateInterface );
                        }
                    }
                }
            }

            // Retrieve ephemeris orientation
            ephemerisFrameOrientation = bodyIterator->second->getEphemeris( )->getReferenceFrameOrientation( );
            // If two are not equal, throw error.
            if( ephemerisFrameOrientation != globalFrameOrientation )
            {
                throw std::runtime_error(
                            "Error, ephemeris orientation of body " + bodyIterator->first
                            + " is not the same as global orientation " + ephemerisFrameOrientation
                            + ", " + globalFrameOrientation );
            }


        }

        // Set global frame origin identifiers
        if( globalFrameOrigin == bodyIterator->first )
        {
            bodyIterator->second->setIsBodyGlobalFrameOrigin( 1 );
        }
        else
        {
            bodyIterator->second->setIsBodyGlobalFrameOrigin( 0 );
        }

        // Check if body has rotational ephemeris.
        if( bodyIterator->second->getRotationalEphemeris( ) != nullptr )
        {
            // Check if rotational ephemeris base frame orienatation is equal to to global orientation.
            rotationModelFrame = bodyIterator->second->getRotationalEphemeris( )->getBaseFrameOrientation( );

            // Throw error if two frames are not equal.
            if( rotationModelFrame != globalFrameOrientation )
            {
                throw std::runtime_error(
                            "Error, rotation base orientation of body " + bodyIterator->first +
                            " is not the same as global orientation " + rotationModelFrame + ", " +
                            globalFrameOrientation );
            }
        }
    }

    // Set body state-dependent environment variables
    for( NamedBodyMap::const_iterator bodyIterator = bodyMap.begin( );
         bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        bodyIterator->second->updateConstantEphemerisDependentMemberQuantities( );
    }

}

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEBODIES_H
