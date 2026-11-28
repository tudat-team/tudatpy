/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GETINITIALSTATEOFBODY_H
#define TUDAT_GETINITIALSTATEOFBODY_H

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace propagators
{
//! Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time.
/*!
 * Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time.
 * \param bodiesToIntegrate List of bodies for which to retrieve state.
 * \param centralBodies Origins w.r.t. which to retrieve states of bodiesToIntegrate.
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve states.
 * \param frameManager OBject with which to calculate frame origin translations.
 * \return Initial state vector (with 6 Cartesian elements per body, in order of bodiesToIntegrate vector).
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStatesOfBodiesFromFrameManager(
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SystemOfBodies& bodies,
        const TimeType initialTime,
        const std::shared_ptr< ephemerides::ReferenceFrameManager > frameManager )
{
    // Set initial states of bodies to integrate.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( bodiesToIntegrate.size( ) * 6, 1 );
    std::shared_ptr< ephemerides::Ephemeris > ephemerisOfCurrentBody;

    // Iterate over all bodies.
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        ephemerisOfCurrentBody = bodies.at( bodiesToIntegrate.at( i ) )->getEphemeris( );

        if( !ephemerisOfCurrentBody )
        {
            throw std::runtime_error( "Could not determine initial state for body " + bodiesToIntegrate.at( i ) +
                                      " because it does not have a valid Ephemeris object." );
        }

        // Get body initial state from ephemeris
        systemInitialState.segment( i * 6, 6 ) =
                ephemerisOfCurrentBody->getTemplatedStateFromEphemeris< StateScalarType, TimeType >( initialTime );

        // Correct initial state if integration origin and ephemeris origin are not equal.
        if( centralBodies.at( i ) != ephemerisOfCurrentBody->getReferenceFrameOrigin( ) )
        {
            std::shared_ptr< ephemerides::Ephemeris > correctionEphemeris =
                    frameManager->getEphemeris( ephemerisOfCurrentBody->getReferenceFrameOrigin( ), centralBodies.at( i ) );
            systemInitialState.segment( i * 6, 6 ) -=
                    correctionEphemeris->getTemplatedStateFromEphemeris< StateScalarType, TimeType >( initialTime );
        }
    }
    return systemInitialState;
}

//! Function to get the rotational states states of a set of bodies, at the requested time.
/*!
 * Function to get the rotational states states of a set of bodies, at the requested time.
 * \param bodiesToIntegrate List of bodies for which to retrieve rotational state.
 * \param baseOrientations Reference base frame orientation
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve states.
 * \return Initial rotational state vector (with 7 elements: 4 for quaternion; 3 for angular velocity) per body,
 * in order of bodiesToIntegrate vector).
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialRotationalStatesOfBodies( const std::vector< std::string >& bodiesToIntegrate,
                                                                                        const std::vector< std::string >& baseOrientations,
                                                                                        const simulation_setup::SystemOfBodies& bodies,
                                                                                        const TimeType initialTime )
{
    // Set initial states of bodies to integrate.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( bodiesToIntegrate.size( ) * 7, 1 );
    std::shared_ptr< ephemerides::RotationalEphemeris > rotationModelOfCurrentBody;

    // Iterate over all bodies.
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        rotationModelOfCurrentBody = bodies.at( bodiesToIntegrate.at( i ) )->getRotationalEphemeris( );

        if( !rotationModelOfCurrentBody )
        {
            throw std::runtime_error( "Could not determine initial state for body " + bodiesToIntegrate.at( i ) +
                                      " because it does not have a valid RotationalEphemeris object." );
        }

        // Get body initial state from ephemeris
        systemInitialState.segment( i * 7, 7 ) =
                rotationModelOfCurrentBody->getRotationStateVector( initialTime ).template cast< StateScalarType >( );

        // Correct initial state if integration origin and rotation model origin are not equal.
        if( baseOrientations.at( i ) != rotationModelOfCurrentBody->getBaseFrameOrientation( ) )
        {
            throw std::runtime_error( "Error, cannot get initial rotational state w.r.t. non-base frame" );
        }
    }
    return systemInitialState;
}

//! Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time.
/*!
 * Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time, creates
 * frameManager from input data.
 * \param bodiesToIntegrate List of bodies for which to retrieve state.
 * \param centralBodies Origins w.r.t. which to retrieve states of bodiesToIntegrate.
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve states.
 * \return Initial state vector (with 6 Cartesian elements per body, in order of bodiesToIntegrate vector).
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStatesOfBodies( const std::vector< std::string >& bodiesToIntegrate,
                                                                              const std::vector< std::string >& centralBodies,
                                                                              const simulation_setup::SystemOfBodies& bodies,
                                                                              const TimeType initialTime )
{
    // Create ReferenceFrameManager and call overloaded function.
    return getInitialStatesOfBodiesFromFrameManager< TimeType, StateScalarType >(
            bodiesToIntegrate, centralBodies, bodies, initialTime, simulation_setup::createFrameManager( bodies.getMap( ) ) );
}

//! Function to get the states of single body, w.r.t. some central body, at the requested time.
/*!
 * Function to get the states of  single body, w.r.t. some central body, at the requested time. This function creates
 * frameManager from input data to perform all required conversions.
 * \param bodyToIntegrate Body for which to retrieve state
 * \param centralBody Origin w.r.t. which to retrieve state of bodyToIntegrate.
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve state.
 * \return Initial state vector of bodyToIntegrate
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStateOfBody( const std::string& bodyToIntegrate,
                                                                           const std::string& centralBody,
                                                                           const simulation_setup::SystemOfBodies& bodies,
                                                                           const TimeType initialTime )
{
    return getInitialStatesOfBodies< TimeType, StateScalarType >( { bodyToIntegrate }, { centralBody }, bodies, initialTime );
}

//! Function to get the rotational states state of a body, at the requested time.
/*!
 * Function to get the rotational states state of a body, at the requested time..
 * \param bodyToIntegrate Body for which to retrieve rotational state.
 * \param baseOrientation Reference base frame orientation
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve states.
 * \return Initial rotational state vector (with 7 elements: 4 for quaternion; 3 for angular velocity)
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialRotationalStateOfBody( const std::string& bodyToIntegrate,
                                                                                     const std::string& baseOrientation,
                                                                                     const simulation_setup::SystemOfBodies& bodies,
                                                                                     const TimeType initialTime )
{
    return getInitialRotationalStatesOfBodies< TimeType, StateScalarType >(
            std::vector< std::string >{ bodyToIntegrate }, std::vector< std::string >{ baseOrientation }, bodies, initialTime );
}

//! Function to get the state of single body, w.r.t. some central body, at a set of requested times, concatanated into one vector.
/*!
 * Function to get the states of  single body, w.r.t. some central body, at a set of requested times, concatanated into one vector.
 * This function creates frameManager from input data to perform all required conversions.
 * \param bodyToIntegrate Body for which to retrieve state
 * \param centralBodies Origins w.r.t. which to retrieve state of bodyToIntegrate (per arc).
 * \param bodies List of bodies to use in simulations.
 * \param arcStartTimes List of times at which to retrieve states.
 * \return Initial state vectosr of bodyToIntegrate at requested times.
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialArcWiseStateOfBody( const std::string& bodyToIntegrate,
                                                                                  const std::vector< std::string >& centralBodies,
                                                                                  const simulation_setup::SystemOfBodies& bodies,
                                                                                  const std::vector< TimeType > arcStartTimes )
{
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 6 * arcStartTimes.size( ), 1 );
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        initialStates.block( 6 * i, 0, 6, 1 ) =
                getInitialStateOfBody< double, StateScalarType >( bodyToIntegrate, centralBodies.at( i ), bodies, arcStartTimes.at( i ) );
    }
    return initialStates;
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_GETINITIALSTATEOFBODY_H
