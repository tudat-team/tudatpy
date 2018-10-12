/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATETHRUSTMODELGUIDANCE_H
#define TUDAT_CREATETHRUSTMODELGUIDANCE_H

#include "Tudat/Astrodynamics/SystemModels/engineModel.h"
#include "Tudat/Astrodynamics/Propulsion/thrustGuidance.h"
#include "Tudat/Astrodynamics/Propagators/environmentUpdateTypes.h"
#include "Tudat/Astrodynamics/Propulsion/thrustMagnitudeWrapper.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"
#include "Tudat/SimulationSetup/PropagationSetup/thrustSettings.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/SystemModels/engineModel.h"

namespace tudat
{

namespace simulation_setup
{


//! Function to retrieve the effective thrust direction from a set of thrust sources.
/*!
 * Function to retrieve the effective thrust direction from a set of thrust sources.
 * \param thrustDirections List of functions returning thrust directions.
 * \param thrustMagnitudes List of functions returning thrust magnitude.
 * \return Effective thrust direction.
 */
Eigen::Vector3d getCombinedThrustDirection(
        const std::vector< std::function< Eigen::Vector3d( )> >& thrustDirections,
        const std::vector< std::function< double( )> >& thrustMagnitudes );

//! Function to create a function that returns the thrust direction in the body-fixed frame.
/*!
 * Function to create a function that returns the thrust direction in the body-fixed frame.
 * \param thrustMagnitudeSettings Settings for the thrust magnitude
 * \param bodyMap List of body objects that comprises the environment
 * \param bodyName Name of body for which thrust is to be created.
 * \return Function that returns the thrust direction in the body-fixed frame.
 */
std::function< Eigen::Vector3d( ) > getBodyFixedThrustDirection(
        const std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string bodyName );

//! Function to create a wrapper object that computes the thrust magnitude
/*!
 * Function to create a wrapper object that computes the thrust magnitude
 * \param thrustMagnitudeSettings Settings for the thrust magnitude
 * \param bodyMap List of body objects that comprises the environment
 * \param nameOfBodyWithGuidance Name of body for which thrust is to be created.
 * \param magnitudeUpdateSettings Environment update settings that are required to compute the thrust direction (updated
 * by function as needed).
 * \return Object used during propagation to compute the thrust direction
 */
std::shared_ptr< propulsion::ThrustMagnitudeWrapper > createThrustMagnitudeWrapper(
        const std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings );

//! Function to update the thrust magnitude and direction to current time.
/*!
 * Function to update the thrust magnitude and direction to current time.
 * \param thrustMagnitudeWrapper Object used during propagation to compute the thrust magnitude
 * \param thrustDirectionGuidance Object used during propagation to compute the body-fixed thrust direction
 * \param currentTime Time to which objects are to be updated.
 */
void updateThrustMagnitudeAndDirection(
        const std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
        const double currentTime );

//! Function to reset the current time variable of the thrust magnitude and direction wrappers
/*!
* Function to reset the current time variable of the thrust magnitude and direction wrappers. This function does not
* update the actual thrust direction and guidance; it is typically used to reset the current time to NaN, thereby signalling
* the need to recompute the magnitude/direction upon next call to update functions
* \param thrustMagnitudeWrapper Object used during propagation to compute the thrust magnitude
* \param thrustDirectionGuidance Object used during propagation to compute the body-fixed thrust direction
* \param currentTime New current time variable that is to be set.
*/
void resetThrustMagnitudeAndDirectionTime(
        const std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
        const double currentTime = TUDAT_NAN );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATETHRUSTMODELGUIDANCE_H
