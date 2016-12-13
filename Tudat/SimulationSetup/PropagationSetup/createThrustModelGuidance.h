/*    Copyright (c) 2010-2016, Delft University of Technology
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
        const std::vector< boost::function< Eigen::Vector3d( )> >& thrustDirections,
        const std::vector< boost::function< double( )> >& thrustMagnitudes );

//! Function to create a function that returns the thrust direction in the body-fixed frame.
/*!
 * Function to create a function that returns the thrust direction in the body-fixed frame.
 * \param thrustMagnitudeSettings Settings for the thrust magnitude
 * \param bodyMap List of body objects that comprises the environment
 * \param bodyName Name of body for which thrust is to be created.
 * \return Function that returns the thrust direction in the body-fixed frame.
 */
boost::function< Eigen::Vector3d( ) > getBodyFixedThrustDirection(
        const boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings,
        const NamedBodyMap& bodyMap,
        const std::string bodyName );

//! Function to create a list of functions that (compute and) return independent variables for thrust
/*!
 * Function to create a list of functions that (compute and) return independent variables for thrust and/or specific impulse.
 * This parameterization is used in the thrust mangitude type is thrust_magnitude_from_dependent_variables. This function
 * retrieves all input functions from the environment and a list of user-defined functions.
 * \param bodyWithGuidance Name of body for which the propulsion settings are to be retrieved.
 * \param independentVariables List of variables for which function returning them are to be created. Note that the number
 * of guidance_input_dependent_thrust entries must be equal to the size of guidanceInputFunctions. No entries of type
 * maximum_thrust_multiplier are allowed.
 * \param guidanceInputFunctions Functions returning user-defined variables on which the thrust/specific impulse depends
 * \return List of functions that (compute and) return independent variables for thrust
 */
std::vector< boost::function< double( ) > > getPropulsionInputVariables(
        const boost::shared_ptr< Body > bodyWithGuidance,
        const std::vector< propulsion::ThrustDependentVariables > independentVariables,
        const std::vector< boost::function< double( ) > > guidanceInputFunctions );

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
boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > createThrustMagnitudeWrapper(
        const boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings,
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
        const boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
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
        const boost::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
        const boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance,
        const double currentTime = TUDAT_NAN );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATETHRUSTMODELGUIDANCE_H
