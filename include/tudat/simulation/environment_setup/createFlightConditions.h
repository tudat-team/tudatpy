/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEFLIGHTCONDITIONS_H
#define TUDAT_CREATEFLIGHTCONDITIONS_H

#include <boost/multi_array.hpp>

#include <vector>

#include "tudat/astro/aerodynamics/aerodynamicGuidance.h"
#include "tudat/math/interpolators/multiLinearInterpolator.h"
#include "tudat/astro/aerodynamics/flightConditions.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createAerodynamicCoefficientInterface.h"


namespace tudat
{

namespace simulation_setup
{

//! Function to create an atmospheric flight conditions object
/*!
 * Function to create an atmospheric flight conditions object, which is responsible for calculating the various
 * dependent variables required for calculation of the aerodynamic acceleration
 * \param bodyWithFlightConditions Body for which flight conditions are to be created.
 * \param centralBody Body in  the atmosphere of which bodyWithFlightConditions is flying
 * \param nameOfBodyUndergoingAcceleration Name of body undergoing acceleration.
 * \param nameOfBodyExertingAcceleration Name of body with the atmosphere causing acceleration.
 * \param angleOfAttackFunction Function returning the current angle of attack (default 0).
 * \param angleOfSideslipFunction Function returning the current angle of sideslip (default 0).
 * \param bankAngleFunction Function returning the current bank angle (default 0).
 * \param angleUpdateFunction Function to update the aerodynamic angles to the current time (default none).
 * \return Flight conditions object for given bodies and settings.
 */
std::shared_ptr< aerodynamics::AtmosphericFlightConditions >  createAtmosphericFlightConditions(
        const std::shared_ptr< Body > bodyWithFlightConditions,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::function< double( ) > angleOfAttackFunction = std::function< double( ) >( ),
        const std::function< double( ) > angleOfSideslipFunction = std::function< double( ) >( ),
        const std::function< double( ) > bankAngleFunction = std::function< double( ) >( ),
        const std::function< void( const double ) > angleUpdateFunction = std::function< void( const double ) >( ) );

//! Function to create a flight conditions object
/*!
 * Function to create a flight conditions object, which is responsible for calculating various
 * dependent variables (altitude, latitude, etc ) for non-atmospheric flight
 * \param bodyWithFlightConditions Body for which flight conditions are to be created.
 * \param centralBody Body in  the atmosphere of which bodyWithFlightConditions is flying
 * \param nameOfBodyUndergoingAcceleration Name of body undergoing acceleration.
 * \param nameOfBodyExertingAcceleration Name of body with the atmosphere causing acceleration.
 * \return Flight conditions object for given bodies and settings.
 */
std::shared_ptr< aerodynamics::FlightConditions >  createFlightConditions(
        const std::shared_ptr< Body > bodyWithFlightConditions,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration );

//! Function to set the angle of attack to trimmed conditions.
/*!
 * Function to set the angle of attack to trimmed conditions. Using this function requires the aerodynamic coefficient
 * interface to be dependent on the angle of attack.
 * \param flightConditions Flight conditions for body that is to have trimmed conditions.
 */
std::shared_ptr< aerodynamics::TrimOrientationCalculator > setTrimmedConditions(
        const std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions );


//! Function to set the angle of attack to trimmed conditions.
/*!
 * Function to set the angle of attack to trimmed conditions. Using this function requires the aerodynamic coefficient
 * interface to be dependent on the angle of attack.
 * \param bodyWithFlightConditions Body for which trimmed conditions are to be imposed.
 */
std::shared_ptr< aerodynamics::TrimOrientationCalculator > setTrimmedConditions(
        const std::shared_ptr< Body > bodyWithFlightConditions );


//! Function that must be called to link the AerodynamicGuidance object to the simulation
/*!
 * Function that must be called to link the AerodynamicGuidance object to the simulation
 * \param aerodynamicGuidance Object computing the current aerodynamic angles.
 * \param angleCalculator Object that handles all aerodynamic angles in the numerical propagation
 */
void setGuidanceAnglesFunctions(
        const std::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const std::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator,
        const bool silenceWarnings = false );

//! Function that must be called to link the AerodynamicGuidance object to the simulation
/*!
 * Function that must be called to link the AerodynamicGuidance object to the simulation
 * \param aerodynamicGuidance Object computing the current aerodynamic angles.
 * \param bodyWithAngles Body for which the orientation is to be controlled.
 */
void setGuidanceAnglesFunctions(
        const std::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const std::shared_ptr< simulation_setup::Body > bodyWithAngles,
        const bool silenceWarnings = false );

void setAerodynamicOrientationFunctions(
        const std::shared_ptr< simulation_setup::Body > body,
        const std::function< double( ) > angleOfAttackFunction = std::function< double( ) >( ),
        const std::function< double( ) > angleOfSideslipFunction = std::function< double( ) >( ),
        const std::function< double( ) > bankAngleFunction =  std::function< double( ) >( ),
        const std::function< void( const double ) > angleUpdateFunction = std::function< void( const double ) >( ) );

void setConstantAerodynamicOrientation(
        const std::shared_ptr< simulation_setup::Body > body,
        const double angleOfAttack,
        const double sideslipAngle,
        const double bankAngle,
        const bool silenceWarnings = false );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEACCELERATIONMODELS_H
