/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <boost/assign/list_of.hpp>

#include <vector>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h"
#include "Tudat/Mathematics/Interpolators/multiLinearInterpolator.h"
#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createAerodynamicCoefficientInterface.h"


namespace tudat
{

namespace simulation_setup
{

//! Function to create a flight conditions object
/*!
 * Function to create a flight conditions object, which is responsible for calculating the various
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
boost::shared_ptr< aerodynamics::FlightConditions > createFlightConditions(
        const boost::shared_ptr< Body > bodyWithFlightConditions,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::function< double( ) > angleOfAttackFunction = boost::function< double( ) >( ),
        const boost::function< double( ) > angleOfSideslipFunction = boost::function< double( ) >( ),
        const boost::function< double( ) > bankAngleFunction = boost::function< double( ) >( ),
        const boost::function< void( const double ) > angleUpdateFunction = boost::function< void( const double ) >( ) );


//! Function to set the angle of attack to trimmed conditions.
/*!
 * Function to set the angle of attack to trimmed conditions. Using this function requires teh aerodynamic coefficient
 * interface to be dependent on the angle of attack.
 * \param flightConditions Flight conditions for body that is to have trimmed conditions.
 */
boost::shared_ptr< aerodynamics::TrimOrientationCalculator > setTrimmedConditions(
        const boost::shared_ptr< aerodynamics::FlightConditions > flightConditions );


//! Function to set the angle of attack to trimmed conditions.
/*!
 * Function to set the angle of attack to trimmed conditions. Using this function requires teh aerodynamic coefficient
 * interface to be dependent on the angle of attack.
 * \param bodyWithFlightConditions Body for which trimmed conditions are to be imposed.
 */
boost::shared_ptr< aerodynamics::TrimOrientationCalculator > setTrimmedConditions(
        const boost::shared_ptr< Body > bodyWithFlightConditions );


//! Function that must be called to link the AerodynamicGuidance object to the simulation
/*!
 * Function that must be called to link the AerodynamicGuidance object to the simulation
 * \param aerodynamicGuidance Object computing the current aerodynamic angles.
 * \param angleCalculator Object that handles all aerodynamic angles in the numerical propagation
 */
void setGuidanceAnglesFunctions(
        const boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator );

//! Function that must be called to link the AerodynamicGuidance object to the simulation
/*!
 * Function that must be called to link the AerodynamicGuidance object to the simulation
 * \param aerodynamicGuidance Object computing the current aerodynamic angles.
 * \param bodyWithAngles Body for which the orientation is to be controlled.
 */
void setGuidanceAnglesFunctions(
        const boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance,
        const boost::shared_ptr< simulation_setup::Body > bodyWithAngles );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEACCELERATIONMODELS_H
