/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FULLPROPAGATIONLAMBERTTARGETER_H
#define TUDAT_FULLPROPAGATIONLAMBERTTARGETER_H

#include <string>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterIzzo.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"

namespace tudat
{

namespace propagators
{

//! Function to directly setup a body map corresponding to the assumptions of the Lambert targeter.
/*!
 * Function to drectly setup Lambert targeter map. The body map only contains the central body and the body to be propagated.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \return
 */
simulation_setup::NamedBodyMap setupBodyMapLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate);



//! Function to directly setup an acceleration map for the Lambert targeter.
/*!
 * Function to directly setup an acceleration map for the Lambert targeter. Only the central body exert a point-mass gravity acceleration
 * upon the body to be propagated.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param bodyMap Body map for the Lambert targeter.
 * \return
 */
basic_astrodynamics::AccelerationMap setupAccelerationMapLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const simulation_setup::NamedBodyMap& bodyMap );


//! Function to determine the cartesian state at a given time for a keplerian orbit, based on the initial state.
/*!
 * Function to determine the cartesian state at a given time for a keplerian orbit, based on the initial state.
 * \param initialState Initial cartesian state on this orbit (x-position coordinate [m], y-position coordinate [m], z-position coordinate [m],
 * x-velocity coordinate [m/s], y-velocity coordinate [m/s], z-velocity coordinate [m/s]).
 * \param finalPropagationTime Final time at which the cartesian state must be given [s].
 * \param gravitationalParameter Gravitation parameter defining the keplerian orbit [m^3 s^-2].
 * \return
 */
Eigen::Vector6d propagateLambertTargeterSolution(
        const Eigen::Vector6d& initialState,
        const double finalPropagationTime,
        const double gravitationalParameter);



//! Function to propagate the full dynamics problem and the Lambert targeter solution.
/*!
 * Function to propagate the full dynamics problem and the Lambert targeter solution. The function computes the cartesian states
 * obtained with the Lambert targeter and after propagation of the full dynamics problem as a function of time.
 * \param cartesianPositionAtDeparture Cartesian position of the body to be propagated at departure [m].
 * \param cartesianPositionAtArrival Cartesian position of the body to be propagated at arrival [m].
 * \param timeOfFlight Time of flight [s].
 * \param bodyMap Body map.
 * \param accelerationModelMap Acceleration map.
 * \param bodiesToPropagate Vector with the name of the bodies to be propagated.
 * \param centralBodies Vector with the name of the central bodies.
 * \param integratorSettings Integrator settings for the propagation.
 * \param lambertTargeterResult Map of the cartesian state obtained with the Lambert targeter as a function of time (modified within
 * the function).
 * \param fullProblemResult Map of the cartesian state obtained after propagation of the full dynamics problem as a function of time (modified
 *  within the function).
 */
void propagateLambertTargeterAndFullProblem(
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        const std::map< double, Eigen::Vector6d >& fullProblemResult);





//! Function to compute the difference in cartesian state between Lambert targeter solution and full dynamics problem, both at departure
//! and at arrival.
/*!
 * Function to compute the difference in cartesian state between Lambert targeter solution and full dynamics problem, both at departure
 * and at arrival. The function returns a pair of vectors, the first one being the difference in state at departure and the second one the
 * difference in state at arrival.
 * \param cartesianPositionAtDeparture Cartesian position of the body to be propagated at departure [m].
 * \param cartesianPositionAtArrival Cartesian position of the body to be propagated at arrival [m].
 * \param timeOfFlight Time of flight [s].
 * \param bodyMap Body map.
 * \param accelerationModelMap Acceleration map.
 * \param bodiesToPropagate Vector with the names of the bodies to be propagated.
 * \param centralBody Vector with the names of the central bodies.
 * \param integratorSettings Integrator settings for the propagation.
 * \return
 */
std::pair< Eigen::Vector6d, Eigen::Vector6d > getDifferenceFullPropagationWrtLambertTargeterAtDepartureAndArrival(
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings);


} // namespace propagators

} // namespace tudat

#endif // TUDAT_FULLPROPAGATIONLAMBERTTARGETER_H
