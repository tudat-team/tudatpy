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

//! Function to setup a body map corresponding to the assumptions of the Lambert targeter,
//! retrieving positions of departure and arrival bodies from ephemerides.
/*!
 * Function to setup Lambert targeter map. The body map only contains the central body and the body to be propagated.
 * The positions of the departure and arrival bodies are directly retrived from ephemerides.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param departureAndArrivalBodies Vector containing the names of the departure and arrival bodies.
 * \return Body map for the Lambert targeter
 */
simulation_setup::NamedBodyMap setupBodyMapFromEphemeridesForLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& departureAndArrivalBodies );



//! Function to setup a body map corresponding to the assumptions of the Lambert targeter,
//! the positions of departure and arrival bodies being provided as inputs.
/*!
 * Function to setup Lambert targeter map. The body map only contains the central body and the body to be propagated.
 * The positions of the departure and arrival bodies are defined by the user and provided as inputs.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param departureAndArrivalBodies Vector containing the names of the departure and arrival bodies.
 * \param cartesianPositionAtDeparture Vector containing the position coordinates of the departure body [m].
 * \param cartesianPositionAtArrival Vector containing the position coordinates of the arrival body [m].
 * \return Body map for the Lambert targeter.
 */
simulation_setup::NamedBodyMap setupBodyMapFromUserDefinedStatesForLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& departureAndArrivalBodies,
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival );


//! Function to directly setup an acceleration map for the Lambert targeter.
/*!
 * Function to directly setup an acceleration map for the Lambert targeter. Only the central body exerts a point-mass gravity acceleration
 * on the body to be propagated.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param bodyMap Body map for the Lambert targeter.
 * \return Acceleration map for the Lambert targeter.
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
 * \param finalPropagationTime Final time at which the cartesian state has to be calculated [s].
 * \param gravitationalParameter Gravitation parameter defining the keplerian orbit [m^3 s^-2].
 * \return Vector containing the cartesian state at a given time for a keplerian orbit.
 */
Eigen::Vector6d computeCartesianStateFromKeplerianOrbit(
        const Eigen::Vector6d& initialState,
        const double finalPropagationTime,
        const double gravitationalParameter);


//! Function to compute the cartesian state at half of the time of flight for a Lambert targeter.
/*!
 * Function to compute the cartesian state at half of the time of flight for a Lambert targeter.
 * \param cartesianStateAtDeparture Vector containing the cartesian state at departure (x-position coordinate [m], y-position coordinate [m],
 * z-position coordinate [m], x-velocity coordinate [m/s], y-velocity coordinate [m/s], z-velocity coordinate [m/s]).
 * \param gravitationalParameterCentralBody Gravitational parameter of the central body of the Lambert targeter [m^3 s^-2].
 * \param timeOfFlight Time of flight of the Lambert Targeter [s].
 * \return Vector containing the cartesian state at half of the time of flight.
 */
Eigen::Vector6d computeCartesianStateHalfTimeOfFlightLambertTargeter(
        const Eigen::Vector6d& cartesianStateAtDeparture,
        const double gravitationalParameterCentralBody,
        const double timeOfFlight);



//! Function to propagate the full dynamics problem and the Lambert targeter solution.
/*!
 * Function to propagate the full dynamics problem and the Lambert targeter solution. The function computes the cartesian state as a function of time in two
 * different ways: from the Lambert targeter and from the propagation of the full dynamics problem.
 * \param cartesianPositionAtDeparture Cartesian position of the body to be propagated at departure [m].
 * \param cartesianPositionAtArrival Cartesian position of the body to be propagated at arrival [m].
 * \param timeOfFlight Time of flight [s].
 * \param initialTime Initial time of the propagation [s].
 * \param bodyMap Body map.
 * \param accelerationModelMap Acceleration map.
 * \param bodyToPropagate Name of the body to be propagated.
 * \param centralBody Name of the central body for the Lambert targeter.
 * \param integratorSettings Integrator settings for the propagation.
 * \param lambertTargeterResult Map of the cartesian state obtained with the Lambert targeter as a function of time (modified within
 * the function).
 * \param fullProblemResult Map of the cartesian state obtained after propagation of the full dynamics problem as a function of time (modified
 *  within the function).
 * \param arrivalAndDepartureInitialisationFromEphemerides Boolean denoting whether the positions of the departure and arrival bodies are retrieved from
 * ephemerides (true) or defined by the input vectors cartesianPositionAtDeparture and cartesianPositionAtArrival (false). The default value of this boolean is
 * false.
 * \param terminationSphereOfInfluence Boolean denoting whether the propagation stops at the position of the arrival body (false) or at the sphere of influence
 * of the arrival body (true). The default value of this boolean is false.
 * \param departureBodyGravitationalParameterParameter Gravitational parameter of the departure body [m^3 s^-2]. If not provided as input, it is retrieved from
 * the body map.
 * \param arrivalBodyGravitationalParameter Gravitational parameter of the arrival body [m^3 s^-2]. If not provided as input, it is retrieved from
 * the body map.
 * \param centralBodyGravitationalParameter Gravitational parameter of the central body [m^3 s^-2]. If not provided as input, it is retrieved from
 * the body map.
 */
void propagateLambertTargeterAndFullProblem(
        Eigen::Vector3d cartesianPositionAtDeparture,
        Eigen::Vector3d cartesianPositionAtArrival,
        const double timeOfFlight,
        const double initialTime,
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        const std::vector<std::string>& departureAndArrivalBodies,
        const bool arrivalAndDepartureInitialisationFromEphemerides = false,
        const bool terminationSphereOfInfluence = false,
        const double departureBodyGravitationalParameter = TUDAT_NAN,
        const double arrivalBodyGravitationalParameter = TUDAT_NAN,
        const double centralBodyGravitationalParameter = TUDAT_NAN);





//! Function to compute the difference in cartesian state between Lambert targeter solution and full dynamics problem, both at departure
//! and at arrival.
/*!
 * Function to compute the difference in cartesian state between Lambert targeter solution and full dynamics problem, both at departure
 * and at arrival. The function returns a pair of vectors, the first one being the difference in state at departure and the second one the
 * difference in state at arrival.
 * \param cartesianPositionAtDeparture Cartesian position of the body to be propagated at departure [m].
 * \param cartesianPositionAtArrival Cartesian position of the body to be propagated at arrival [m].
 * \param timeOfFlight Time of flight [s].
 * \param initialTime Initial time of the propagation [s].
 * \param bodyMap Body map.
 * \param accelerationModelMap Acceleration map.
 * \param bodiesToPropagate Vector with the names of the bodies to be propagated.
 * \param centralBody Vector with the names of the central bodies.
 * \param integratorSettings Integrator settings for the propagation.
 * \param arrivalAndDepartureInitialisationFromEphemerides Boolean denoting whether the positions of the departure and arrival bodies are retrieved from
 * ephemerides (true) or defined by the input vectors cartesianPositionAtDeparture and cartesianPositionAtArrival (false). The default value of this boolean is
 * false.
 * \param terminationSphereOfInfluence Boolean denoting whether the propagation stops at the position of the arrival body (false) or at the sphere of influence
 * of the arrival body (true). The default value of this boolean is false.
 * \return Pair of vectors containing the difference between the Lambert targeter and the full problem cartesian states
 * (at departure and at arrival respectively).
 */
std::pair< Eigen::Vector6d, Eigen::Vector6d > getDifferenceFullPropagationWrtLambertTargeterAtDepartureAndArrival(
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double timeOfFlight,
        const double initialTime,
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::vector< std::string >& departureAndArrivalBodies,
        const bool arrivalAndDepartureInitialisationFromEphemerides = false,
        const bool terminationSphereOfInfluence = false);



} // namespace propagators

} // namespace tudat

#endif // TUDAT_FULLPROPAGATIONLAMBERTTARGETER_H
