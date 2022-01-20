/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;


#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/astro/mission_segments/lambertTargeter.h"
#include "tudat/astro/mission_segments/lambertTargeterIzzo.h"
#include "tudat/astro/mission_segments/lambertRoutines.h"

namespace tudat
{

namespace propagators
{

//! Function to setup a system of bodies corresponding to the assumptions of the Lambert targeter,
//! using default ephemerides for the central, departure and arrival bodies.
/*!
 * Function to setup Lambert targeter map. The system of bodies contains the central body, the body to be propagated and the departure and
 * arrival bodies. The positions of the departure and arrival bodies are directly retrived from ephemerides.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param departureAndArrivalBodies Vector containing the names of the departure and arrival bodies.
 * \return Body map for the Lambert targeter
 */
simulation_setup::SystemOfBodies setupBodyMapFromEphemeridesForLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const  std::pair< std::string, std::string >& departureAndArrivalBodies );


simulation_setup::SystemOfBodies setupBodyMapFromUserDefinedEphemeridesForLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const  std::pair< std::string, std::string >& departureAndArrivalBodies,
        const std::vector< ephemerides::EphemerisPointer >& ephemerisVectorDepartureAndArrivalBodies);

//! Function to setup a system of bodies corresponding to the assumptions of the Lambert targeter,
//! using default ephemerides for the central body only, while the positions of departure and arrival bodies are provided as inputs.
/*!
 * Function to setup a Lambert targeter map. The system of bodies contains the central, departure and arrival bodies and the body to be propagated.
 * The positions of the departure and arrival bodies are defined by the user and provided as inputs.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param departureAndArrivalBodies Vector containing the names of the departure and arrival bodies.
 * \param cartesianPositionAtDeparture Vector containing the position coordinates of the departure body [m].
 * \param cartesianPositionAtArrival Vector containing the position coordinates of the arrival body [m].
 * \return Body map for the Lambert targeter.
 */
simulation_setup::SystemOfBodies setupBodyMapFromUserDefinedStatesForLambertTargeter(
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
 * \param bodies Body map for the Lambert targeter.
 * \return Acceleration map for the Lambert targeter.
 */
basic_astrodynamics::AccelerationMap setupAccelerationMapLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const simulation_setup::SystemOfBodies& bodies );


std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
std::shared_ptr< propagators::PropagationTerminationSettings > > getLambertTargeterTerminationSettings(
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::SystemOfBodies& bodyMap,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::pair< std::string, std::string >& departureAndArrivalBodies,
        const bool setSphereOfInfluenceTermination = true,
        const std::pair< double, double > distanceTerminationAsSoiFraction = std::make_pair( 1.0, 1.0 ),
        const double timeTerminationInSynodicPeriods = 2.0 );

void propagateLambertTargeterAndFullProblem(
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::SystemOfBodies& bodyMap,
        const std::string& centralBody,
        const std::pair< std::string, std::string >& departureAndArrivalBodies ,
        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult );

//! Function to propagate the full dynamics problem and the Lambert targeter solution.
/*!
 * Function to propagate the full dynamics problem and the Lambert targeter solution. The function computes the cartesian state as a function of time in two
 * different ways: from the Lambert targeter and from the propagation of the full dynamics problem. The propagator settings for the full problem
 * propagation are directly provided as inputs.
 * \param timeOfFlight Time of flight [s].
 * \param initialTime Initial time of the propagation [s].
 * \param bodies Body map.
 * \param centralBody Name of the central body for the Lambert targeter.
 * \param propagatorSettings Propagator settings for the full problem.
 * \param integratorSettings Integrator settings for the propagation.
 * \param lambertTargeterResult Map of the cartesian state obtained with the Lambert targeter as a function of time (modified within
 * the function).
 * \param fullProblemResult Map of the cartesian state obtained after propagation of the full dynamics problem as a function of time (modified
 *  within the function).
 * \param cartesianPositionAtDeparture Cartesian position of the body to be propagated at departure [m].
 * \param cartesianPositionAtArrival Cartesian position of the body to be propagated at arrival [m].
 * \param centralBodyGravitationalParameter Gravitational parameter of the central body [m^3 s^-2]. If not provided as input, it is retrieved from
 * the system of bodies.
 */
void propagateLambertTargeterAndFullProblem(
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& centralBody,
        const std::pair< std::string, std::string >& departureAndArrivalBodies ,
        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        std::map< double, Eigen::VectorXd >& dependentVariableResult);

//! Function to propagate the full dynamics problem and the Lambert targeter solution.
/*!
 * Function to propagate the full dynamics problem and the Lambert targeter solution. The function computes the cartesian state as a function of time in two
 * different ways: from the Lambert targeter and from the propagation of the full dynamics problem. The propagator settings are defined inside
 * the function from the propagator type and dependent variables to save provided as inputs.
 * \param timeOfFlight Time of flight [s].
 * \param initialTime Initial time of the propagation [s].
 * \param bodies Body map.
 * \param accelerationModelMap Acceleration map.
 * \param bodyToPropagate Name of the body to be propagated.
 * \param centralBody Name of the central body for the Lambert targeter.
 * \param integratorSettings Integrator settings for the propagation.
 * \param lambertTargeterResult Map of the cartesian state obtained with the Lambert targeter as a function of time (modified within
 * the function).
 * \param fullProblemResult Map of the cartesian state obtained after propagation of the full dynamics problem as a function of time (modified
 *  within the function).
 * \param terminationSphereOfInfluence Boolean denoting whether the propagation stops at the position of the departure and arrival body (false) or at the sphere of influence
 * of the departure and arrival body (true).
 * \param cartesianPositionAtDeparture Cartesian position of the body to be propagated at departure [m].
 * \param cartesianPositionAtArrival Cartesian position of the body to be propagated at arrival [m].
 * \param departureBodyGravitationalParameter Gravitational parameter of the departure body [m^3 s^-2]. If not provided as input, it is retrieved from
 * the system of bodies.
 * \param arrivalBodyGravitationalParameter Gravitational parameter of the arrival body [m^3 s^-2]. If not provided as input, it is retrieved from
 * the system of bodies.
 * \param centralBodyGravitationalParameter Gravitational parameter of the central body [m^3 s^-2]. If not provided as input, it is retrieved from
 * the system of bodies.
 * \param dependentVariablesToSave List of dependent variables to be saved during the full problem propagation.
 * \param propagator Type of propagator to be used for the propagation of the full dynamics problem.
 */
void propagateLambertTargeterAndFullProblem(
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::SystemOfBodies& bodies,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::pair< std::string, std::string >& departureAndArrivalBodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        std::map< double, Eigen::VectorXd >& dependentVariableResult,
        const bool terminationSphereOfInfluence = true,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = std::shared_ptr< DependentVariableSaveSettings >( ),
        const TranslationalPropagatorType propagator = cowell );


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
 * \param bodies Body map.
 * \param accelerationModelMap Acceleration map.
 * \param bodyToPropagate Names of the body to be propagated.
 * \param centralBody  Names of the central body of propagation
 * \param integratorSettings Integrator settings for the propagation.
 * \param terminationSphereOfInfluence Boolean denoting whether the propagation stops at the position of the departure and arrival body (false)
 * or at the sphere of influence of the departure and arrival body (true).
 * \return Pair of vectors containing the difference between the Lambert targeter and the full problem cartesian states
 * (at departure and at arrival respectively).
 */
std::pair< Eigen::Vector6d, Eigen::Vector6d > getDifferenceFullPropagationWrtLambertTargeterAtDepartureAndArrival(
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::SystemOfBodies& bodies,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::vector< std::string >& departureAndArrivalBodies,
        const bool terminationSphereOfInfluence,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = std::shared_ptr< DependentVariableSaveSettings >( ),
        const TranslationalPropagatorType propagator = cowell);



} // namespace propagators

} // namespace tudat

#endif // TUDAT_FULLPROPAGATIONLAMBERTTARGETER_H
