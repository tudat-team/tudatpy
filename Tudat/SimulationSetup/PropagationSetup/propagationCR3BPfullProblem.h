/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATION_CR3BP_FULL
#define TUDAT_PROPAGATION_CR3BP_FULL

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"


namespace tudat
{

namespace propagators
{

//! Setup CR3BP body map.
/*!
 * Setup CR3BP body map. The two primaries are declared in the body map. They are in circular orbit about their barycenter, orbiting it with
 * the same mean motion and thus stay aligned during the propagation. The third, smaller body to be propagated in the CR3BP is defined as well.
 * \param distancePrimarySecondary Distance between primaries.   [m]
 * \param namePrimaryBody Name of the primary body.
 * \param nameSecondaryBody Name of the secondary body.
 * \param nameBodyToPropagate Name of the third, massless body to be propagated.
 * \return Body Map modelling the CR3BP.
 */
simulation_setup::NamedBodyMap setupBodyMapCR3BPBodyMap(
        const double distancePrimarySecondary,
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate);



//! Setup CR3BP acceleration map.
/*!
 * Setup CR3BP acceleration map. Define the acceleration map for the CR3BP, from the corresponding body map.
 * \param namePrimaryBody Name of the primary body.
 * \param nameSecondaryBody Name of the secondary body.
 * \param nameBodyToPropagate Name of the third, massless body to be propagated.
 * \param bodiesToPropagate Bodies to be propagated.
 * \param centralBodies Central bodies for the propagation.
 * \param bodyMap CR3BP body map.
 * \return Acceleration map for the CR3BP.
 */
basic_astrodynamics::AccelerationMap setupAccelerationMapCR3BP(
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::NamedBodyMap& bodyMap );



//! Propagate CR3BP from CR3BP environment
/*!
 * Propagate CR3BP from CR3BP environment.
 * \param initialTime Initial time for the propagation.  [s]
 * \param finalPropagationTime Final time at which the propagation ends. [s]
 * \param initialState Initial state of the third, massless body to be propagated (x-position coordinate [m], y-position coordinate [m],
 *   z-position coordinate [m], x-velocity coordinate [m/s], y-velocity coordinate [m/s], z-velocity coordinate [m/s]).
 * \param integratorSettings Integrator settings for the propagation (the initial time and time-step have to be defined as dimensional times,
 * their conversion to adimensional times is implemented within the function).
 * \param bodyMap Body map for the CR3BP.
 * \param bodiesCR3BP Name of the two primaries defining the CR3BP.
 * \return Third body state after propagation in the CR3BP, converted to inertial cartesian coordinates.
 *          (x-position coordinate                [m],
 *          y-position coordinate                 [m],
 *          z-position coordinate                 [m],
 *          x-velocity coordinate               [m/s],
 *          y-velocity coordinate               [m/s],
 *          z-velocity coordinate               [m/s]).
 */
Eigen::Vector6d propagateCR3BPFromEnvironment(
        const double initialTime,
        const double finalPropagationTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        simulation_setup::NamedBodyMap& bodyMap,
        const std::vector < std::string >& bodiesCR3BP );



//! Propagate the CR3BP and the full dynamics problem and compute the state difference at the end of the propagation.
/*!
 * Propagate the CR3BP and the full dynamics problem, and compute the state difference at the end of the propagation.
 * \param initialTime Initial time for the propagation [s].
 * \param finalTime Final time at which the propagation ends [s].
 * \param initialState Initial state of the third, massless body to be propagated (x-position coordinate [m], y-position coordinate [m],
 *   z-position coordinate [m], x-velocity coordinate [m/s], y-velocity coordinate [m/s], z-velocity coordinate [m/s]).
 * \param integratorSettings Integrator settings for the propagation (the initial time and time-step have to be defined as dimensional times,
 * their conversion to adimensional times is implemented within the function).
 * \param accelerationModelMap Acceleration map for the CR3BP.
 * \param bodiesToPropagate Bodies to be propagated.
 * \param centralBodies Central bodies for the propagation.
 * \param bodyMap Body Map for the CR3BP.
 * \param bodiesCR3BP Name of the two primaries defining the CR3BP.
 * \return State difference between the full dynamics problem and the CR3BP at the final propagation time (expressed in inertial cartesian coordinates).
 */
Eigen::Vector6d propagateCR3BPandFullDynamicsProblem(
        const double initialTime,
        const double finalTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        simulation_setup::NamedBodyMap& bodyMap,
        const std::vector < std::string >& bodiesCR3BP );

}

}


#endif // TUDAT_PROPAGATION_CR3BP_FULL
