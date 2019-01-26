/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATION_PATCHED_CONIC_FULL
#define TUDAT_PROPAGATION_PATCHED_CONIC_FULL

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"

namespace tudat
{

namespace propagators
{

//! Function to setup a body map corresponding to the assumptions of the patched conics trajectory,
//! retrieving positions of the transfer bodies from ephemerides.
simulation_setup::NamedBodyMap setupBodyMapFromEphemeridesForPatchedConicsTrajectory(const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& nameTransferBodies );

//! Function to setup a body map corresponding to the assumptions of the patched conics trajectory,
//! the positions of the transfer bodies being provided as inputs.
simulation_setup::NamedBodyMap setupBodyMapFromUserDefinedStatesForPatchedConicsTrajectory(const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& nameTransferBodies,
        const std::vector<ephemerides::EphemerisPointer> &ephemerisVector, const std::vector<double>& gravitationalParametersTransferBodies);


//! Function to create the trajectory from the body map.
/*!
 * Function to create the trajectory from the body map.
 * \param bodyMap Body map from which the trajectory is to be defined.
 * \param transferBodyOrder Vector containing the names of the transfer bodies involved in the trajectory.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param transferLegTypes Vector containing the leg types.
 * \param trajectoryIndependentVariables Vector containing all the defining variables for the whole trajectory.
 * \param minimumPericenterRadii Vector containing the minimum distance between the spacecraft and the body.
 * \param includeDepartureDeltaV Boolean denoting whether to include the Delta V at departure.
 * \param departureSemiMajorAxis Semi-major axis of the departure leg.
 * \param departureEccentricity Eccentricity of the departure leg.
 * \param includeArrivalDeltaV Boolean denoting whether to include the Delta V at arrival.
 * \param arrivalSemiMajorAxis Semi-major axis of the arrival leg.
 * \param arrivalEccentricity Eccentricity of the arrival leg.
 * \return Trajectory defined from the body map.
 */
transfer_trajectories::Trajectory createTransferTrajectoryObject(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::vector< transfer_trajectories::TransferLegType >& transferLegTypes,
        const std::vector< double > trajectoryIndependentVariables,
        const std::vector< double > minimumPericenterRadii,
        const bool includeDepartureDeltaV,
        const double departureSemiMajorAxis,
        const double departureEccentricity,
        const bool includeArrivalDeltaV,
        const double arrivalSemiMajorAxis,
        const double arrivalEccentricity );



//! Function to calculate the patched conics trajectory and to propagate the corresponding full problem.
/*!
 * Function to calculate the patched conics trajectory and to propagate the corresponding full problem.
 * \param bodyMap Body map for the patched conics trajectory.
 * \param accelerationMap Vector of acceleration maps for each leg to propagate the full problem.
 * \param transferBodyOrder Vector containing the names of the transfer bodies involved in the trajectory.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param bodyToPropagate Name of the body to be propagated.
 * \param legTypeVector Vector containing the leg types.
 * \param trajectoryVariableVector Vector containing all the defining variables for the whole trajectory.
 * \param minimumPericenterRadiiVector Vector containing the minimum distance between the spacecraft and the body.
 * \param semiMajorAxesVector Vector containing the semi-major axes of the departure and arrival legs.
 * \param eccentricitiesVector Vector containing the eccentricities of the departure and arrival legs.
 * \param integratorSettings Integrator settings for the propagation of the full problem.
 * \param lambertTargeterResultForEachLeg Lambert Targeter results along each leg.
 * \param fullProblemResultForEachLeg Full problem propagation results along each leg.
 */
void fullPropagationPatchedConicsTrajectory(
        simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const std::vector<double>& trajectoryVariableVector,
        const std::vector<double>& minimumPericenterRadiiVector,
        const std::vector<double>& semiMajorAxesVector,
        const std::vector<double>& eccentricitiesVector,
        const std::shared_ptr<numerical_integrators::IntegratorSettings<double> >& integratorSettings,
        std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg);


//! Function to calculate the patched conics trajectory and to propagate the corresponding full problem
//! with the same acceleration map for every leg.
/*!
 * Function to calculate the patched conics trajectory and to propagate the corresponding full problem with the same acceleration map for every leg.
 * \param bodyMap Body map for the patched conics trajectory.
 * \param accelerationMap Acceleration map (to be used for every leg) to propagate the full problem.
 * \param transferBodyOrder Vector containing the names of the transfer bodies involved in the trajectory.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param bodyToPropagate Name of the body to be propagated.
 * \param legTypeVector Vector containing the leg types.
 * \param trajectoryVariableVector Vector containing all the defining variables for the whole trajectory.
 * \param minimumPericenterRadiiVector Vector containing the minimum distance between the spacecraft and the body.
 * \param semiMajorAxesVector Vector containing the semi-major axes of the departure and arrival legs.
 * \param eccentricitiesVector Vector containing the eccentricities of the departure and arrival legs.
 * \param integratorSettings Integrator settings for the propagation of the full problem.
 * \param lambertTargeterResultForEachLeg Lambert Targeter results along each leg.
 * \param fullProblemResultForEachLeg Full problem propagation results along each leg.
 */
void fullPropagationPatchedConicsTrajectorySingleAccelerationMap(
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const std::vector< double >& trajectoryVariableVector,
        const std::vector< double >& minimumPericenterRadiiVector,
        const std::vector< double >& semiMajorAxesVector,
        const std::vector< double >& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg);



//! Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem,
//! at both departure and arrival positions for each leg.
/*!
 * Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem, at both departure and arrival positions for each leg.
 * \param bodyMap Body map for the patched conics trajectory.
 * \param accelerationMap Vector of acceleration maps for each leg to propagate the full problem.
 * \param transferBodyOrder Vector containing the names of the transfer bodies involved in the trajectory.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param bodyToPropagate Name of the body to be propagated.
 * \param legTypeVector Vector containing the leg types.
 * \param trajectoryVariableVector Vector containing all the defining variables for the whole trajectory.
 * \param minimumPericenterRadiiVector Vector containing the minimum distance between the spacecraft and the body.
 * \param semiMajorAxesVector Vector containing the semi-major axes of the departure and arrival legs.
 * \param eccentricitiesVector Vector containing the eccentricities of the departure and arrival legs.
 * \param integratorSettings Integrator settings for the propagation of the full problem.
 * \return Map of vector pairs. Each vector pair contains the difference in cartesian state between patched conics trajectory and full problem for a given leg,
 * at departure and arrival respectively.
 */
std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullProblemWrtPatchedConicsTrajectory(
        simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const std::vector< double >& trajectoryVariableVector,
        const std::vector< double >& minimumPericenterRadiiVector,
        const std::vector< double >& semiMajorAxesVector,
        const std::vector< double >& eccentricitiesVector,
        const std::shared_ptr<numerical_integrators::IntegratorSettings<double> > &integratorSettings);



//! Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem,
//! at both departure and arrival positions for each leg when using the same accelerations for each leg.
/*!
 * Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem, at both departure and arrival positions for each leg when using the same accelerations for each leg.
 * \param bodyMap Body map for the patched conics trajectory.
 * \param accelerationMap Acceleration map (to be used for every leg) to propagate the full problem.
 * \param transferBodyOrder Vector containing the names of the transfer bodies involved in the trajectory.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param bodyToPropagate Name of the body to be propagated.
 * \param legTypeVector Vector containing the leg types.
 * \param trajectoryVariableVector Vector containing all the defining variables for the whole trajectory.
 * \param minimumPericenterRadiiVector Vector containing the minimum distance between the spacecraft and the body.
 * \param semiMajorAxesVector Vector containing the semi-major axes of the departure and arrival legs.
 * \param eccentricitiesVector Vector containing the eccentricities of the departure and arrival legs.
 * \param integratorSettings Integrator settings for the propagation of the full problem.
 * \return Map of vector pairs. Each vector pair contains the difference in cartesian state between patched conics trajectory and full problem for a given leg,
 * at departure and arrival respectively.
 */
std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullProblemWrtPatchedConicsTrajectoryWithSingleAccelerationMap(
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const std::vector<double>& trajectoryVariableVector,
        const std::vector<double>& minimumPericenterRadiiVector,
        const std::vector<double>& semiMajorAxesVector,
        const std::vector<double>& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings);


}

}


#endif // TUDAT_PROPAGATION_PATCHED_CONIC_FULL
