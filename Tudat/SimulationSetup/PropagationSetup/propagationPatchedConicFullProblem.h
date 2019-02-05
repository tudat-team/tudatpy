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


std::vector< double > getDefaultMinimumPericenterRadii( const std::vector< std::string >& bodyNames );

//! Function to setup a body map corresponding to the assumptions of a patched conics trajectory,
//! using default ephemerides for the central and transfer bodies.
/*!
 * Function to setup a body map for the patched conics trajectory. The body map contains the central body, the transfer
 * bodies and the body to be propagated. The positions of the central and transfer bodies are directly retrieved from ephemerides.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param nameTransferBodies Vector containing the names of the transfer bodies.
 * \return Body map for the patched conics trajectory.
 */
simulation_setup::NamedBodyMap setupBodyMapFromEphemeridesForPatchedConicsTrajectory(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& nameTransferBodies );


//! Function to setup a body map corresponding to the assumptions of the patched conics trajectory,
//! the ephemerides of the transfer bodies being provided as inputs.
/*!
 * Function to setup a body map for the patched conics trajectory. The body map contains the central body, the transfer
 * bodies and the body to be propagated. Default ephemeris is used for the central body but the ephemerides of the transfer
 * bodies need to be provided as inputs.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param nameTransferBodies Vector containing the names of the transfer bodies.
 * \param ephemerisVectorTransferBodies Vector containing the ephemeris pointers of the different transfer bodies.
 * \param gravitationalParametersTransferBodies Vector containing the gravitational parameters of the transfer bodies [m^3 s^-2].
 * \return Body map for the patched conics trajectory.
 */
simulation_setup::NamedBodyMap setupBodyMapFromUserDefinedEphemeridesForPatchedConicsTrajectory(const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& nameTransferBodies,
        const std::vector<ephemerides::EphemerisPointer> &ephemerisVectorTransferBodies,
        const std::vector<double>& gravitationalParametersTransferBodies);


//! Function to directly setup a vector of acceleration maps for a patched conics trajectory.
/*!
 * Function to directly setup a vector of acceleration maps for a patched conics trajectory. For each leg, only the central body
 * exerts a point-mass gravity acceleration on the body to be propagated.
 * \param numberOfLegs Number of legs of the patched conics trajectory.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param bodyMap Body map for the Lambert targeter.
 * \return Acceleration map for the Lambert targeter.
 */
std::vector < basic_astrodynamics::AccelerationMap > setupAccelerationMapPatchedConicsTrajectory(
        const double numberOfLegs,
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const simulation_setup::NamedBodyMap& bodyMap );


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


void fullPropagationMGA(
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const int numberOfLegs,
        const std::vector< std::string >& transferBodyOrder,
        const std::vector< std::string >& bodiesAndManoeuvresOrder,
        const std::vector< std::string >& centralBody,
        const std::vector< std::string >& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const Eigen::VectorXd& trajectoryVariableVector,
        const Eigen::VectorXd& minimumPericenterRadiiVector,
        const Eigen::VectorXd& semiMajorAxesVector,
        const Eigen::VectorXd& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg);


std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullPropagationWrtLambertTargeterMGA(
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const int numberOfLegs,
        const std::vector< std::string >& transferBodyOrder,
        const std::vector< std::string >& bodiesAndManoeuvresOrder,
        const std::vector< std::string >& centralBody,
        const std::vector< std::string >& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const Eigen::VectorXd& trajectoryVariableVector,
        const Eigen::VectorXd& minimumPericenterRadiiVector,
        const Eigen::VectorXd& semiMajorAxesVector,
        const Eigen::VectorXd& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings);



}

}


#endif // TUDAT_PROPAGATION_PATCHED_CONIC_FULL
