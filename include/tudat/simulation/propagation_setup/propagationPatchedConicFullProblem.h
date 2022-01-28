/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/simulation/simulation.h"
#include "tudat/astro/mission_segments/transferTrajectory.h"

namespace tudat
{

namespace propagators
{

//! Function to get default minimum pericenter radii for a list of bodies
/*!
 * Function to get default minimum pericenter radii for a list of bodies
 * \param bodyNames List of names of bodies for which periapsis radii are to be returned (only the eight planets + Pluto supported)
 * \return List of minimum pericenter radii
 */
std::vector< double > getDefaultMinimumPericenterRadii( const std::vector< std::string >& bodyNames );

//! Function to setup a system of bodies corresponding to the assumptions of a patched conics trajectory,
//! using default ephemerides for the central and transfer bodies.
/*!
 * Function to setup a system of bodies for the patched conics trajectory. The system of bodies contains the central body, the transfer
 * bodies and the body to be propagated. The positions of the central and transfer bodies are directly retrieved from ephemerides.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param nameTransferBodies Vector containing the names of the transfer bodies.
 * \return Body map for the patched conics trajectory.
 */
simulation_setup::SystemOfBodies setupBodyMapFromEphemeridesForPatchedConicsTrajectory(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& nameTransferBodies );


//! Function to setup a system of bodies corresponding to the assumptions of the patched conics trajectory,
//! the ephemerides of the transfer bodies being provided as inputs.
/*!
 * Function to setup a system of bodies for the patched conics trajectory. The system of bodies contains the central body, the transfer
 * bodies and the body to be propagated. Default ephemeris is used for the central body but the ephemerides of the transfer
 * bodies need to be provided as inputs.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param nameTransferBodies Vector containing the names of the transfer bodies.
 * \param ephemerisVectorTransferBodies Vector containing the ephemeris pointers of the different transfer bodies.
 * \param gravitationalParametersTransferBodies Vector containing the gravitational parameters of the transfer bodies [m^3 s^-2].
 * \return Body map for the patched conics trajectory.
 */
simulation_setup::SystemOfBodies setupBodyMapFromUserDefinedEphemeridesForPatchedConicsTrajectory(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& nameTransferBodies,
        const std::vector<ephemerides::EphemerisPointer> &ephemerisVectorTransferBodies,
        const std::vector<double>& gravitationalParametersTransferBodies,
        const std::string& frameOrientation = "J2000" );


//! Function to directly setup a vector of acceleration maps for a patched conics trajectory.
/*!
 * Function to directly setup a vector of acceleration maps for a patched conics trajectory. For each leg, only the central body
 * exerts a point-mass gravity acceleration on the body to be propagated.
 * \param numberOfLegs Number of legs of the patched conics trajectory.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param bodies Body map for the patched conics trajectory.
 * \return Acceleration map for the patched conics trajectory.
 */
std::vector < basic_astrodynamics::AccelerationMap > setupAccelerationMapPatchedConicsTrajectory(
        const double numberOfLegs,
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const simulation_setup::SystemOfBodies& bodies );


//! Function to create the trajectory from the system of bodies.
/*!
 * Function to create the trajectory from the system of bodies.
 * \param bodies Body map from which the trajectory is to be defined.
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
 * \return Trajectory defined from the system of bodies.
 */
transfer_trajectories::Trajectory createTransferTrajectoryObject(
        const simulation_setup::SystemOfBodies& bodies,
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

//! Function to both calculate a patched conics leg without DSM and propagate the full dynamics problem.
/*!
 * Function to both calculate a patched conics leg without DSM and propagate the full dynamics problem.
 * \param bodies Body map for the patched conics leg.
 * \param departureAndArrivalBodies Vector containing the names of the departure and arrival bodies of the leg.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param cartesianPositionAtDeparture Cartesian position of the body to be propagated at the leg departure [m].
 * \param cartesianPositionAtArrival Cartesian position of the body to be propagated at the leg arrival [m].
 * \param initialTime Time at departure [s].
 * \param timeOfFlight Time of flight for the leg [s].
 * \param propagatorSettings Propagator settings for the propagation of the full dynamics problem.
 * \param integratorSettings Integration settings for the full problem propagation.
 * \param patchedConicsResult Patched conics solution for the leg.
 * \param fullProblemResult propagation results of the full problem over the leg.
 */
void propagateMgaWithoutDsmAndFullProblem(
        simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string > departureAndArrivalBodies,
        const std::string centralBody,
        const Eigen::Vector3d cartesianPositionAtDeparture,
        const Eigen::Vector3d cartesianPositionAtArrival,
        const double initialTime,
        const double timeOfFlight,
        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        std::map< double, Eigen::Vector6d >& patchedConicsResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        std::map< double, Eigen::VectorXd >& dependentVariableResultCurrentLeg );

//! Function to both calculate a patched conics leg including a DSM and propagate the corresponding full dynamics problem.
/*!
 * Function to both calculate a patched conics leg including a DSM and propagate the corresponding full dynamics problem. The patched
 * conics leg with DSM is calculated using the velocity formulation.
 * \param bodies Body map for the patched conics leg.
 * \param departureAndArrivalBodies Vector containing the names of the departure and arrival bodies of the leg.
 * \param dsm Name of the DSM.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param cartesianPositionAtDeparture Cartesian position of the body to be propagated at the leg departure [m].
 * \param cartesianPositionDSM Cartesian position of the body to be propagated at the DSM location [m].
 * \param cartesianPositionAtArrival Cartesian position of the body to be propagated at the leg arrival [m].
 * \param initialTime Time at departure [s].
 * \param timeDsm Time at which the DSM is to be performed [s].
 * \param timeArrival Time at arrival [s].
 * \param legType Type of the leg.
 * \param trajectoryVariableVector Trajectory variable vector characterising the leg.
 * \param semiMajorAxis Semi-major axis at trajectory departure (only used for a departure leg and not a swing-by one) [m].
 * \param eccentricity Eccentricity at trajectory departure (only used for a departure leg and not a swing-by one).
 * \param velocityAfterDeparture Velocity coordinates of the body to be propagated just after the swing-by it has performed about the
 * departure body of the leg [m/s].
 * \param velocityBeforeArrival  Velocity coordinates of the body to be propagated just before it reaches the arrival body of the leg [m/s].
 * \param propagatorSettingsBeforeDsm Propagator settings for the full problem propagation from the departure body to the DSM location.
 * \param propagatorSettingsAfterDsm propagators settings for the full problem propagation from the DSM location to the arrival body.
 * \param integratorSettings Integrator settings for the propagation of the full dynamics problem.
 * \param patchedConicsResultFromDepartureToDsm Patched conics solution for the first part of the leg (from departure body to DSM).
 * \param fullProblemResultFromDepartureToDsm propagation results of the full dynamics problem for the first part of the leg (from
 * departure body to DSM).
 * \param patchedConicsResultFromDsmToArrival Patched conics solution for the second part of the leg (from DSM to arrival body).
 * \param fullProblemResultFromDsmToArrival propagation results of the full dynamics problem for the second part of the leg (from DSM
 * to arrival body).
 */
void propagateMga1DsmVelocityAndFullProblem(
        simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string > departureAndArrivalBodies,
        const std::string& dsm,
        const std::string& centralBody,
        const Eigen::Vector3d cartesianPositionAtDeparture,
        const Eigen::Vector3d cartesianPositionDSM,
        const Eigen::Vector3d cartesianPositionAtArrival,
        const double initialTime,
        const double timeDsm,
        const double timeArrival,
        const transfer_trajectories::TransferLegType& legType,
        const std::vector< double >& trajectoryVariableVector,
        const double semiMajorAxis,
        const double eccentricity,
        Eigen::Vector3d& velocityAfterDeparture,
        Eigen::Vector3d& velocityBeforeArrival,
        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettingsBeforeDsm,
        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettingsAfterDsm,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        std::map< double, Eigen::Vector6d >& patchedConicsResultFromDepartureToDsm,
        std::map< double, Eigen::Vector6d >& fullProblemResultFromDepartureToDsm,
        std::map< double, Eigen::VectorXd >& dependentVariablesFromDepartureToDsm,
        std::map< double, Eigen::Vector6d >& patchedConicsResultFromDsmToArrival,
        std::map< double, Eigen::Vector6d >& fullProblemResultFromDsmToArrival,
        std::map< double, Eigen::VectorXd >& dependentVariablesFromDsmToArrival );



//! Function to both calculate a patched conics leg including a DSM and propagate the corresponding full dynamics problem.
/*!
 * Function to both calculate a patched conics leg including a DSM and propagate the corresponding full dynamics problem. The patched
 * conics leg with DSM is calculated using the position formulation.
 * \param bodies bodies Body map for the patched conics leg.
 * \param departureAndArrivalBodies Vector containing the names of the departure and arrival bodies of the leg.
 * \param dsm Name of the DSM.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param cartesianPositionAtDeparture Cartesian position of the body to be propagated at the leg departure [m].
 * \param cartesianPositionDSM Cartesian position of the body to be propagated at the DSM location [m].
 * \param cartesianPositionAtArrival Cartesian position of the body to be propagated at the leg arrival [m].
 * \param initialTime Time at departure [s].
 * \param timeDsm Time at which the DSM is to be performed [s].
 * \param timeArrival Time at arrival [s].
 * \param legType Type of the leg.
 * \param trajectoryVariableVector Trajectory variable vector characterising the leg.
 * \param minimumPericenterRadius Minimum pericenter radius (only used for a swing-by leg and not a trajectory departure one).
 * \param semiMajorAxis Semi-major axis at trajectory departure (only used for a departure leg and not a swing-by one) [m].
 * \param eccentricity Eccentricity at trajectory departure (only used for a departure leg and not a swing-by one).
 * \param velocityAfterDeparture Velocity coordinates of the body to be propagated just after the swing-by it has performed about the
 * departure body of the leg [m/s].
 * \param velocityBeforeArrival Velocity coordinates of the body to be propagated just before it reaches the arrival body of the leg [m/s].
 * \param propagatorSettingsBeforeDsm Propagator settings for the full problem propagation from the departure body to the DSM location.
 * \param propagatorSettingsAfterDsm propagators settings for the full problem propagation from the DSM location to the arrival body.
 * \param integratorSettings Integrator settings for the propagation of the full dynamics problem.
 * \param patchedConicsResultFromDepartureToDsm Patched conics solution for the first part of the leg (from departure body to DSM).
 * \param fullProblemResultFromDepartureToDsm propagation results of the full dynamics problem for the first part of the leg (from
 * departure body to DSM).
 * \param patchedConicsResultFromDsmToArrival Patched conics solution for the second part of the leg (from DSM to arrival body).
 * \param fullProblemResultFromDsmToArrival propagation results of the full dynamics problem for the second part of the leg (from DSM
 * to arrival body).
 */
void propagateMga1DsmPositionAndFullProblem(
        simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string > departureAndArrivalBodies,
        const std::string& dsm,
        const std::string& centralBody,
        const Eigen::Vector3d cartesianPositionAtDeparture,
        const Eigen::Vector3d cartesianPositionDSM,
        const Eigen::Vector3d cartesianPositionAtArrival,
        const double initialTime,
        const double timeDsm,
        const double timeArrival,
        const transfer_trajectories::TransferLegType& legType,
        const std::vector< double >& trajectoryVariableVector,
        const double minimumPericenterRadius,
        const double semiMajorAxis,
        const double eccentricity,
        Eigen::Vector3d& velocityAfterDeparture,
        Eigen::Vector3d& velocityBeforeArrival,
        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettingsBeforeDsm,
        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettingsAfterDsm,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        std::map< double, Eigen::Vector6d >& patchedConicsResultFromDepartureToDsm,
        std::map< double, Eigen::Vector6d >& fullProblemResultFromDepartureToDsm,
        std::map< double, Eigen::VectorXd >& dependentVariablesFromDepartureToDsm,
        std::map< double, Eigen::Vector6d >& patchedConicsResultFromDsmToArrival,
        std::map< double, Eigen::Vector6d >& fullProblemResultFromDsmToArrival,
        std::map< double, Eigen::VectorXd >& dependentVariablesFromDsmToArrival );

std::shared_ptr< propagators::PropagationTerminationSettings > getSingleLegPartSphereOfInfluenceTerminationSettings(
        simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::string& departureBody,
        const std::string& arrivalBody,
        const double initialTimeCurrentLeg,
        const double finalTimeCurrentLeg,
        const bool useBackwardLegToDepartureBody,
        const double terminationDistanceScaler = 1.0 );

std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
std::shared_ptr< propagators::PropagationTerminationSettings > > getSingleLegSphereOfInfluenceTerminationSettings(
        simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::string& departureBody,
        const std::string& arrivalBody,
        const double initialTimeCurrentLeg,
        const double finalTimeCurrentLeg,
        const double terminationDistanceScaler = 1.0 );

//! Function to calculate the patched conics trajectory and to propagate the corresponding full problem.
std::vector< std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > > getPatchedConicPropagatorSettings(
        simulation_setup::SystemOfBodies& bodies,
        const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const std::vector< double >& trajectoryVariableVector,
        const std::vector< double >& minimumPericenterRadiiVector,
        const std::vector< double >& semiMajorAxesVector,
        const std::vector< double >& eccentricitiesVector,
        const std::vector< std::shared_ptr< DependentVariableSaveSettings > > dependentVariablesToSave,
        const TranslationalPropagatorType propagator,
        const bool terminationSphereOfInfluence,
        const double terminationDistanceScaler = 1.0 );


//! Function to propagate the motion of a body over a trajectory leg, both along a keplerian orbit and in a full dynamics problem.
/*!
 * Function to propagate the motion of a body over a trajectory leg, both along a keplerian orbit and in a full dynamics problem.
 * \param timeOfFlight Time of flight during which the body has to be propagated [s]
 * \param initialTime Initial time [s].
 * \param bodies Body map defining the problem.
 * \param centralBody Name of the central body of the keplerian trajectory.
 * \param departureAndArrivalBodies Name of the departure and arrival bodies of the leg.
 * \param velocityAfterDeparture Velocity coordinates of the body to be propagated just after the swing-by it has performed about the
 * departure body of the leg [m/s].
 * \param propagatorSettings Propagator settings for the propagation of the full dynamics problem.
 * \param integratorSettings Integrator settings for the propagation of the full dynamics problem.
 * \param keplerianOrbitResult Keplerian orbit solution.
 * \param fullProblemResult propagation results of the full problem.
 * \param centralBodyGravitationalParameter Gravitational parameter of the central body [m^3 s^-2]. If not provided as input, it is retrieved
 * from the system of bodies.
 * \param cartesianPositionAtDeparture Cartesian position of the body to the propagated at the leg departure [m]. If not provided as input, it
 * is retrieved from the ephemerides.
 */
void propagateKeplerianOrbitLegAndFullProblem(
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& centralBody,
        const std::vector<std::string>& departureAndArrivalBodies,
        const Eigen::Vector3d& velocityAfterDeparture,
        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& keplerianOrbitResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        std::map< double, Eigen::VectorXd >& dependentVariables,
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& cartesianPositionAtDeparture);

//! Function to calculate the patched conics trajectory and to propagate the corresponding full problem.
/*!
 * Function to calculate the patched conics trajectory and to propagate the corresponding full problem. The propagator settings to be used
 * for the propagation of the full dynamics problem are directly provided as inputs.
 * \param bodies Body map for the patched conics trajectory.
 * \param transferBodyOrder Vector containing the names of the transfer bodies involved in the trajectory.
 * \param patchedConicCentralBody Name of the central body of the patched conics trajectory.
 * \param legTypeVector Vector containing the leg types.
 * \param trajectoryVariableVector Vector containing all the defining variables for the whole trajectory.
 * \param minimumPericenterRadiiVector Vector containing the minimum distance between the spacecraft and the body.
 * \param semiMajorAxesVector Vector containing the semi-major axes of the departure and arrival legs.
 * \param eccentricitiesVector Vector containing the eccentricities of the departure and arrival legs.
 * \param propagatorSettings Vector containing a pair of propagator settings for each leg (one for the backward and the other for the
 * forward propagation)
 * \param integratorSettings Integrator settings for the propagation of the full problem.
 * \param patchedConicsResultForEachLeg Patched conics solution along each leg.
 * \param fullProblemResultForEachLeg Full problem propagation results along each leg.
 */
void fullPropagationPatchedConicsTrajectory(
        simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& patchedConicCentralBody,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const std::vector< double >& trajectoryVariableVector,
        const std::vector< double >& minimumPericenterRadiiVector,
        const std::vector< double >& semiMajorAxesVector,
        const std::vector< double >& eccentricitiesVector,
        const std::vector< std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > > propagatorSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        std::map< int, std::map< double, Eigen::Vector6d > >& patchedConicsResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg,
        std::map< int, std::map< double, Eigen::VectorXd > >& dependentVariableResultForEachLeg);



//! Function to calculate the patched conics trajectory and to propagate the corresponding full problem.
/*!
 * Function to calculate the patched conics trajectory and to propagate the corresponding full problem. The propagator settings for the
 * full problem propagation are defined inside the function from the propagator type and dependent variables to save
 * which are provided as inputs.
 * \param bodies Body map for the patched conics trajectory.
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
 * \param patchedConicsResultForEachLeg Patched conics solution along each leg.
 * \param fullProblemResultForEachLeg Full problem propagation results along each leg.
 * \param terminationSphereOfInfluence Boolean denoting whether the propagation stops at the exact position (false) or at the sphere of
 * influence (true) of the departure and arrival body of each leg of the trajectory. The default value is false.
 * \param dependentVariablesToSave Vector containing the dependent variables to be saved during the full problem propagation for each leg.
 * \param propagator Type of propagator to be used for the full problem propagation.
 */
void fullPropagationPatchedConicsTrajectory(
        simulation_setup::SystemOfBodies& bodies,
        const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const std::vector< double >& trajectoryVariableVector,
        const std::vector< double >& minimumPericenterRadiiVector,
        const std::vector< double >& semiMajorAxesVector,
        const std::vector< double >& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        std::map< int, std::map< double, Eigen::Vector6d > >& patchedConicsResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg,
        std::map< int, std::map< double, Eigen::VectorXd > >& dependentVariableResultForEachLeg,
        const bool terminationSphereOfInfluence = false,
        const std::vector< std::shared_ptr< DependentVariableSaveSettings > > dependentVariablesToSave =
        std::vector < std::shared_ptr< DependentVariableSaveSettings > >( ),
        const TranslationalPropagatorType propagator = cowell);



//! Function to calculate the patched conics trajectory and to propagate the corresponding full problem
//! with the same acceleration map for every leg.
/*!
 * Function to calculate the patched conics trajectory and to propagate the corresponding full problem with the same acceleration map for
 * every leg. The propagator settings for the full problem propagation are defined inside the function from the propagator type and
 * dependent variables to save provided as inputs.
 * \param bodies Body map for the patched conics trajectory.
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
 * \param patchedConicsResultForEachLeg Patched conics solution along each leg.
 * \param fullProblemResultForEachLeg Full problem propagation results along each leg.
 * \param terminationSphereOfInfluence Boolean denoting whether the propagation stops at the exact position (false) or at the sphere of
 * influence (true) of the departure and arrival body of each leg of the trajectory. The default value is false.
 * \param dependentVariablesToSave Vector containing the dependent variables to be saved during the full problem propagation for each leg.
 * \param propagator Type of propagator to be used for the full problem propagation.
 */
void fullPropagationPatchedConicsTrajectory(
        simulation_setup::SystemOfBodies& bodies,
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
        std::map< int, std::map< double, Eigen::Vector6d > >& patchedConicsResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg,
        std::map< int, std::map< double, Eigen::VectorXd > >& dependentVariableResultForEachLeg,
        const bool terminationSphereOfInfluence = false,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = std::shared_ptr< DependentVariableSaveSettings > ( ),
        const TranslationalPropagatorType propagator = cowell );


//! Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem,
//! at both departure and arrival positions for each leg.
/*!
 * Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem, at both
 * departure and arrival positions for each leg. The propagator settings to be used for the propagation of the full dynamics
 * problem are directly provided as inputs.
 * \param bodies Body map for the patched conics trajectory.
 * \param transferBodyOrder Vector containing the names of the transfer bodies involved in the trajectory.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param legTypeVector Vector containing the leg types.
 * \param trajectoryVariableVector Vector containing all the defining variables for the whole trajectory.
 * \param minimumPericenterRadiiVector Vector containing the minimum distance between the spacecraft and the body.
 * \param semiMajorAxesVector Vector containing the semi-major axes of the departure and arrival legs.
 * \param eccentricitiesVector Vector containing the eccentricities of the departure and arrival legs.
 * \param propagatorSettings Vector containing a pair of propagator settings for each leg (one for the backward and the other for the
 * forward propagation)
 * \param integratorSettings Integrator settings for the propagation of the full problem.
 * \return Map of vector pairs. Each vector pair contains the difference in cartesian state between patched conics trajectory and full problem for a given leg,
 * at departure and arrival respectively.
 */
std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullProblemWrtPatchedConicsTrajectory(
        simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const std::vector<double>& trajectoryVariableVector,
        const std::vector<double>& minimumPericenterRadiiVector,
        const std::vector<double>& semiMajorAxesVector,
        const std::vector<double>& eccentricitiesVector,
        const std::vector< std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > > propagatorSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings);




//! Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem,
//! at both departure and arrival positions for each leg.
/*!
 * Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem, at both
 * departure and arrival positions for each leg. The propagator settings for the full problem propagation are defined inside the
 * function from the propagator type and dependent variables to save provided as inputs.
 * \param bodies Body map for the patched conics trajectory.
 * \param accelerationMap Vector of acceleration maps for each leg, to be used to propagate the full problem.
 * \param transferBodyOrder Vector containing the names of the transfer bodies involved in the trajectory.
 * \param centralBody Name of the central body of the patched conics trajectory.
 * \param bodyToPropagate Name of the body to be propagated.
 * \param legTypeVector Vector containing the leg types.
 * \param trajectoryVariableVector Vector containing all the defining variables for the whole trajectory.
 * \param minimumPericenterRadiiVector Vector containing the minimum distance between the spacecraft and the body.
 * \param semiMajorAxesVector Vector containing the semi-major axes of the departure and arrival legs.
 * \param eccentricitiesVector Vector containing the eccentricities of the departure and arrival legs.
 * \param integratorSettings Integrator settings for the propagation of the full problem.
 * \param terminationSphereOfInfluence Boolean denoting whether the propagation stops at the exact position (false) or at the sphere of
 * influence (true) of the departure and arrival body of each leg of the trajectory. The default value is false.
 * \param dependentVariablesToSave Vector containing the dependent variables to be saved during the full problem propagation for each leg.
 * \param propagator Type of propagator to be used for the full problem propagation.
 * \return Map of vector pairs. Each vector pair contains the difference in cartesian state between patched conics trajectory and full problem for a given leg,
 * at departure and arrival respectively.
 */
std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullProblemWrtPatchedConicsTrajectory(
        simulation_setup::SystemOfBodies& bodies,
        const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const std::vector<double>& trajectoryVariableVector,
        const std::vector<double>& minimumPericenterRadiiVector,
        const std::vector<double>& semiMajorAxesVector,
        const std::vector<double>& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        const bool terminationSphereOfInfluence = false,
        const std::vector< std::shared_ptr< DependentVariableSaveSettings > > dependentVariablesToSave =
        std::vector < std::shared_ptr< DependentVariableSaveSettings > >( ),
        const TranslationalPropagatorType propagator = cowell);



//! Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem,
//! at both departure and arrival positions for each leg, using the same accelerations for each leg.
/*!
 * Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem, at both departure
 * and arrival positions for each leg when using the same accelerations for each leg. The propagator settings for the full problem
 * propagation are defined inside the function from the propagator type and dependent variables to save provided as inputs.
 * \param bodies Body map for the patched conics trajectory.
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
 * \param terminationSphereOfInfluence Boolean denoting whether the propagation stops at the exact position (false) or at the sphere of
 * influence (true) of the departure and arrival body of each leg of the trajectory. The default value is false.
 * \param dependentVariablesToSave Vector containing the dependent variables to be saved during the full problem propagation for each leg.
 * \param propagator Type of propagator to be used for the full problem propagation.
 * \return Map of vector pairs. Each vector pair contains the difference in cartesian state between patched conics trajectory and full problem for a given leg,
 * at departure and arrival respectively.
 */
std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullProblemWrtPatchedConicsTrajectory(
        simulation_setup::SystemOfBodies& bodies,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const std::vector<double>& trajectoryVariableVector,
        const std::vector<double>& minimumPericenterRadiiVector,
        const std::vector<double>& semiMajorAxesVector,
        const std::vector<double>& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        const bool terminationSphereOfInfluence,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = std::shared_ptr< DependentVariableSaveSettings > ( ) ,
        const TranslationalPropagatorType propagator = cowell);

}

}


#endif // TUDAT_PROPAGATION_PATCHED_CONIC_FULL
