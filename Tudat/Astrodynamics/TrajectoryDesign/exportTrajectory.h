/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      More information on the trajectory design code and how the quantities
 *      are calculated can be found in [Musegaas, 2012], who is also the author
 *      of this code
 *
*/

// This contains some functions to facilitate making plots of trajectories.

#include <vector>

#include <Eigen/Core>

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Return a vector of positions and times corresponding to a trajectory at a certain epoch.
/*!
 * \brief returns the trajectory of the spacecraft.
 * \param initialCartesianState initial cartesian state of the spacecraft.
 * \param centralBodyGravitationalParameter central body gravitational paremeter.
 * \param duration duration of the exported trajectory
 * \param maximumTimeStep maximum time step between the positions that are returned
 * \param positionVector vector of positions along the trajectory.
 * \param timeVector vector of times corresponding to the positions.
 * \param startingTime initial time of the trajectory.
 */
void returnTrajectory( Eigen::VectorXd initialCartesianState,
                       double centralBodyGravitationalParameter,
                       double duration,
                       double maximumTimeStep,
                       std::vector < Eigen::Vector3d >& positionVector,
                       std::vector < double >& timeVector,
                       double startingTime = 0. );

//! Return a vector of positions and times corresponding to a 2D circular trajectory.
/*!
 * Return a vector of positions and times corresponding to a 2D circular trajectory.
 * \param maximumTimeStep maximum time step between the positions that are returned.
 * \param positionVector vector of positions along the trajectory.
 * \param timeVector vector of times corresponding to the positions.
 * \param semiMajorAxis semi-major axis of the circular orbit.
 * \param centralGravitationalParamater central body gravitational paremeter.
 * \param startingTime initial time of the circular trajectory.
 */
void returnCircularTrajectory( double maximumTimeStep,
                               std::vector < Eigen::Vector3d >& positionVector,
                               std::vector < double >& timeVector,
                               double semiMajorAxis,
                               double centralGravitationalParamater,
                               double startingTime = 0. );

//! Write a trajectory to a data file.
/*!
 * Write a trajectory to a data file.
 * \param positionVector vector of positions along the trajectory.
 * \param timeVector vector of times corresponding to the positions.
 * \param fileName name of the exported file
 */
void writeTrajectoryToFile(std::vector < Eigen::Vector3d > positionVector,
                            std::vector < double > timeVector,
                            const std::string fileName );

//! Write a trajectory to a data file.
/*!
 * Write a trajectory to a data file.
 * \param positionVector vector of positions of the manuevers
 * \param timeVector vector of times corresponding to the positions.
 * \param deltaVVector vector of delta V's of the maneuvers.
 * \param fileName name of the exported file.
 */
void writeManeuversToFile( std::vector < Eigen::Vector3d > positionVector,
                           std::vector < double > timeVector,
                           std::vector < double > deltaVVector,
                           const char * fileName );


} // namespace transfer_trajectories
} // namespace tudat

