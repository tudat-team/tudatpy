/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120606    P. Musegaas       First creation of code.
 *
 *    References
 *
 */

// This contains some functions to facilitate making plots of trajectories.

#include <vector>

#include <Eigen/Core>

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace spaceTrajectories
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
void writeTrajectoryToFile( std::vector < Eigen::Vector3d > positionVector,
                            std::vector < double > timeVector,
                            const char * fileName );

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


} // namespace spaceTrajectories
} // namespace tudat

