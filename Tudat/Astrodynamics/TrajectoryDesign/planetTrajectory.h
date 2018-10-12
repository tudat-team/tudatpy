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

// This contains some functions to facilitate making plots of planet trajectories.

#include <vector>

#include <Eigen/Core>

#include <memory>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Return a vector of positions and times from ephemeris data for a certain epoch and duration.
/*!
 * Return a vector of positions and times from ephemeris data for a certain epoch and duration.
 * \param ephemerisPtr pointer to the ephemeris of the planet.
 * \param duration duration for which the planet trajectory is given.
 * \param maximumTimeStep maximum time between points along the trajectory.
 * \param positionVector vector of positions along the trajectory.
 * \param timeVector times corresponding to the positions.
 * \param startingTime initial time of the trajectory.
 */
void returnPlanetTrajectory(const ephemerides::EphemerisPointer& ephemerisPtr,
                             const double duration,
                             const double maximumTimeStep,
                             std::vector < Eigen::Vector3d >& positionVector,
                             std::vector < double >& timeVector,
                             const double startingTime );

//! Return a propagated vector of positions and times from ephemeris data for one revolution.
/*!
 * Return a propagated vector of positions and times from ephemeris data for one revolution.
 * \param ephemerisPtr pointer to the ephemeris of the planet.
 * \param centralBodyGravitationalParameter central body gravitational parameter.
 * \param startingEpochMJD2000 starting epoch in MJD2000.
 * \param maximumTimeStep maximum time between points along the trajectory.
 * \param positionVector vector of positions along the trajectory.
 * \param timeVector times corresponding to the positions.
 * \param startingTime initial time of the trajectory.
 */
void returnSingleRevolutionPlanetTrajectory(
        const ephemerides::EphemerisPointer& ephemerisPtr,
        const double centralBodyGravitationalParameter,
        const double startingEpochMJD2000,
        const double maximumTimeStep,
        std::vector < Eigen::Vector3d >& positionVector,
        std::vector < double >& timeVector,
        const double startingTime = 0. );

} // namespace transfer_trajectories
} // namespace tudat

