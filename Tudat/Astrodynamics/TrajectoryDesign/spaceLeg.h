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

#ifndef TUDAT_SPACE_LEG_H
#define TUDAT_SPACE_LEG_H

#include <vector>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/TrajectoryDesign/missionLeg.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Space Leg base class.
/*!
 * A base class for calculating the required impulses for a high-thrust, patched-conics space leg.
 * It inherits from the mission leg base class. The space leg class is inherited by different
 * trajectory models, for departure legs and swing-by legs.
 * One additional variable is declared here.
 */
class SpaceLeg : public MissionLeg
{
public:
    //! Constructor with immediate definition of parameters.
    /*!
     *  Constructor, sets objects and functions from which relevant environment and state variables
     *  are retrieved.
     *  \param departureBodyPosition location of the departure body.
     *  \param arrivalBodyPosition position of the target body.
     *  \param timeOfFlight Length of the leg.
     *  \param departureBodyVelocity velocity of the departure body.
     *  \param centralBodyGravitationalParameter gravitational parameter of the cebtral body (most cases the Sun).
     */
    SpaceLeg( const Eigen::Vector3d& departureBodyPosition,
              const Eigen::Vector3d& arrivalBodyPosition,
              const double timeOfFlight,
              const Eigen::Vector3d& departureBodyVelocity,
              const double centralBodyGravitationalParameter ):
        MissionLeg( departureBodyPosition,
                    timeOfFlight,
                    departureBodyVelocity,
                    centralBodyGravitationalParameter),
        arrivalBodyPosition_( arrivalBodyPosition ){ }

    virtual ~SpaceLeg( ){ }

    //! Update the ephemeris.
    /*!
     * Sets the positions and the velocities to the newly specified values. Required for re-using
     * the class, without re-initializing it.
     *  \param departureBodyPosition sets the new departure body position.
     *  \param arrivalBodyPosition sets the new arrival body position.
     *  \param departureBodyVelocity sets the new departure body velocity.
     */
    void updateEphemeris( const Eigen::Vector3d& departureBodyPosition,
                          const Eigen::Vector3d& arrivalBodyPosition,
                          const Eigen::Vector3d& departureBodyVelocity )
    {
        departureBodyPosition_ = departureBodyPosition;
        departureBodyVelocity_ = departureBodyVelocity;
        arrivalBodyPosition_ = arrivalBodyPosition;
    }

protected:

    //! The arrival body position.
    /*!
     * The position of the arrival body at the arrival time.
     */
    Eigen::Vector3d arrivalBodyPosition_;

private:


};
} // namespace transfer_trajectories
} // namespace tudat

#endif // TUDAT_SPACE_LEG_H
