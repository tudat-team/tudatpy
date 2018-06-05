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

#ifndef TUDAT_MISSION_LEG_H
#define TUDAT_MISSION_LEG_H

#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace transfer_trajectories
{

//! Mission Leg base class.
/*!
 * Abstract base class for a mission leg.
 * This base class is inherited by two types of legs: a space leg, which is a transfer from one
 * body to another, and a capture leg, in which the spacecraft remains at a certain body for a
 * specified time in a specified orbit.
 * This base class contains the basic methods required in a trajectory class. These are a calculate
 * leg function, an intermediate points function for plotting purposes and finally a maneuvers
 * method to extract information about the maneuvers in the leg.
 * Also the variables required for all variables are declared here.
 */
class MissionLeg
{
public:
    //! Constructor with immediate definition of parameters.
    /*!
     *  Constructor, sets objects and functions from which relevant environment and state variables
     *  are retrieved.
     *  \param departureBodyPosition location of the departure body.
     *  \param timeOfFlight Length of the leg.
     *  \param departureBodyVelocity velocity of the departure body.
     *  \param centralBodyGravitationalParameter gravitational parameter of the cebtral body (most cases the Sun).
     */
    MissionLeg(const Eigen::Vector3d& departureBodyPosition,
               const double timeOfFlight,
               const Eigen::Vector3d& departureBodyVelocity,
               const double centralBodyGravitationalParameter ):
            departureBodyPosition_( departureBodyPosition ),
            timeOfFlight_( timeOfFlight ),
            departureBodyVelocity_( departureBodyVelocity ),
            centralBodyGravitationalParameter_( centralBodyGravitationalParameter ){ }

    //! virtual destructor.
    /*!
     * virtual destructor.
     */
          virtual ~MissionLeg( ){ }

    //! Calculate the leg
    /*!
     * Performs all calculations required for this leg, in this class it is pure virtual.
     *  \param velocityBeforeArrivalBody the velocity of the spacecraft before it arrives at the target body.
     *  \param deltaV the delta V required to perform the leg.
     */
    virtual void calculateLeg( Eigen::Vector3d& velocityBeforeArrivalBody,
                               double& deltaV ) = 0;

    //! Calculate intermediate positions and their corresponding times.
    /*!
     * Calculates intermediate positions and their corresponding times in the leg, based on a
     * maximum time between two points.
     *  \param maximumTimeStep the maximum time between two points along the trajectory.
     *  \param positionVector Vector of positions along the orbit, space according to the maximum time step.
     *  \param timeVector The times corresponding to the positions.
     *  \param startingTime the initial time from which the intermediate points are given.
     */
    virtual void intermediatePoints( const double maximumTimeStep,
                                     std::vector < Eigen::Vector3d >& positionVector,
                                     std::vector < double >& timeVector,
                                     const double startingTime = 0. ) = 0;

    //! Return maneuvres along the leg.
    /*!
     * Returns the maneuver points, times and sizes along the trajectory.
     *  \param positionVector Vector of the positions of the maneuvers.
     *  \param timeVector The times corresponding to the positions.
     *  \param deltaVVector the delta V required for each maneuver.
     *  \param startingTime the initial time from which the maneuvers are given.
     */
    virtual void maneuvers( std::vector < Eigen::Vector3d >& positionVector,
                            std::vector < double >& timeVector,
                            std::vector < double >& deltaVVector,
                            const double startingTime = 0. ) = 0;

    //! Update the ephemeris.
    /*!
     * Sets the positions and the velocities to the newly specified values. Required for re-using
     * the class, without re-initializing it.
     *  \param departureBodyPosition sets the new departure body position.
     *  \param arrivalBodyPosition sets the new arrival body position.
     *  \param departureBodyVelocity sets the new departure body velocity.
     */
    virtual void updateEphemeris( const Eigen::Vector3d& departureBodyPosition,
                                  const Eigen::Vector3d& arrivalBodyPosition,
                                  const Eigen::Vector3d& departureBodyVelocity ) = 0;

    //! Update the defining variables.
    /*!
     * Sets the trajectory defining variables to the newly specified values. Required for re-using
     * the class, without re-initializing it.
     * \param variableVector the new variable vector.
     */
    virtual void updateDefiningVariables( const Eigen::VectorXd& variableVector ) = 0;

    //! Return Departure Variables.
    /*!
     * Returns the departure body position, velocity and the velocity after departure. Mainly used
     * for getting the launch conditions for TandEM like problems.
     *  \param departureBodyPosition the departure body position.
     *  \param departureBodyVelocity the velocity of the departure body.
     *  \param velocityAfterDeparture the velocity of the spacecraft after the departure maneuver.
     */
    void returnDepartureVariables( Eigen::Vector3d& departureBodyPosition,
                                   Eigen::Vector3d& departureBodyVelocity,
                                   Eigen::Vector3d& velocityAfterDeparture )
    {
        departureBodyPosition = departureBodyPosition_;
        departureBodyVelocity = departureBodyVelocity_;
        velocityAfterDeparture = velocityAfterDeparture_;
    }

protected:

    //! The departure body position.
    /*!
     * The position of the departure body at the departure time.
     */
    Eigen::Vector3d departureBodyPosition_;

    //! The time of flight.
    /*!
     * The time of flight for this leg.
     */
    double timeOfFlight_;

    //! The departure body velocity.
    /*!
     * The velocity of the departure body at the departure time.
     */
    Eigen::Vector3d departureBodyVelocity_;

    //! The central body gravitational parameter.
    /*!
     * The gravitational parameter of the central body in the leg.
     */
    double centralBodyGravitationalParameter_;

    //! The velocity after departure.
    /*!
     * The heliocentric velocity of the spacecraft right after the departure maneuver.
     */
    Eigen::Vector3d velocityAfterDeparture_;

    //! The total deltaV of the leg.
    /*!
     * The total deltaV that is required for this leg.
     */
    double deltaV_;

private:


};
} // namespace transfer_trajectories
} // namespace tudat

#endif // TUDAT_MISSION_LEG_H
