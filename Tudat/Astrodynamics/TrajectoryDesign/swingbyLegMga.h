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

#ifndef TUDAT_SWINGBY_LEG_MGA_H
#define TUDAT_SWINGBY_LEG_MGA_H

#include <vector>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLeg.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Swingby Leg class of an MGA trajectory model.
/*!
 * A class that calculates the required impulses for a swingby leg for an MGA trajectory model.
 * It extends the swingby leg base class and uses a lambert targeter for finding the required
 * departure and arrival velocities of the spacecraft.
 */
class SwingbyLegMga : public SwingbyLeg
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
     *  \param swingbyBodyGravitationalParameter gravitational parameter of the swing-by body.
     *  \param velocityBeforeDepartureBodyPtr pointer to the velocity before the swing-by.
     *  \param minimumPericenterRadius the minimum pericenter radius for the swing-by body.
    */
    SwingbyLegMga( const Eigen::Vector3d& departureBodyPosition,
                   const Eigen::Vector3d& arrivalBodyPosition,
                   const double timeOfFlight,
                   const Eigen::Vector3d& departureBodyVelocity,
                   const double centralBodyGravitationalParameter,
                   double swingbyBodyGravitationalParameter,
                   std::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr,
                   const double minimumPericenterRadius
                   ):
        SwingbyLeg( departureBodyPosition,
                    arrivalBodyPosition,
                    timeOfFlight,
                    departureBodyVelocity,
                    centralBodyGravitationalParameter,
                    swingbyBodyGravitationalParameter,
                    velocityBeforeDepartureBodyPtr),
        minimumPericenterRadius_( minimumPericenterRadius )
    {
        velocityAfterDeparture_( 0 ) = TUDAT_NAN;
    }

    ~SwingbyLegMga(){}

    //! Calculate the leg
    /*!
     * Performs all calculations required for this leg.
     *  \param velocityBeforeArrivalBody the velocity of the spacecraft before it arrives at the target body.
     *  \param deltaV the delta V required to perform the leg.
     */
    void calculateLeg( Eigen::Vector3d& velocityBeforeArrivalBody,
                       double& deltaV );

    //! Calculate intermediate positions and their corresponding times.
    /*!
     * Calculates intermediate positions and their corresponding times in the leg, based on a
     * maximum time between two points.
     *  \param maximumTimeStep the maximum time between two points along the trajectory.
     *  \param positionVector Vector of positions along the orbit, space according to the maximum time step.
     *  \param timeVector The times corresponding to the positions.
     *  \param startingTime the initial time from which the intermediate points are given.
     */
    void intermediatePoints( const double maximumTimeStep,
                             std::vector < Eigen::Vector3d >& positionVector,
                             std::vector < double >& timeVector,
                             const double startingTime = 0. );

    //! Return maneuvres along the leg.
    /*!
     * Returns the maneuver points, times and sizes along the trajectory.
     *  \param positionVector Vector of the positions of the maneuvers.
     *  \param timeVector The times corresponding to the positions.
     *  \param deltaVVector the delta V required for each maneuver.
     *  \param startingTime the initial time from which the maneuvers are given.
     */
    void maneuvers( std::vector < Eigen::Vector3d >& positionVector,
                    std::vector < double >& timeVector,
                    std::vector < double >& deltaVVector,
                    const double startingTime = 0. );

    //! Update the defining variables.
    /*!
     * Sets the trajectory defining variables to the newly specified values. Required for re-using
     * the class, without re-initializing it. For this leg: time of flight.
     *  \param variableVector the new variable vector.
     */
    void updateDefiningVariables( const Eigen::VectorXd& variableVector );

protected:

private:

    //! The minimum pericenter radius
    /*!
     * The minimum allowable pericenter radius for the swing-by.
     */
    double minimumPericenterRadius_;

};
} // namespace transfer_trajectories
} // namespace tudat

#endif // TUDAT_SWINGBY_LEG_MGA_H
