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

#ifndef TUDAT_SWINGBY_LEG_MGA_1DSM_VELOCITY_H
#define TUDAT_SWINGBY_LEG_MGA_1DSM_VELOCITY_H

#include <vector>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLeg.h"


namespace tudat
{
namespace transfer_trajectories
{

//! Swingby Leg class of an MGA-1DSM velocity formulation trajectory model.
/*!
 * A class that calculates the required impulses for a swinbgy leg for an MGA-1DSM velocity
 * formulation trajectory model. This model specifies the trajectory by defining the initial
 * hyperbolic excess velocity (in 3D) and the moment of application of the DSM. The leg is split
 * in two sublegs. The first subleg is analyzed using a module that propagates the trajectory
 * using the given gravity assist parameters and subsequently by a kepler propagator. The second
 * leg is analyzed using a lambert targeter. Typically the user will be primarily interested in the
 * total delta V and the arrival velocity of the spacecraft before the arrival planet.
 * Note that often a MGA-1DSM velocity formulation trajectory model is used in which the gravity
 * assist is assumed to be unpowered. Both the powered and unpowered versions can be simulated
 * using this module.
 */
class SwingbyLegMga1DsmVelocity : public SwingbyLeg
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
     *  \param dsmTimeOfFlightFraction the fraction of the TOF at which the DSM is performed.
     *  \param rotationAngle The rotation angle of the gravity assist.
     *  \param pericenterRadius the minimum pericenter radius for the swing-by body.
     *  \param swingbyDeltaV the delta V for the swing-by manuever.
    */
    SwingbyLegMga1DsmVelocity( const Eigen::Vector3d& departureBodyPosition,
                               const Eigen::Vector3d& arrivalBodyPosition,
                               const double timeOfFlight,
                               const Eigen::Vector3d& departureBodyVelocity,
                               const double centralBodyGravitationalParameter,
                               const double swingbyBodyGravitationalParameter,
                               std::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr,
                               const double dsmTimeOfFlightFraction,
                               const double rotationAngle,
                               const double pericenterRadius,
                               const double swingbyDeltaV ):
        SwingbyLeg( departureBodyPosition,
                    arrivalBodyPosition,
                    timeOfFlight,
                    departureBodyVelocity,
                    centralBodyGravitationalParameter,
                    swingbyBodyGravitationalParameter,
                    velocityBeforeDepartureBodyPtr),
        dsmTimeOfFlightFraction_( dsmTimeOfFlightFraction ),
        rotationAngle_( rotationAngle ),
        pericenterRadius_( pericenterRadius ),
        swingbyDeltaV_( swingbyDeltaV )
    {
        velocityAfterDeparture_( 0 ) = TUDAT_NAN;
    }

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
     * the class, without re-initializing it. For this leg: time of flight, DSM time of flight
     * fraction, rotation angle, pericenter radius and the swing-by deltaV.
     * \param variableVector the new variable vector.
     */
    void updateDefiningVariables( const Eigen::VectorXd& variableVector );

protected:

private:
    //! The fraction of the time of flight of the DSM
    /*!
     * The fraction of the time of flight of the corresponding leg at which the DSM is performed.
     */
    double dsmTimeOfFlightFraction_;

    //! The rotation angle
    /*!
     * The rotation angle of the gravity assist. This angle defines the 3D orientation of the
     * gravity assist.
     */
    double rotationAngle_;

    //! The pericenter radius
    /*!
     * The pericenter radius of the gravity assist
     */
    double pericenterRadius_;

    //! The swing-by deltaV
    /*!
     * The swing-by deltaV that is used at pericenter of the gravity assist. (The Velocity effect
     * deltaV)
     */
    double swingbyDeltaV_;

    //! The DSM location
    /*!
     * The position at which the deep space maneuver is performed
     */
    Eigen::Vector3d dsmLocation_;

    //! The DSM time
    /*!
     * The time at which the deep space maneuver is performed
     */
    double dsmTime_;

    //! The velocity before the DSM
    /*!
     * The velocity of the spacecraft just before the deep space maneuver is performed
     */
    Eigen::Vector3d velocityBeforeDsm_;

    //! The velocity after the DSM
    /*!
     * The velocity of the spacecraft just after the deep space maneuver is performed
     */
    Eigen::Vector3d velocityAfterDsm_;

    //! The deltaV of the DSM
    /*!
     * The deltaV of the deep Space Maneuver
     */
    double deltaVDsm_;

};
} // namespace transfer_trajectories
} // namespace tudat

#endif // TUDAT_SWINGBY_LEG_MGA_1DSM_VELOCITY_H
