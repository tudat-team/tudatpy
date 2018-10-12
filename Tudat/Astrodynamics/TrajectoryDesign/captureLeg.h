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

#ifndef TUDAT_CAPTURE_LEG_H
#define TUDAT_CAPTURE_LEG_H

#include <boost/make_shared.hpp>
#include <memory>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/TrajectoryDesign/missionLeg.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Capture leg class.
/*!
 * A class that calculates capture maneuver of a spacecraft entering a specified orbit around a
 * body. It also contains the same functionality as the space legs, meaning this class can also
 * be used for a captured phase in a trajectory. For example to calculate a sample & return mission
 * in one go.
 */
class CaptureLeg : public MissionLeg
{
public:

    //! Constructor with immediate definition of parameters..
    /*!
     *  Constructor, sets objects and functions from which relevant environment and state variables
     *  are retrieved.
     *  \param departureBodyPosition Location of the departure body.
     *  \param timeOfFlight Length of the leg.
     *  \param departureBodyVelocity Velocity of the departure body.
     *  \param centralBodyGravitationalParameter Gravitational parameter of the cebtral body (most cases the Sun).
     *  \param captureBodyGravitationalParameter Gravitational parameter of the capture body.
     *  \param velocityBeforeDepartureBodyPtr Pointer to the velocity before arriving at the departure body.
     *  \param semiMajorAxis Semi-major axis of the orbit after the capture is performed.
     *  \param eccentricity Eccentricity of the orbit after the capture is performed.
     *  \param includeArrivalDeltaV Boolean denoting whether to include the Delta V of arrival.
     */
    CaptureLeg( const Eigen::Vector3d& departureBodyPosition,
                const double timeOfFlight,
                const Eigen::Vector3d& departureBodyVelocity,
                const double centralBodyGravitationalParameter,
                const double captureBodyGravitationalParameter,
                std::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr,
                const double semiMajorAxis,
                const double eccentricity,
                const bool includeArrivalDeltaV = true ):
            MissionLeg( departureBodyPosition,
                        timeOfFlight,
                        departureBodyVelocity,
                        centralBodyGravitationalParameter),
            captureBodyGravitationalParameter_( captureBodyGravitationalParameter ),
            velocityBeforeDepartureBodyPtr_( velocityBeforeDepartureBodyPtr ),
            semiMajorAxis_( semiMajorAxis ),
            eccentricity_( eccentricity ),
            includeArrivalDeltaV_( includeArrivalDeltaV )
    {
        velocityAfterDeparture_( 0 ) = TUDAT_NAN;
    }

    //! Calculate the leg
    /*!
     * Performs all calculations required for this leg.
     *  \param velocityBeforeArrivalBody the velocity of the spacecraft before it arrives at the target body.
     *  \param deltaV the delta V required to perform the capture manuever.
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
    }

    //! Update the defining variables.
    /*!
     *  Sets the trajectory defining variables to the newly specified values. Required for re-using
     *  the class, without re-initializing it. For this leg: time of flight.
     *  \param variableVector the new variable vector.
     */
    void updateDefiningVariables( const Eigen::VectorXd& variableVector );

    //! Function to retrieve the value of the capture Delta V.
    /*!
     *  Function to retrieve the value of the capture Delta V.
     *  \param captureDeltaV Double denoting the value of the capture Delta V.
     *  \return Double denoting the value of the capture Delta V (returned by reference ).
     */
    void getCaptureDeltaV( double& captureDeltaV )
    {
        captureDeltaV = deltaV_;
    }

protected:

private:

    //! The capture body gravitational parameter.
    /*!
     *  The gravitational parameter of the capture body in the leg.
     */
    double captureBodyGravitationalParameter_;

    //! The velocity of the spacecraft before capture.
    /*!
     *  The heliocentric velocity before the capture maneuver. This is passed using a pointer,
     *  because this parameter is dependent on the previous leg. The result of the previous leg is
     *  not necessarily known when this leg is initiated.
     */
    std::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr_;

    //! The capture orbit semi-major axis.
    /*!
     *  The semi-major axis of the intended capture orbit.
     */
    double semiMajorAxis_;

    //! The capture orbit eccentricity.
    /*!
     * The eccentricity of the intended capture orbit.
     */
    double eccentricity_;

    //! Boolean denoting whether to include the Delta V of arrival.
    bool includeArrivalDeltaV_;

};

} // namespace transfer_trajectories

} // namespace tudat

#endif // TUDAT_CAPTURE_LEG_H
