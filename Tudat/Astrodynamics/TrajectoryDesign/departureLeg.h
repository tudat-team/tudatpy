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

#ifndef TUDAT_DEPARTURE_LEG_H
#define TUDAT_DEPARTURE_LEG_H

#include "Tudat/Astrodynamics/TrajectoryDesign/spaceLeg.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/missionLeg.h"

namespace tudat
{

namespace transfer_trajectories
{

//! Departure Leg base class.
/*!
 * A base class that calculates the required impulses for a high-thrust, patched-conics departure
 * leg.
 * It inherits from the space leg base class. It is inherited by different trajectory models.
 * The basic variables and functions are similar to the interplanetary leg base class and hence
 * this is an empty base class. It is however included to point out the distinction with a swing-by
 * leg, which does declare new variables with respect to the interplantary leg base class.
 */
class DepartureLeg : public SpaceLeg
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
     *  \param departureBodyGravitationalParameter gravitational parameter of the departure body.
     *  \param semiMajorAxis semi-major axis of the orbit after the capture is performed.
     *  \param eccentricity eccentricity of the orbit after the capture is performed.
     *  \param includeDepartureDeltaV Boolean denoting whether to include the Delta V of departure.
     */
    DepartureLeg( const Eigen::Vector3d& departureBodyPosition,
                  const Eigen::Vector3d& arrivalBodyPosition,
                  const double timeOfFlight,
                  const Eigen::Vector3d& departureBodyVelocity,
                  const double centralBodyGravitationalParameter,
                  const double departureBodyGravitationalParameter,
                  const double semiMajorAxis,
                  const double eccentricity,
                  const bool includeDepartureDeltaV = true ):
        SpaceLeg( departureBodyPosition,
                  arrivalBodyPosition,
                  timeOfFlight,
                  departureBodyVelocity,
                  centralBodyGravitationalParameter),
        departureBodyGravitationalParameter_( departureBodyGravitationalParameter ),
        semiMajorAxis_( semiMajorAxis ),
        eccentricity_( eccentricity ),
        includeDepartureDeltaV_( includeDepartureDeltaV ){ }

    //! Function to retrieve the value of the escape Delta V.
    /*!
     *  Function to retrieve the value of the escape Delta V.
     *  \param escapeDeltaV Double denoting the value of the escape Delta V.
     *  \return Double denoting the value of the escape Delta V (returned by reference ).
     */
    void getEscapeDeltaV( double& escapeDeltaV )
    {
        escapeDeltaV = escapeDeltaV_;
    }

protected:

    //! The departure body gravitational parameter.
    /*!
     *  The gravitational parameter of the departure body in the leg.
     */
    double departureBodyGravitationalParameter_;

    //! The semi-major axis of the departure orbit.
    /*!
     * The semi-major axis of the departure orbit, where the escape maneuver is applied.
     */
    double semiMajorAxis_;

    //! The eccentricity of the departure orbit.
    /*!
     *  The eccentricity of the departure orbit, where the escape maneuver is applied.
     */
    double eccentricity_;

    //! Boolean denoting whether to include the Delta V of arrival.
    bool includeDepartureDeltaV_;

    //! Double denoting the value of the escape Delta V.
    double escapeDeltaV_;

private:

};

} // namespace transfer_trajectories

} // namespace tudat

#endif // TUDAT_DEPARTURE_LEG_H
