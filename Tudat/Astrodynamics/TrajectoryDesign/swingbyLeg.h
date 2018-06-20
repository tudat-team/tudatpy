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

#ifndef TUDAT_SWINGBY_LEG_H
#define TUDAT_SWINGBY_LEG_H

#include <boost/make_shared.hpp>
#include <memory>

#include "Tudat/Astrodynamics/TrajectoryDesign/spaceLeg.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Swingby Leg base class.
/*!
 * A base class that calculates the required impulses for the second and subsequent interplanetary
 * legs in a hight-thrust, patched-conics space trajectory.
 * It inherits from the interplanetary leg base class. It is inherited by different trajectory
 * models. In the base class the basic functionalities that are foreseen for the different
 * trajectory models are declared and implemented.
 * Apart from those declared in the space leg base class, these include: the velocity before the
 * departure planet and the gravitational parameter of the swing-by body.
 */
class SwingbyLeg : public SpaceLeg
{

public:

protected:
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
    */
    SwingbyLeg( const Eigen::Vector3d& departureBodyPosition,
                const Eigen::Vector3d& arrivalBodyPosition,
                const double timeOfFlight,
                const Eigen::Vector3d& departureBodyVelocity,
                const double centralBodyGravitationalParameter,
                double swingbyBodyGravitationalParameter,
                std::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr ):
        SpaceLeg( departureBodyPosition,
                  arrivalBodyPosition,
                  timeOfFlight,
                  departureBodyVelocity,
                  centralBodyGravitationalParameter),
        swingbyBodyGravitationalParameter_( swingbyBodyGravitationalParameter ),
        velocityBeforeDepartureBodyPtr_( velocityBeforeDepartureBodyPtr ){ }

    virtual ~SwingbyLeg( ){ }

    //! The swing-by body gravitational parameter.
    /*!
     * The gravitational parameter of the swing-by body in the leg.
     */
    double swingbyBodyGravitationalParameter_;


    //! The velocity of the spacecraft at departure.
    /*!
     * The heliocentric velocity at the departure time (before gravity assist). This is passed
     * using a pointer, because this parameter is dependent on the previous leg. The result of
     * the previous leg is not necessarily known when this leg is initiated.
     */
    std::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr_;


private:

};
} // namespace transfer_trajectories
} // namespace tudat

#endif // TUDAT_SWINGBY_LEG_H
