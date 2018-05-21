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
 *      120509    P. Musegaas       First creation of code.
 *      120611    P. Musegaas       Adaptation to new mission segments functions and update of
 *                                  of functionality.
 *
 *    References
 *
 */

#ifndef TUDAT_SWINGBY_LEG_H
#define TUDAT_SWINGBY_LEG_H

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "spaceLeg.h"

namespace tudat
{
namespace spaceTrajectories
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

    SwingbyLeg( const Eigen::Vector3d& departureBodyPosition,
                const Eigen::Vector3d& arrivalBodyPosition,
                const double timeOfFlight,
                const Eigen::Vector3d& departureBodyVelocity,
                const double centralBodyGravitationalParameter,
                double swingbyBodyGravitationalParameter,
                boost::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr ):
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
    boost::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr_;


private:

};
} // namespace spaceTrajectories
} // namespace tudat

#endif // TUDAT_SWINGBY_LEG_H
