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

#ifndef TUDAT_DEPARTURE_LEG_H
#define TUDAT_DEPARTURE_LEG_H

#include "spaceLeg.h"
#include "missionLeg.h"

namespace tudat
{
namespace spaceTrajectories
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
    DepartureLeg( const Eigen::Vector3d& departureBodyPosition,
                  const Eigen::Vector3d& arrivalBodyPosition,
                  const double timeOfFlight,
                  const Eigen::Vector3d& departureBodyVelocity,
                  const double centralBodyGravitationalParameter,
                  const double departureBodyGravitationalParameter,
                  const double semiMajorAxis,
                  const double eccentricity):
                SpaceLeg( departureBodyPosition,
                          arrivalBodyPosition,
                          timeOfFlight,
                          departureBodyVelocity,
                          centralBodyGravitationalParameter),
                departureBodyGravitationalParameter_( departureBodyGravitationalParameter ),
                semiMajorAxis_( semiMajorAxis ),
                eccentricity_( eccentricity ){ }
protected:

    //! The departure body gravitational parameter.
    /*!
     * The gravitational parameter of the departure body in the leg.
     */
    double departureBodyGravitationalParameter_;

    //! The semi major axis of the departure orbit.
    /*!
     * The semi major axis of the departure orbit, where the escape maneuver is applied.
     */
    double semiMajorAxis_;

    //! The eccentricity of the departure orbit.
    /*!
     * The eccentricity of the departure orbit, where the escape maneuver is applied.
     */
    double eccentricity_;

private:

};
} // namespace spaceTrajectories
} // namespace tudat

#endif // TUDAT_DEPARTURE_LEG_H
