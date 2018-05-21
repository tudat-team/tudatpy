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

#ifndef TUDAT_SPACE_LEG_H
#define TUDAT_SPACE_LEG_H

#include <vector>

#include <Eigen/Core>

#include "missionLeg.h"

namespace tudat
{
namespace spaceTrajectories
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
} // namespace spaceTrajectories
} // namespace tudat

#endif // TUDAT_SPACE_LEG_H
