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
 *      120611    P. Musegaas       First creation of code.
 *
 *    References
 *
 */

#ifndef TUDAT_CAPTURE_LEG_H
#define TUDAT_CAPTURE_LEG_H

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "missionLeg.h"

namespace tudat
{
namespace spaceTrajectories
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

    //! Constructor with immediate definition of parameters.
    /*!
     * Constructor with immediate definition of parameters.
     */
    CaptureLeg( const Eigen::Vector3d& departureBodyPosition,
                const double timeOfFlight,
                const Eigen::Vector3d& departureBodyVelocity,
                const double centralBodyGravitationalParameter,
                const double captureBodyGravitationalParameter,
                boost::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr,
                const double semiMajorAxis,
                const double eccentricity ):
            MissionLeg( departureBodyPosition,
                        timeOfFlight,
                        departureBodyVelocity,
                        centralBodyGravitationalParameter),
            captureBodyGravitationalParameter_( captureBodyGravitationalParameter ),
            velocityBeforeDepartureBodyPtr_( velocityBeforeDepartureBodyPtr ),
            semiMajorAxis_( semiMajorAxis ),
            eccentricity_( eccentricity )

    {
        velocityAfterDeparture_( 0 ) = TUDAT_NAN;
    }

    //! Calculate the leg
    /*!
     * Performs all calculations required for this leg.
     */
    void calculateLeg( Eigen::Vector3d& velocityBeforeArrivalBody,
                       double& deltaV );

    //! Calculate intermediate positions and their corresponding times.
    /*!
     * Calculates intermediate positions and their corresponding times in the leg, based on a
     * maximum time between two points.
     */
    void intermediatePoints( const double maximumTimeStep,
                             std::vector < Eigen::Vector3d >& positionVector,
                             std::vector < double >& timeVector,
                             const double startingTime = 0. );

    //! Return maneuvres along the leg.
    /*!
     * Returns the maneuver points, times and sizes along the trajectory.
     */
    void maneuvers( std::vector < Eigen::Vector3d >& positionVector,
                    std::vector < double >& timeVector,
                    std::vector < double >& deltaVVector,
                    const double startingTime = 0. );

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
    }

    //! Update the defining variables.
    /*!
     * Sets the trajectory defining variables to the newly specified values. Required for re-using
     * the class, without re-initializing it. For this leg: time of flight.
     */
    void updateDefiningVariables( const Eigen::VectorXd& variableVector );

protected:

private:

    //! The capture body gravitational parameter.
    /*!
     * The gravitational parameter of the capture body in the leg.
     */
    double captureBodyGravitationalParameter_;

    //! The velocity of the spacecraft before capture.
    /*!
     * The heliocentric velocity before the capture maneuver. This is passed using a pointer,
     * because this parameter is dependent on the previous leg. The result of the previous leg is
     * not necessarily known when this leg is initiated.
     */
    boost::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr_;

    //! The capture orbit semi major axis.
    /*!
     * The semi major axis of the intended capture orbit.
     */
    double semiMajorAxis_;

    //! The capture orbit eccentricity.
    /*!
     * The eccentricity of the intended capture orbit.
     */
    double eccentricity_;

};
} // namespace spaceTrajectories
} // namespace tudat

#endif // TUDAT_CAPTURE_LEG_H
