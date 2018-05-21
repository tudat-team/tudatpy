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

#ifndef TUDAT_DEPARTURE_LEG_MGA_1DSM_VELOCITY_H
#define TUDAT_DEPARTURE_LEG_MGA_1DSM_VELOCITY_H

#include <vector>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "departureLeg.h"

namespace tudat
{
namespace spaceTrajectories
{

//! Departure Leg class of an MGA-1DSM velocity formulation trajectory model.
/*!
 * A class that calculates the required impulses for a departure leg for an MGA-1DSM velocity
 * formulation trajectory model.
 * This model specifies the trajectory by defining the initial hyperbolic excess velocity (in 3D)
 * and the moment of application of the DSM. The leg is split in two sublegs. The first subleg is
 * analyzed using a kepler propagator, whereas the second leg is analyzed using a lambert targeter.
 * Typically the user will be primarily interested in the total delta V and the arrival velocity of
 * the spacecraft before the arrival planet.
 */
class DepartureLegMga1DsmVelocity : public DepartureLeg
{
public:

    //! Constructor with immediate definition of parameters.
    /*!
     * Constructor with immediate definition of parameters.
     */
    DepartureLegMga1DsmVelocity( const Eigen::Vector3d& departureBodyPosition,
                                 const Eigen::Vector3d& arrivalBodyPosition,
                                 const double timeOfFlight,
                                 const Eigen::Vector3d& departureBodyVelocity,
                                 const double centralBodyGravitationalParameter,
                                 const double departureBodyGravitationalParameter,
                                 const double semiMajorAxis,
                                 const double eccentricity,
                                 const double dsmTimeOfFlightFraction,
                                 const double excessVelocityMagnitude,
                                 const double excessVelocityInPlaneAngle,
                                 const double excessVelocityOutOfPlaneAngle ):
        DepartureLeg( departureBodyPosition,
                      arrivalBodyPosition,
                      timeOfFlight,
                      departureBodyVelocity,
                      centralBodyGravitationalParameter,
                      departureBodyGravitationalParameter,
                      semiMajorAxis,
                      eccentricity),
        dsmTimeOfFlightFraction_( dsmTimeOfFlightFraction ),
        excessVelocityMagnitude_( excessVelocityMagnitude ),
        excessVelocityInPlaneAngle_( excessVelocityInPlaneAngle ),
        excessVelocityOutOfPlaneAngle_( excessVelocityOutOfPlaneAngle )
    {
        velocityAfterDeparture_( 0 ) = TUDAT_NAN;
    }

    //! Calculates the leg.
    /*!
     * Performs all calculations required for this leg by the associated trajectory model.
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

    //! Update the defining variables.
    /*!
     * Sets the trajectory defining variables to the newly specified values. Required for re-using
     * the class, without re-initializing it. For this leg: time of flight, DSM time of flight
     * fraction, excess velocity magnitude, excess velocity in plane angle and the excess velocity
     * out of plane angle.
     */
    void updateDefiningVariables( const Eigen::VectorXd& variableVector );

protected:

private:
    //! The fraction of the time of flight of the DSM.
    /*!
     * The fraction of the time of flight of the corresponding leg at which the DSM is performed.
     */
    double dsmTimeOfFlightFraction_;

    //! The excess velocity magnitude.
    /*!
     * The magnitude of the relative velocity of the spacecraft leaving the departure planet with
     * the departure planet.
     */
    double excessVelocityMagnitude_;

    //! The excess velocity in plane angle.
    /*!
     * The rotation angle of the excess velocity in the xy-plane of the central reference frame.
     */
    double excessVelocityInPlaneAngle_;

    //! The excess velocity out of plane angle.
    /*!
     * The rotation angle of the excess velocity in the z-direction of the central reference frame.
     */
    double excessVelocityOutOfPlaneAngle_;

    //! The DSM location.
    /*!
     * The position at which the deep space maneuver is performed.
     */
    Eigen::Vector3d dsmLocation_;

    //! The DSM time.
    /*!
     * The time at which the deep space maneuver is performed.
     */
    double dsmTime_;

    //! The velocity before the DSM.
    /*!
     * The velocity of the spacecraft just before the deep space maneuver is performed.
     */
    Eigen::Vector3d velocityBeforeDsm_;

    //! The velocity after the DSM.
    /*!
     * The velocity of the spacecraft just after the deep space maneuver is performed.
     */
    Eigen::Vector3d velocityAfterDsm_;

    //! The deltaV of the DSM.
    /*!
     * The deltaV of the deep Space Maneuver.
     */
    double deltaVDsm_;

    //! The deltaV of the departure maneuver.
    /*!
     * The deltaV of the departure maneuver.
     */
    double deltaVDeparture_;

};

} // namespace spaceTrajectories
} // namespace tudat

#endif // TUDAT_DEPARTURE_LEG_MGA_1DSM_VELOCITY_H
