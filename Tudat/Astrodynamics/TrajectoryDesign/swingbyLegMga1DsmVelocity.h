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

#ifndef TUDAT_SWINGBY_LEG_MGA_1DSM_VELOCITY_H
#define TUDAT_SWINGBY_LEG_MGA_1DSM_VELOCITY_H

#include <vector>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "swingbyLeg.h"


namespace tudat
{
namespace spaceTrajectories
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
     * Constructor with immediate definition of parameters.
     */
    SwingbyLegMga1DsmVelocity( const Eigen::Vector3d& departureBodyPosition,
                               const Eigen::Vector3d& arrivalBodyPosition,
                               const double timeOfFlight,
                               const Eigen::Vector3d& departureBodyVelocity,
                               const double centralBodyGravitationalParameter,
                               const double swingbyBodyGravitationalParameter,
                               boost::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr,
                               const double dsmTimeOfFlightFraction,
                               const double rotationAngle,
                               const double pericenterRadius,
                               const double swingbyDeltaV
                               ):
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

    //! Calculates the leg
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
     * fraction, rotation angle, pericenter radius and the swing-by deltaV.
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
} // namespace spaceTrajectories
} // namespace tudat

#endif // TUDAT_SWINGBY_LEG_MGA_1DSM_VELOCITY_H
