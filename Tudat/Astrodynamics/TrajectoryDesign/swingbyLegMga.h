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

#ifndef TUDAT_SWINGBY_LEG_MGA_H
#define TUDAT_SWINGBY_LEG_MGA_H

#include <vector>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "swingbyLeg.h"

namespace tudat
{
namespace spaceTrajectories
{

//! Swingby Leg class of an MGA trajectory model.
/*!
 * A class that calculates the required impulses for a swingby leg for an MGA trajectory model.
 * It extends the swingby leg base class and uses a lambert targeter for finding the required
 * departure and arrival velocities of the spacecraft.
 */
class SwingbyLegMga : public SwingbyLeg
{
public:

    //! Constructor with immediate definition of parameters.
    /*!
     * Constructor with immediate definition of parameters.
     */
    SwingbyLegMga( const Eigen::Vector3d& departureBodyPosition,
                   const Eigen::Vector3d& arrivalBodyPosition,
                   const double timeOfFlight,
                   const Eigen::Vector3d& departureBodyVelocity,
                   const double centralBodyGravitationalParameter,
                   double swingbyBodyGravitationalParameter,
                   boost::shared_ptr< Eigen::Vector3d > velocityBeforeDepartureBodyPtr,
                   const double minimumPericenterRadius
                   ):
        SwingbyLeg( departureBodyPosition,
                    arrivalBodyPosition,
                    timeOfFlight,
                    departureBodyVelocity,
                    centralBodyGravitationalParameter,
                    swingbyBodyGravitationalParameter,
                    velocityBeforeDepartureBodyPtr),
        minimumPericenterRadius_( minimumPericenterRadius )
    {
        velocityAfterDeparture_( 0 ) = TUDAT_NAN;
    }

    ~SwingbyLegMga(){}

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
     * the class, without re-initializing it. For this leg: time of flight.
     */
    void updateDefiningVariables( const Eigen::VectorXd& variableVector );

protected:

private:

    //! The minimum pericenter radius
    /*!
     * The minimum allowable pericenter radius for the swing-by.
     */
    double minimumPericenterRadius_;

};
} // namespace spaceTrajectories
} // namespace tudat

#endif // TUDAT_SWINGBY_LEG_MGA_H
