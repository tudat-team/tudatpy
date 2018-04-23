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

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/Astrodynamics/MissionSegments/gravityAssist.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"

#include "exportTrajectory.h"
#include "swingbyLegMga.h"

namespace tudat
{
namespace spaceTrajectories
{

//! Calculate the leg and update the Delta V and the velocity before the next body.
void SwingbyLegMga::calculateLeg( Eigen::Vector3d& velocityBeforeArrivalBody,
                                  double& deltaV )
{
    // Calculate and set the spacecraft velocities after departure and before arrival.
    mission_segments::solveLambertProblemIzzo( departureBodyPosition_, arrivalBodyPosition_,
                                               timeOfFlight_, centralBodyGravitationalParameter_,
                                               velocityAfterDeparture_, velocityBeforeArrivalBody );

    // Perform a gravity assist at the current body, using the previously obtained velocity and
    // the properties of the swing-by body. Store the deltaV required.
    deltaV_ = mission_segments::gravityAssist( swingbyBodyGravitationalParameter_,
                                               departureBodyVelocity_,
                                               ( *velocityBeforeDepartureBodyPtr_ ),
                                               velocityAfterDeparture_,
                                               minimumPericenterRadius_ );

    // Return the deltaV
    deltaV = deltaV_;
}

//! Calculate intermediate positions and their corresponding times.
void SwingbyLegMga::intermediatePoints( const double maximumTimeStep,
                                        std::vector < Eigen::Vector3d >& positionVector,
                                        std::vector < double >& timeVector,
                                        const double startingTime )
{
    // Test if the trajectory has already been calculated.
    if ( std::isnan( velocityAfterDeparture_( 0 ) ) )
    {
        // If the velocity after departure has not been set yet, the trajectory has not been
        // calculated yet and hence still needs to be calculated, which is done below.
        Eigen::Vector3d tempVelocityBeforeArrivalBody;
        double tempDeltaV;
        calculateLeg( tempVelocityBeforeArrivalBody, tempDeltaV );
    }

    // Store the initial state.
    Eigen::VectorXd initialState ( 6 );
    initialState.segment( 0, 3 ) = departureBodyPosition_;
    initialState.segment( 3, 3 ) = velocityAfterDeparture_;

    // Call the trajectory return method to obtain the intermediate points along the trajectory.
    returnTrajectory( initialState, centralBodyGravitationalParameter_, timeOfFlight_,
                      maximumTimeStep, positionVector, timeVector, startingTime );
}

//! Return maneuvres along the leg.
void SwingbyLegMga::maneuvers( std::vector < Eigen::Vector3d >& positionVector,
                               std::vector < double >& timeVector,
                               std::vector < double >& deltaVVector,
                               const double startingTime )
{
    // Test if the trajectory has already been calculated.
    if ( std::isnan( velocityAfterDeparture_( 0 ) ) )
    {
        // If the velocity after departure has not been set yet, the trajectory has not been
        // calculated yet and hence still needs to be calculated, which is done below.
        Eigen::Vector3d tempVelocityBeforeArrivalBody;
        double tempDeltaV;
        calculateLeg( tempVelocityBeforeArrivalBody, tempDeltaV );
    }

    // Resize vectors to the correct size.
    positionVector.resize( 1 );
    timeVector.resize( 1 );
    deltaVVector.resize( 1 );

    // Assign correct values to the vectors.
    positionVector[ 0 ] = departureBodyPosition_;
    timeVector[ 0 ] = 0.0 + startingTime;
    deltaVVector[ 0 ] = deltaV_;
}

//! Update the defining variables.
void SwingbyLegMga::updateDefiningVariables( const Eigen::VectorXd& variableVector )
{
    timeOfFlight_ = variableVector[ 0 ];
}

} // namespace spaceTrajectories
} // namespace tudat

