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
 *       Vinko, T. and Izzo, D. (2008). Global Optimisation Heuristics and Test Problems for
 *          Preliminary Spacecraft Trajectory Design. Technical Report, Advanced Concept Team, ESA.
 *
 */

#include "Tudat/Basics/basicTypedefs.h"

#include "Eigen/Dense"

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"

#include "departureLegMga1DsmVelocity.h"
#include "exportTrajectory.h"

namespace tudat
{
namespace spaceTrajectories
{

//! Calculate the leg and update the Delta V and the velocity before the next body.
void DepartureLegMga1DsmVelocity::calculateLeg( Eigen::Vector3d& velocityBeforeArrivalBody,
                                                double& deltaV )
{
    // Calculate the DSM time of application from the time of flight fraction
    dsmTime_ = dsmTimeOfFlightFraction_ * timeOfFlight_;

    // Calculate unit vectors as described in [Vinko and Izzo, 2008].
    const Eigen::Vector3d unitVector1 = departureBodyVelocity_ / departureBodyVelocity_.norm( );
    const Eigen::Vector3d unitVector3 = departureBodyPosition_.cross( departureBodyVelocity_ ) /
                                  departureBodyPosition_.cross( departureBodyVelocity_ ).norm( );
    const Eigen::Vector3d unitVector2 = unitVector3.cross( unitVector1 );

    // Calculate the velocity after departure as described in [Vinko and Izzo, 2008].
    // First add the departure body velocity
    velocityAfterDeparture_ = departureBodyVelocity_ +
    // Then add the velocity in the direction of the departure velocity
            excessVelocityMagnitude_ * cos( excessVelocityInPlaneAngle_ ) *
            cos( excessVelocityOutOfPlaneAngle_ ) * unitVector1 +
    // Then add the velocity in the other direction of the 2D plane of the departure velocity
            excessVelocityMagnitude_ * sin( excessVelocityInPlaneAngle_ ) *
            cos( excessVelocityOutOfPlaneAngle_ ) * unitVector2 +
    // Finally add the 3D component
            excessVelocityMagnitude_ * sin( excessVelocityOutOfPlaneAngle_ ) * unitVector3;

    // Transfer the initial position and velocity into a vectorXd object with cartesian coordinates.
    Eigen::Vector6d cartesianElements (6), keplerianElements (6);
    cartesianElements.segment( 0, 3 ) = departureBodyPosition_;
    cartesianElements.segment( 3, 3 ) = velocityAfterDeparture_;

    // Convert the cartesian elements into keplerian elements.
    keplerianElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                cartesianElements, centralBodyGravitationalParameter_ );

    // Propagate the keplerian elements until the moment of application of the DSM.
    keplerianElements = orbital_element_conversions::propagateKeplerOrbit( keplerianElements,
                dsmTime_, centralBodyGravitationalParameter_ );
    // Convert the keplerian elements back into Cartesian elements.
    cartesianElements = orbital_element_conversions::convertKeplerianToCartesianElements(
                keplerianElements, centralBodyGravitationalParameter_ );

    // Set the corresponding position and velocity vectors.
    dsmLocation_ = cartesianElements.segment( 0, 3 );
    velocityBeforeDsm_ = cartesianElements.segment( 3, 3 );

    // Calculate the velocities after the DSM and before the arrival body.
    mission_segments::solveLambertProblemIzzo( dsmLocation_, arrivalBodyPosition_, timeOfFlight_ -
                                               dsmTime_, centralBodyGravitationalParameter_,
                                               velocityAfterDsm_, velocityBeforeArrivalBody );

    // Calculate the deltaV originating from the departure maneuver and the DSM.
    deltaVDeparture_ = mission_segments::computeEscapeOrCaptureDeltaV(
                departureBodyGravitationalParameter_, semiMajorAxis_, eccentricity_,
                excessVelocityMagnitude_ );

    deltaVDsm_ = ( velocityAfterDsm_ - velocityBeforeDsm_ ).norm( );

    //Calculate the total deltaV.
    deltaV_ = deltaVDeparture_ + deltaVDsm_;

    // Return the deltaV
    deltaV = deltaV_;
}

//! Calculate intermediate positions and their corresponding times.
void DepartureLegMga1DsmVelocity::intermediatePoints( const double maximumTimeStep,
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

    // Call the trajectory return method to obtain the intermediate points along the first part of
    // the leg.
    returnTrajectory( initialState, centralBodyGravitationalParameter_, dsmTime_,
                      maximumTimeStep, positionVector, timeVector, startingTime );

    // Make vectors that can be used to call the trajectory return method for the second part of
    // the leg.
    std::vector < Eigen::Vector3d > temporaryPositions;
    std::vector < double > temporaryTimes;

    // Store the state at the DSM.
    initialState.segment( 0, 3 ) = dsmLocation_;
    initialState.segment( 3, 3 ) = velocityAfterDsm_;

    // Call the trajectory return method to obtain the intermediate points along the first part of
    // the leg.
    returnTrajectory( initialState, centralBodyGravitationalParameter_, timeOfFlight_ - dsmTime_,
                      maximumTimeStep, temporaryPositions, temporaryTimes,
                      startingTime + dsmTime_);

    // Add the vectors of this part to those of the entire leg.
    positionVector.insert( positionVector.end( ), temporaryPositions.begin( ) + 1,
                           temporaryPositions.end( ) );
    timeVector.insert( timeVector.end( ), temporaryTimes.begin( ) + 1, temporaryTimes.end( ) );
}

//! Return maneuvres along the leg.
void DepartureLegMga1DsmVelocity::maneuvers( std::vector < Eigen::Vector3d >& positionVector,
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
    positionVector.resize( 2 );
    timeVector.resize( 2 );
    deltaVVector.resize( 2 );

    // Assign correct values to the vectors.
    positionVector[ 0 ] = departureBodyPosition_;
    timeVector[ 0 ] = 0.0 + startingTime;
    deltaVVector[ 0 ] = deltaVDeparture_;
    positionVector[ 1 ] = dsmLocation_;
    timeVector[ 1 ] = dsmTime_ + startingTime;
    deltaVVector[ 1 ] = deltaVDsm_;
}

//! Update the defining variables.
void DepartureLegMga1DsmVelocity::updateDefiningVariables( const Eigen::VectorXd& variableVector )
{
    timeOfFlight_ = variableVector[ 0 ];
    dsmTimeOfFlightFraction_ = variableVector[ 1 ];
    excessVelocityMagnitude_ = variableVector[ 2 ];
    excessVelocityInPlaneAngle_ = variableVector[ 3 ];
    excessVelocityOutOfPlaneAngle_ = variableVector[ 4 ];
}

} // namespace spaceTrajectories
} // namespace tudat
