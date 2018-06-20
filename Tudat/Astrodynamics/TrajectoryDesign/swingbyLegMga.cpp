#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/Astrodynamics/MissionSegments/gravityAssist.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"

#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLegMga.h"

namespace tudat
{
namespace transfer_trajectories
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

} // namespace transfer_trajectories
} // namespace tudat

