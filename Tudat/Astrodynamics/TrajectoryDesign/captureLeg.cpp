#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"

#include "Tudat/Astrodynamics/TrajectoryDesign/captureLeg.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Calculate the leg and update the Delta V and the velocity before the next body.
void CaptureLeg::calculateLeg( Eigen::Vector3d& velocityBeforeArrivalBody,
                               double& deltaV )
{
    // This velocity does not have physical meaning in this leg. (Should be programmed differently)
    velocityBeforeArrivalBody << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN;

    // Set the velocity after departure equal to that of the departure body. This is maily done to
    // flag the fact that the leg has been calculated. (Should be programmed differently)
    velocityAfterDeparture_ = departureBodyVelocity_;

    // Calculate the required deltaV for capture.
    if( includeArrivalDeltaV_ )
    {
        deltaV_ = mission_segments::computeEscapeOrCaptureDeltaV( captureBodyGravitationalParameter_,
                                                                  semiMajorAxis_, eccentricity_,
                                                                  ( *velocityBeforeDepartureBodyPtr_ -
                                                                    departureBodyVelocity_ ).norm( ) );    }
    else
    {
        deltaV_ = 0.0;
    }

    // Return the deltaV
    deltaV = deltaV_;
}

//! Calculate intermediate positions and their corresponding times.
void CaptureLeg::intermediatePoints( const double maximumTimeStep,
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
void CaptureLeg::maneuvers( std::vector < Eigen::Vector3d >& positionVector,
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

void CaptureLeg::updateDefiningVariables( const Eigen::VectorXd& variableVector )
{
    timeOfFlight_ = variableVector[ 0 ];
}

} // namespace transfer_trajectories
} // namespace tudat
