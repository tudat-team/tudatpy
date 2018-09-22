// The ephemeris extractor now takes planet ephemeris as only means of extracting the data. This
//      could be improved.

#include <cmath>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/planetTrajectory.h"

//#include "Temporary/temporaryFunctions.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Return a vector of positions and times from ephemeris data for a certain epoch and duration.
void returnPlanetTrajectory( const ephemerides::EphemerisPointer & ephemerisPtr,
                             const double centralBodyGravitationalParameter,
                             const double startingEpochMJD2000,
                             const double duration,
                             const double maximumTimeStep,
                             std::vector < Eigen::Vector3d >& positionVector,
                             std::vector < double >& timeVector,
                             const double startingTime )
{

    double MJD2000 = 51544.5;

    // Calculate the number of intermediate steps required and the corresponding time step.
    int numberOfSteps = std::ceil( duration / maximumTimeStep );
    double timeStep = duration / numberOfSteps;

    // Resize the position and time vectors. Note that the "+ 1" comes from the end points.
    positionVector.resize( numberOfSteps + 1 );
    timeVector.resize( numberOfSteps + 1 );

    // Initiate variables to temporarily store the ephemeris data.
    Eigen::Vector6d temporaryCartesianElements;

    // Initiate time variable in JD for extracting data, and a time variable in seconds to store
    // the time in the vector.
    double timeMJD2000 = startingEpochMJD2000 + MJD2000;
    double time = 0.0;

//    std::shared_ptr<ephemerides::ApproximatePlanetPositions> approxEphemerisPtr =
//            std::dynamic_pointer_cast<ephemerides::ApproximatePlanetPositions>(ephemerisPtr);
//    double timeJD2000;
//    double timeSecondsSinceEpoch;

    // Start the for loop to obtain positions along the trajectory.
    for ( int counter = 0; counter < numberOfSteps + 1; counter++ )
    {
        // Get Keplerian state of the planet at the corresponding time.
        temporaryCartesianElements = ephemerisPtr->getCartesianState( time );


        // Store the position and time
        positionVector[ counter ] = temporaryCartesianElements.segment( 0, 3 );
        timeVector[ counter ] = time + startingTime;

        // Update the times.
        time += timeStep;
        timeMJD2000 += timeStep / physical_constants::JULIAN_DAY;
    }
}


void returnSingleRevolutionPlanetTrajectory(const ephemerides::EphemerisPointer &ephemerisPtr,
        const double centralBodyGravitationalParameter,
        const double startingEpochMJD2000,
        const double maximumTimeStep,
        std::vector < Eigen::Vector3d >& positionVector,
        std::vector < double >& timeVector,
        const double startingTime )
{
    // Initiate variables to temporarily store the ephemeris data.
    Eigen::Vector6d initialCartesianElements;

    // Convert MJD2000 to seconds since epoch
    double startingEpochMJD = startingEpochMJD2000 + 51544.5;
    double timeJD2000 = basic_astrodynamics::convertModifiedJulianDayToJulianDay(startingEpochMJD);
    double timeSecondsSinceEpoch = basic_astrodynamics::convertJulianDayToSecondsSinceEpoch(timeJD2000);

    // Get Cartesian state of the planet at the corresponding time.
    initialCartesianElements = ephemerisPtr->getCartesianState( timeSecondsSinceEpoch );

    // Compute duration (one revolution).
    double period = 2 * mathematical_constants::PI / std::sqrt( centralBodyGravitationalParameter /
            std::pow( initialCartesianElements.norm(), 3 ) );

    // Compute position and time vector, by pasing computed data to other function.
    returnTrajectory( initialCartesianElements, centralBodyGravitationalParameter, period,
                      maximumTimeStep, positionVector, timeVector, startingTime );

}

} // namespace transfer_trajectories
} // namespace tudat

