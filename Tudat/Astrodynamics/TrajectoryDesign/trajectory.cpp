#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/Astrodynamics/TrajectoryDesign/captureLeg.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/departureLegMga.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/departureLegMga1DsmPosition.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/departureLegMga1DsmVelocity.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/planetTrajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLegMga.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLegMga1DsmPosition.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLegMga1DsmVelocity.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Calculate the legs
void Trajectory::calculateTrajectory( double& totalDeltaV )
{
    // Set the deltaV equal to zero.
    totalDeltaV = 0.;

    // Loop through all the interplanetary legs and update the deltaV.
    for ( int counter = 0; counter < numberOfLegs_; counter++ )
    {
        missionLegPtrVector_[ counter ]->calculateLeg(
                    *spacecraftVelocityPtrVector_[ counter ], deltaVVector_[ counter ] );

        totalDeltaV += deltaVVector_[ counter ];
    }
}

//! Returns intermediate points along the trajectory.
void Trajectory::intermediatePoints( double maximumTimeStep,
                                     std::vector < Eigen::Vector3d >& positionVector,
                                     std::vector < double >& timeVector )
{
    // Initiate vectors in which which the position and time vectors from the space leg classes
    // can be stored temporarily.
    std::vector < Eigen::Vector3d > temporaryPositions;
    std::vector < double > temporaryTimes;

    // Set the initial position and time.
    positionVector.push_back( planetPositionVector_[ 0 ] );
    timeVector.push_back( trajectoryVariableVector_[ 0 ] );

    // Set the time variable, which keeps track of the overall time in the trajectory. Necessary
    // because all the legs only keep track of the time relative to the start of the leg.
    double time = 0;

    // Loop through all the interplanetary legs and add the trajectory plots of each part of the
    // trajectory to the total position and time vectors.
    for ( int counter = 0; counter < numberOfLegs_; counter++ )
    {
        // Update the time at the start of the leg.
        time += trajectoryVariableVector_[ counter ];

        // Obtain the intermediate points of the leg.
        missionLegPtrVector_[ counter ]->intermediatePoints( maximumTimeStep,
                                                                 temporaryPositions,
                                                                 temporaryTimes, time );

        // Add the vectors of this leg to tose of the entire trajectory.
        positionVector.insert( positionVector.end( ), temporaryPositions.begin( ) + 1,
                               temporaryPositions.end( ) );
        timeVector.insert( timeVector.end( ), temporaryTimes.begin( ) + 1, temporaryTimes.end( ) );
    }
}

//! Return maneuvres along the trajectory.
void Trajectory::maneuvers( std::vector < Eigen::Vector3d >& positionVector,
                            std::vector < double >& timeVector,
                            std::vector < double >& deltaVVector )
{
    // Empty the vectors.
    positionVector.resize( 0 );
    timeVector.resize( 0 );
    deltaVVector.resize( 0 );

    // Initiate vectors for temporarily storing the variables from the interplanetary legs.
    std::vector < Eigen::Vector3d > temporaryPositions;
    std::vector < double > temporaryTimes;
    std::vector < double > temporaryDeltaVs;

    // Set the time variable, which keeps track of the overall time in the trajectory. Necessary
    // because all the legs only keep track of the time relative to the start of the leg.
    double time = 0;

    // Loop through all the interplanetary legs and add the maneuver information to the
    // corresponding vectors.
    for ( unsigned int counter = 0; counter < numberOfLegs_; counter++ )
    {
        // Update the time at the start of the leg.
        time += trajectoryVariableVector_[ counter ];

        // Obtain the intermediate points of the leg.
        missionLegPtrVector_[ counter ]->maneuvers( temporaryPositions, temporaryTimes,
                                                        temporaryDeltaVs, time );

        // Add the vectors of this leg to tose of the entire trajectory.
        positionVector.insert( positionVector.end( ), temporaryPositions.begin( ),
                               temporaryPositions.end( ) );
        timeVector.insert( timeVector.end( ), temporaryTimes.begin( ), temporaryTimes.end( ) );
        deltaVVector.insert( deltaVVector.end( ), temporaryDeltaVs.begin( ),
                             temporaryDeltaVs.end( ) );
    }
}

//! Return planetary orbits.
void Trajectory::planetaryOrbits( double maximumTimeStep,
                                  std::vector< std::vector < Eigen::Vector3d > >&
                                        positionVectorVector,
                                  std::vector< std::vector < double > >& timeVectorVector )
{
    // Resize the corresponding vectors.
    positionVectorVector.resize( numberOfLegs_ + 1 );
    timeVectorVector.resize( numberOfLegs_ + 1 );

    // Initiate time variable in MJD to extract the ephemeris.
    double timeMJD = 0.;

    for ( int counter = 0; counter < numberOfLegs_ ; counter++ )
    {
        // Update the time to visit at a planet.
        timeMJD = timeMJD + trajectoryVariableVector_[ counter ] / physical_constants::JULIAN_DAY;

        // Obtain position and time vectors for this planet.
        returnSingleRevolutionPlanetTrajectory( ephemerisVector_[ counter ],
                                                centralBodyGravitationalParameter_, timeMJD,
                                                maximumTimeStep, positionVectorVector[ counter ],
                                                timeVectorVector[ counter ],
                                                timeMJD * physical_constants::JULIAN_DAY );
    }
}

//! Return planetary encounters.
void Trajectory::planetaryEncounters( std::vector < Eigen::Vector3d >& positionVector,
                          std::vector < double >& timeVector )
{
    // A check should be added to see if the trajectory has been calculated already.

    // Resize vectors.
    positionVector.resize( numberOfLegs_ + 1 );
    timeVector.resize( numberOfLegs_ + 1 );

    // Initiate a time variable.
    double time = 0.;

    // Start a loop in which the positions and vectors are obtained.
    for ( int counter = 0; counter < numberOfLegs_ + 1; counter++ )
    {
        time += trajectoryVariableVector_[ counter ];
        positionVector[ counter ] = planetPositionVector_[ counter ];
        timeVector[ counter ] = time;
    }
}

//! Update the ephemeris.
void Trajectory::updateEphemeris( )
{
    // Calculate the ephemeris and store it in the corresponding variables in this class.
    extractEphemeris( );

    // Loop through all the mission legs and update their ephemeris variables.
    for ( int counter = 0; counter < numberOfLegs_; counter++ )
    {
        missionLegPtrVector_[ counter ]->updateEphemeris( planetPositionVector_[ counter ],
                                                              planetPositionVector_[ counter + 1],
                                                              planetVelocityVector_[ counter ] );
    }
}

    //! Update the variable vector.
void Trajectory::updateVariableVector( const Eigen::VectorXd& trajectoryVariableVector )
{
    // Store the new vector in the trajectory class.
    trajectoryVariableVector_ = trajectoryVariableVector;

    // Variable that counts the number of additional variables (apart from the timing variables)
    // have been called so far. This is required to keep track of the locations of the correct
    // variables describing trajectories including DSMs.
    int additionalVariableCounter = 0;

    // Initialize a vector in which the variables that have to be passed will be stored temporarily.
    Eigen::VectorXd tempVector;

    // Loop through all the mission legs and update their defining variables.
    for ( int counter = 0; counter < numberOfLegs_; counter++ )
    {
        switch ( legTypeVector_[ counter ] )
        {
            case mga_Departure: case mga_Swingby: case capture:
                tempVector.resize( 1 );
                tempVector << trajectoryVariableVector_[ 1 /*jump over t_0*/ + counter ];
                break;
            case mga1DsmPosition_Departure: case mga1DsmPosition_Swingby:
            case mga1DsmVelocity_Departure: case mga1DsmVelocity_Swingby:
                tempVector.resize( 5 );
                tempVector << trajectoryVariableVector_[ 1 + counter ],
                              trajectoryVariableVector_.segment( 1 + numberOfLegs_ +
                                                                 additionalVariableCounter, 4 );
                additionalVariableCounter += 4;
                break;
        }
        missionLegPtrVector_[ counter ]->updateDefiningVariables( tempVector );
    }
}

//! Test to check if the size of all input parameters is correct.
bool Trajectory::incorrectSize( )
{
    // Check the legTypeVector is of the correct size.
    if ( legTypeVector_.size( ) != numberOfLegs_ )
    {
        std::cerr << "\nIncorrect size of legtypeVector.";
        return true;
    }

    // Check the planetVector is of the correct size.
    if ( ephemerisVector_.size( ) != numberOfLegs_  )
    {
        std::cerr << "\nIncorrect size of planetVector.";
        return true;
    }

    // Check the gravitationalParameterVector is of the correct size.
    if ( gravitationalParameterVector_.size( ) != numberOfLegs_ )
    {
        std::cerr << "\nIncorrect size of gravitationalParameterVector.";
        return true;
    }

    // Check the size of the trajectoryDefiningParameterVector.
    if ( trajectoryVariableVector_.size( ) != checkTrajectoryVariableVectorSize( ) )
    {
        std::cerr << "\nIncorrect size of trajectoryDefiningParameterVector. "<<checkTrajectoryVariableVectorSize( )<<" "<<
                     trajectoryVariableVector_.size( )<<std::endl;;
        return true;
    }
    return false;
}

//! Prepare velocity and position vectors.
void Trajectory::prepareVelocityAndPositionVectors( )
{
    planetPositionVector_.resize( numberOfLegs_ + 1 );
    planetVelocityVector_.resize( numberOfLegs_ + 1 );
    spacecraftVelocityPtrVector_.resize( numberOfLegs_ );
    deltaVVector_.resize( numberOfLegs_ );

    // Prepare empty contents for the spacecraft velocity vector.
    Eigen::Vector3d temp ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
    for ( int counter = 0; counter < numberOfLegs_; counter++)
    {
        spacecraftVelocityPtrVector_[ counter ] = std::make_shared< Eigen::Vector3d > ( temp );
    }
}

int Trajectory::checkTrajectoryVariableVectorSize( )
{
    // The size is always 1, which is the departure epoch.
    int size = 1;

    // Go through all the legs and add the appropriate amount of additional variables.
    for ( int counter = 0; counter < numberOfLegs_; counter++ )
    {
        switch ( legTypeVector_[ counter ] )
        {
            case mga_Departure: case mga_Swingby: case capture:
                size += 1;
                break;
            case mga1DsmPosition_Departure: case mga1DsmPosition_Swingby:
            case mga1DsmVelocity_Departure: case mga1DsmVelocity_Swingby:
                size += 5;
                break;
        }
    }
    return size;
}

//! Prepare the legs and link the variables
void Trajectory::prepareLegs( )
{
    missionLegPtrVector_.resize( numberOfLegs_ );

    // Variable that counts the number of additional variables (apart from the timing
    // variables) have been called so far. This is required to keep track of the locations of
    // the correct variables describing trajectories including DSMs.
    int additionalVariableCounter = 0;

    // A counter that keeps track of the number of departures and captures that have been performed.
    int departureOrCaptureCounter = 0;

    // Start main loop, which adds a new interplanetary leg in each iteration.
    for ( int counter = 0; counter < numberOfLegs_; counter++)
    {
        // Variable that saves the current leg.
        std::shared_ptr< MissionLeg > missionLeg;

        // Depending on the leg type, a different interplanetary leg is to be added. In this switch
        // structure this leg is added and the accompanying variables are linked.
        switch ( legTypeVector_[ counter ] )
        {
            case mga_Departure:
            {
                // Initialize leg with the corresponding variables/pointers.
                std::shared_ptr< DepartureLegMga > departureLegMga =
                        std::make_shared < DepartureLegMga >(
                             planetPositionVector_[ counter ],
                             planetPositionVector_[ counter + 1],
                             trajectoryVariableVector_[ counter + 1 ],
                             planetVelocityVector_[ counter ],
                             centralBodyGravitationalParameter_,
                             gravitationalParameterVector_[ counter ],
                             semiMajorAxesVector_[ departureOrCaptureCounter ],
                             eccentricityVector_[ departureOrCaptureCounter ],
                        includeDepartureDeltaV_ );

                missionLeg = departureLegMga;
                departureLeg_ = departureLegMga;

                // Update the departure and capture counter.
                departureOrCaptureCounter++;
                break;
            }

            case mga_Swingby:
            {
                // Initialize leg with the corresponding variables/pointers.
                std::shared_ptr< SwingbyLegMga > swingbyLegMga =
                        std::make_shared< SwingbyLegMga >(
                           planetPositionVector_[ counter ],
                           planetPositionVector_[ counter + 1],
                           trajectoryVariableVector_[ counter + 1 ],
                           planetVelocityVector_[ counter ],
                           centralBodyGravitationalParameter_,
                           gravitationalParameterVector_[ counter ],
                           spacecraftVelocityPtrVector_[ counter - 1],
                           minimumPericenterRadiiVector_[ counter ] );

                missionLeg = swingbyLegMga;

                break;
            }

            case mga1DsmPosition_Departure:
            {
                // Initialize leg with the corresponding variables/pointers.
                std::shared_ptr< DepartureLegMga1DsmPosition > departureLegMga1DsmPosition =
                        std::make_shared< DepartureLegMga1DsmPosition >(
                             planetPositionVector_[ counter ],
                             planetPositionVector_[ counter + 1],
                             trajectoryVariableVector_[ counter + 1 ],
                             planetVelocityVector_[ counter ],
                             centralBodyGravitationalParameter_,
                             gravitationalParameterVector_[ counter ],
                             semiMajorAxesVector_[ departureOrCaptureCounter ],
                             eccentricityVector_[ departureOrCaptureCounter ],
                             trajectoryVariableVector_[ numberOfLegs_ +
                                   additionalVariableCounter + 1 ],
                             trajectoryVariableVector_[ numberOfLegs_ +
                                   additionalVariableCounter + 2 ],
                             trajectoryVariableVector_[ numberOfLegs_ +
                                   additionalVariableCounter + 3 ],
                             trajectoryVariableVector_[ numberOfLegs_ +
                                   additionalVariableCounter + 4 ],
                        includeDepartureDeltaV_  );


                missionLeg = departureLegMga1DsmPosition;
                departureLeg_ = departureLegMga1DsmPosition;

                // Update the additional variable counter
                additionalVariableCounter += 4;

                // Update the departure and capture counter.
                departureOrCaptureCounter++;
                break;
            }

            case mga1DsmPosition_Swingby:
            {
                // Initialize leg with the corresponding variables/pointers.
                std::shared_ptr< SwingbyLegMga1DsmPosition > swingbyLegMga1DsmPosition =
                        std::make_shared< SwingbyLegMga1DsmPosition >(
                           planetPositionVector_[ counter ],
                           planetPositionVector_[ counter + 1],
                           trajectoryVariableVector_[ counter + 1 ],
                           planetVelocityVector_[ counter ],
                           centralBodyGravitationalParameter_,
                           gravitationalParameterVector_[ counter ],
                           spacecraftVelocityPtrVector_[ counter - 1],
                           minimumPericenterRadiiVector_[ counter ],
                           trajectoryVariableVector_[ numberOfLegs_ +
                                   additionalVariableCounter + 1 ],
                           trajectoryVariableVector_[ numberOfLegs_ +
                                   additionalVariableCounter + 2 ],
                           trajectoryVariableVector_[ numberOfLegs_ +
                                   additionalVariableCounter + 3 ],
                           trajectoryVariableVector_[ numberOfLegs_ +
                                   additionalVariableCounter + 4 ] );

                missionLeg = swingbyLegMga1DsmPosition;

                // Update the additional variable counter
                additionalVariableCounter += 4;
                break;
            }

            case mga1DsmVelocity_Departure:
            {
                // Initialize leg with the corresponding variables/pointers.
                std::shared_ptr< DepartureLegMga1DsmVelocity > departureLegMga1DsmVelocity =
                        std::make_shared< DepartureLegMga1DsmVelocity >(
                             planetPositionVector_[ counter ],
                             planetPositionVector_[ counter + 1],
                             trajectoryVariableVector_[ counter + 1 ],
                             planetVelocityVector_[ counter ],
                             centralBodyGravitationalParameter_,
                             gravitationalParameterVector_[ counter ],
                             semiMajorAxesVector_[ departureOrCaptureCounter ],
                             eccentricityVector_[ departureOrCaptureCounter ],
                             trajectoryVariableVector_[ numberOfLegs_ +
                                    additionalVariableCounter + 1 ],
                             trajectoryVariableVector_[ numberOfLegs_ +
                                    additionalVariableCounter + 2 ],
                             trajectoryVariableVector_[ numberOfLegs_ +
                                    additionalVariableCounter + 3 ],
                             trajectoryVariableVector_[ numberOfLegs_ +
                                    additionalVariableCounter + 4 ],
                        includeDepartureDeltaV_  );

                missionLeg = departureLegMga1DsmVelocity;
                departureLeg_ = departureLegMga1DsmVelocity;

                // Update the additional variable counter
                additionalVariableCounter += 4;

                // Update the departure and capture counter.
                departureOrCaptureCounter++;
                break;
            }

            case mga1DsmVelocity_Swingby:
            {
                // Initialize leg with the corresponding variables/pointers.
                std::shared_ptr< SwingbyLegMga1DsmVelocity > swingbyLegMga1DsmVelocity =
                        std::make_shared< SwingbyLegMga1DsmVelocity >(
                           planetPositionVector_[ counter ],
                           planetPositionVector_[ counter + 1],
                           trajectoryVariableVector_[ counter + 1 ],
                           planetVelocityVector_[ counter ],
                           centralBodyGravitationalParameter_,
                           gravitationalParameterVector_[ counter ],
                           spacecraftVelocityPtrVector_[ counter - 1],
                           trajectoryVariableVector_[ numberOfLegs_ +
                                  additionalVariableCounter + 1 ],
                           trajectoryVariableVector_[ numberOfLegs_ +
                                  additionalVariableCounter + 2 ],
                           trajectoryVariableVector_[ numberOfLegs_ +
                                  additionalVariableCounter + 3 ],
                           trajectoryVariableVector_[ numberOfLegs_ +
                                  additionalVariableCounter + 4 ] );

                missionLeg = swingbyLegMga1DsmVelocity;

                // Update the additional variable counter
                additionalVariableCounter += 4;
                break;
            }

            case capture:
            {
                // Initialize leg with the corresponding variables/pointers.
                std::shared_ptr< CaptureLeg > captureLeg =
                        std::make_shared< CaptureLeg >(
                             planetPositionVector_[ counter ],
                             trajectoryVariableVector_[ counter + 1 ],
                             planetVelocityVector_[ counter ],
                             centralBodyGravitationalParameter_,
                             gravitationalParameterVector_[ counter ],
                             spacecraftVelocityPtrVector_[ counter - 1],
                             semiMajorAxesVector_[ departureOrCaptureCounter ],
                             eccentricityVector_[ departureOrCaptureCounter ],
                        includeArrivalDeltaV_ );

                missionLeg = captureLeg;
                captureLeg_ = captureLeg;

                // Update the departure and capture counter.
                departureOrCaptureCounter++;
                break;
            }
            default:
            {
                std::cerr<<"Trajectory Model "<<counter<<" does not exist.";
            }

        }
        missionLegPtrVector_[ counter ] = missionLeg;

    }
}

//! Extract the ephemeris data.
void Trajectory::extractEphemeris( )
{
    // Initiate a timing variable.
    double time = 0.0;

    // Obtain positions and velocities at the corresponding times, by extracting the ephemeris data
    // at the visitation times.
    for ( int counter = 0; counter < numberOfLegs_ ; counter++ )
    {
        // Update the time to visit at a planet.
        time = time + trajectoryVariableVector_[ counter ];

        // Calculate Cartesian state of Planets
        temporaryCartesianElements_ = ephemerisVector_[ counter ]->getCartesianState( time );

        // Set planet position and velocity from the Cartesian elements.
        planetPositionVector_[ counter ] = temporaryCartesianElements_.segment( 0, 3 );
        planetVelocityVector_[ counter ] = temporaryCartesianElements_.segment( 3, 3 );
    }
}

//! Return launch conditions.
void Trajectory::getLaunchConditions( Eigen::Vector3d& departureBodyPosition,
                                      Eigen::Vector3d& departureBodyVelocity,
                                      Eigen::Vector3d& velocityAfterDeparture )
{
    missionLegPtrVector_[ 0 ]->returnDepartureVariables( departureBodyPosition,
                                                             departureBodyVelocity,
                                                             velocityAfterDeparture );
}

} // namespace transfer_trajectories
} // namespace tudat
