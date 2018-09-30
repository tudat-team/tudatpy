/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      More information on the trajectory design code and how the quantities
 *      are calculated can be found in [Musegaas, 2012], who is also the author
 *      of this code
 *
*/

// Some improvements are foreseen:
// * The application of legs departure vs swingby vs capture is insecure. This could be improved.
// * The minimum pericenter radii are not necessary for the velocity formulation model. The vector
//      that contains these radii does contain a value for it though. This could be improved.
// * The checks that are performed are not sufficient and by inputting incomplete information,
//      the function may crash/behave strangely.
// * Variables could be passed alternatively by data containers of some kind.
// * When returning the maneuvers, no check is present to check if the trajectory has been
//      calculated already. This could be improved.
// * Default constructor may be regarded as unsafe. Other implementations may be possible.

#ifndef TUDAT_TRAJECTORY_H
#define TUDAT_TRAJECTORY_H

#include <vector>

#include <boost/make_shared.hpp>
#include <memory>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/spaceLeg.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/departureLeg.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/captureLeg.h"

namespace tudat
{

namespace transfer_trajectories
{

// Enumeration containing the different leg types that can be part of the trajectory.
enum legTypes
{
    mga_Departure = 1,
    mga_Swingby,
    mga1DsmPosition_Departure,
    mga1DsmPosition_Swingby,
    mga1DsmVelocity_Departure,
    mga1DsmVelocity_Swingby,
    capture
};

//! Base class for computation of trajectories
/*!
 * This class can compute entire space trajectories of various kinds by concatenating different
 * mission legs. It can compute the total deltaV, but also contains functions that can be used
 * for plotting the trajectory and that of the bodies involved in the trajectory.
 */
class Trajectory
{
public:

    //! Constructor with immediate definition of parameters.
    /*!
     *  Constructor with immediate definition of parameters.
     *  \param numberOfLegs the number of legs in the trajectory.
     *  \param legTypeVector vector containing the leg types.
     *  \param ephemerisVector vector of ephemeris pointers to the different planets.
     *  \param gravitationalParameterVector vector of the gravitational parameters of the visited planets.
     *  \param trajectoryVariableVector vector containing all the defining variables for the whole trajectory.
     *  \param centralBodyGravitationalParameter gravitational parameter of the central body.
     *  \param minimumPericenterRadiiVector vector containing the minimum distance between the spacecraft and body.
     *  \param semiMajorAxesVector vector containing the semi-major axes for the departure and capture leg.
     *  \param eccentricityVector vector containing the eccentricities for the departure and capture leg.
     *  \param includeDepartureDeltaV Boolean denoting whether to include the Delta V of departure.
     *  \param includeArrivalDeltaV Boolean denoting whether to include the Delta V of arrival.
     */
    Trajectory( const double numberOfLegs,
                const std::vector< int >& legTypeVector,
                const std::vector< ephemerides::EphemerisPointer >&
                        ephemerisVector,
                const Eigen::VectorXd& gravitationalParameterVector,
                const Eigen::VectorXd& trajectoryVariableVector,
                const double centralBodyGravitationalParameter,
                const Eigen::VectorXd& minimumPericenterRadiiVector,
                const Eigen::VectorXd& semiMajorAxesVector,
                const Eigen::VectorXd& eccentricityVector,
                const bool includeDepartureDeltaV = true,
                const bool includeArrivalDeltaV = true )
    {
        includeDepartureDeltaV_ = includeDepartureDeltaV;
        includeArrivalDeltaV_ = includeArrivalDeltaV;

        // Set the number of legs.
        numberOfLegs_ = numberOfLegs;

        // Set the leg type vector, the planet vector and the gravitational parameter vector.
        legTypeVector_ = legTypeVector;
        ephemerisVector_ = ephemerisVector;
        gravitationalParameterVector_ = gravitationalParameterVector;

        // Set the trajectory defining parameter vector.
        trajectoryVariableVector_ = trajectoryVariableVector;

        // Set the minimum pericenter radii.
        minimumPericenterRadiiVector_ = minimumPericenterRadiiVector;

        // Set the departure and capture orbit properties.
        semiMajorAxesVector_ = semiMajorAxesVector;
        eccentricityVector_ = eccentricityVector;

        // Set the central body gravitational parameter
        centralBodyGravitationalParameter_ = centralBodyGravitationalParameter;

        // Does not fully function yet. Should be improved.
        if ( incorrectSize( ) ) { std::cerr << "\nReturning..."; return; }

        // Prepare all variables and legs.
        prepareVelocityAndPositionVectors( );
        extractEphemeris( );
        prepareLegs( );
    }

    //! Calculate the legs
    /*!
     * Performs all the calculations required for the trajectory.
     * \param totalDeltaV the total delta V needed for the trajectory.
     */
    void calculateTrajectory( double& totalDeltaV );

    //! Function to retrieve the value of the capture Delta V.
    /*!
     *  Function to retrieve the value of the capture Delta V.
     *  \param captureDeltaV Double denoting the value of the capture Delta V.
     *  \return Double denoting the value of the capture Delta V (returned by reference ).
     */
    void getCaptureDeltaV( double& captureDeltaV )
    {
        captureLeg_->getCaptureDeltaV( captureDeltaV );
    }

    //! Function to retrieve the value of the departure Delta V.
    /*!
     *  Function to retrieve the value of the departure Delta V.
     *  \param departureDeltaV Double denoting the value of the departure Delta V.
     *  \return Double denoting the value of the departure Delta V (returned by reference ).
     */
    void getDepartureDeltaV( double& departureDeltaV )
    {
        departureLeg_->getEscapeDeltaV( departureDeltaV );
    }

    //! Return intermediate points along the trajectory.
    /*!
     * Returns intermediate points along the trajectory, which can for instance be used to plot the
     * trajectory.
     *  \param maximumTimeStep the maximum time between two points along the trajectory.
     *  \param positionVector Vector of positions along the orbit, space according to the maximum time step.
     *  \param timeVector The times corresponding to the positions.
     */
     void intermediatePoints( double maximumTimeStep,
                              std::vector < Eigen::Vector3d >& positionVector,
                              std::vector < double >& timeVector );

     //! Return maneuvres along the trajectory.
     /*!
      * Returns the maneuver points, times and sizes along the trajectory.
      *  \param positionVector Vector of the positions of the maneuvers.
      *  \param timeVector The times corresponding to the positions.
      *  \param deltaVVector the delta V required for each maneuver.
     */
     void maneuvers( std::vector < Eigen::Vector3d >& positionVector,
                     std::vector < double >& timeVector,
                     std::vector < double >& deltaVVector );

     //! Return planetary orbits.
     /*!
      * Returns vectors containing planetary orbits, which are simulated using the ephemeris at the
      * time of visitation.
      *  \param maximumTimeStep the maximum time between two points along the planet orbit.
      *  \param positionVectorVector Vector of positions along the orbit, space according to the maximum time step.
      *  \param timeVectorVector The times corresponding to the positions.
      */
     void planetaryOrbits( double maximumTimeStep,
                           std::vector< std::vector < Eigen::Vector3d > >& positionVectorVector,
                           std::vector< std::vector < double > >& timeVectorVector );

     //! Return planetary encounters.
     /*!
       * Returns vectors containing planetary encounters.
       * \param positionVector Vector of positions of the encounters.
       * \param timeVector The times corresponding to the positions.
       */
     void planetaryEncounters( std::vector < Eigen::Vector3d >& positionVector,
                               std::vector < double >& timeVector );

    //! Update the ephemeris.
    /*!
     * Sets all the positions and the velocities of the trajectory class and the underlying mission
     * leg classes to new values, that are obtained by extracting the ephemeris at the times that
     * are now in the trajectory variable vector. This is required for re-using the class, without
     * re-initializing it.
     */
    void updateEphemeris( );

    //! Update the variable vector.
    /*!
     * Sets the trajectory defining variable vector to the newly specified values. Also sets all
     * the defining variables in the underlying mission leg classes to these new values. This is
     * required for re-using the class, without re-initializing it.
     * \param trajectoryVariableVector the new variable vector.
     */
    void updateVariableVector( const Eigen::VectorXd& trajectoryVariableVector );

    //! Return launch conditions.
    /*!
     * Returns the launch conditions, useful if additional information is needed regarding the
     * launch of the spacecraft. This is for instance required for the TandEM problems of GTOP.
     *  \param departureBodyPosition the departure body position.
     *  \param departureBodyVelocity the velocity of the departure body.
     *  \param velocityAfterDeparture the velocity of the spacecraft after the departure maneuver.
     */
    void getLaunchConditions( Eigen::Vector3d& departureBodyPosition,
                              Eigen::Vector3d& departureBodyVelocity,
                              Eigen::Vector3d& velocityAfterDeparture );

protected:

private:

    //! The number of legs in the trajectory.
    /*!
     * The number of legs present in this trajectory.
     */
    double numberOfLegs_;

    //! The leg type vector.
    /*!
     * The different interplanetary leg types for this trajectory. These are stored in this vector,
     * such that the appropriate leg types can be created for the different legs.
     */
    std::vector< int > legTypeVector_;

    //! The vector containing the Ephemeris objects.
    /*!
     * This vector contains all ephemeris objects.
     */
    std::vector< ephemerides::EphemerisPointer > ephemerisVector_;

    //! The variable vector for defining the trajectory.
    /*!
     * All the specific variables for defining the transfers are stored in this vector. Its size
     * is determined using the number of legs and the leg types.
     */
    Eigen::VectorXd trajectoryVariableVector_;

    //! The planet positions vector.
    /*!
     * In this vector the positions of the planets at the visitation times are stored.
     */
    std::vector< Eigen::Vector3d > planetPositionVector_;

    //! The planet velocities vector.
    /*!
     * In this vector the velocities of the planets at the visitation times are stored.
     */
    std::vector< Eigen::Vector3d > planetVelocityVector_;

    //! The spacecraft velocity vector.
    /*!
     * In this vector the velocity of the spacecraft before visiting a new planet are stored.
     */
    std::vector< std::shared_ptr< Eigen::Vector3d > > spacecraftVelocityPtrVector_;

    //! The deltaV vector.
    /*!
     * In this vector the deltaV's associated with each leg are stored.
     */
    std::vector < double > deltaVVector_;

    //! The gravitational parameter vector.
    /*!
     * In this vector the gravitational parameters associated with each leg are stored.
     */
    Eigen::VectorXd gravitationalParameterVector_;

    //! The minimum pericenter radii vector.
    /*!
     * In this vector the minimum pericenter radii associated with each leg are stored.
     */
    Eigen::VectorXd minimumPericenterRadiiVector_;

    //! The semi-major axes vector.
    /*!
     * In this vector the semi-major axes of all the departure and capture orbits are stored.
     */
    Eigen::VectorXd semiMajorAxesVector_;

    //! The eccentricity vector.
    /*!
     * In this vector the eccentricity of all the departure and capture orbits are stored.
     */
    Eigen::VectorXd eccentricityVector_;

    //! The central body gravitational parameter.
    /*!
     * The central body gravitational parameter.
     */
    double centralBodyGravitationalParameter_;

    //! Boolean denoting whether to include the Delta V of departure.
    bool includeDepartureDeltaV_;

    //! Pointer to the departure leg.
    std::shared_ptr< DepartureLeg > departureLeg_;

    //! Boolean denoting whether to include the Delta V of arrival.
    bool includeArrivalDeltaV_;

    //! Pointer to the capture leg.
    std::shared_ptr< CaptureLeg > captureLeg_;

    //! The interplanetary leg vector.
    /*!
     * The vector containing all the different interplanetary legs.
     */
    std::vector< std::shared_ptr< MissionLeg > > missionLegPtrVector_;

    //! Temporary Keplerian elements vector.
    /*!
     * A temporary vector to store Keplerian elements for (repeated) use in the extraction of the
     * ephemeris data.
     */
    Eigen::Vector6d temporaryKeplerianElements_;

    //! Temporary Cartesian elements object.
    /*!
     * A temporary object to store the Cartesian elements for (repeated) use in the extraction of
     * the ephemeris data.
     */
    Eigen::Vector6d temporaryCartesianElements_;

    //! Test to check if the size of all input parameters is correct.
    /*!
     * Computes the correct size for all input parameters. Returns true if they do not correspond
     * to the expected values.
     */
    bool incorrectSize( );

    //! Prepare velocity and position vectors.
    /*!
     * Prepares the velocity and position vectors required to store and pass all the planetary
     * positions and velocities, spacecraft velocities and deltaV sizes. Does not set any values
     * yet, but sizes them correctly based on the number of legs.
     */
    void prepareVelocityAndPositionVectors( );

    //! Check the size of the trajectory defining parameter vector.
    /*!
     * Checks the size of the trajectory defining parameter vector.
     */
    int checkTrajectoryVariableVectorSize( );

    //! Prepare the legs and link the variables
    /*!
     * A vector containing all the interplanetary legs is created here. The variables are linked
     * to their corresponding legs and all these legs are prepared such that they can be
     * calculated through. This only has to happen once, since all variables are pointers.
     */
    void prepareLegs( );

    //! Extract the ephemeris data.
    /*!
     * Extracts the ephemeris data and stores it into the associated position and velocity vectors.
     */
    void extractEphemeris( );

};

} // namespace transfer_trajectories

} // namespace tudat

#endif // TUDAT_TRAJECTORY_H
