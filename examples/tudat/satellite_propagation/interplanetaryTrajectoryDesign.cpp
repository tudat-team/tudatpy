/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <tudat/astro/basic/physicalConstants.h>
#include <tudat/basics/testMacros.h>
#include <tudat/math/basic/mathematicalConstants.h>
#include "tudat/astro/basic/unitConversions.h"
#include <tudat/astro/basic/orbitalElementConversions.h>

#include <tudat/io/basicInputOutput.h>
#include <tudat/io/applicationOutput.h>


#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/trajectory_design/trajectory.h"
#include "tudat/astro/trajectory_design/exportTrajectory.h"
#include "tudat/astro/trajectory_design/planetTrajectory.h"

int main( )
{
    using namespace tudat;
    using namespace tudat::input_output;
    using namespace tudat::input_output::parsed_data_vector_utilities;
    using namespace tudat::transfer_trajectories;

    //////////////////////////////////////////////////////////////////////////
    ////////////////////////// CASSSINI //////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    // Specify required parameters
    // Specify the number of legs and type of legs.
    int numberOfLegs = 6;
    std::vector< TransferLegType > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[ 0 ] = mga_Departure;
    legTypeVector[ 1 ] = mga_Swingby;
    legTypeVector[ 2 ] = mga_Swingby;
    legTypeVector[ 3 ] = mga_Swingby;
    legTypeVector[ 4 ] = mga_Swingby;
    legTypeVector[ 5 ] = capture;

    // Create the ephemeris vector.
    std::vector< ephemerides::EphemerisPointer >
            ephemerisVector( numberOfLegs );
    ephemerisVector[ 0 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 1 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 2 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 3 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 4 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter );
    ephemerisVector[ 5 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn );

    // Create gravitational parameter vector
    Eigen::VectorXd gravitationalParameterVector( numberOfLegs );
    gravitationalParameterVector << 3.9860119e14, 3.24860e14, 3.24860e14, 3.9860119e14, 1.267e17, 3.79e16;

    // Create variable vector.
    Eigen::VectorXd variableVector( numberOfLegs + 1 );
    variableVector << -789.8117, 158.302027105278, 449.385873819743, 54.7489684339665,
            1024.36205846918, 4552.30796805542, 1/*dummy*/;
    variableVector *= physical_constants::JULIAN_DAY;

    // Create departure and capture variables.
    Eigen::VectorXd semiMajorAxes( 2 ), eccentricities( 2 );
    semiMajorAxes << std::numeric_limits< double >::infinity( ), 1.0895e8 / 0.02;
    eccentricities << 0.0, 0.98;

    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;

    // Create minimum pericenter radii vector
    Eigen::VectorXd minimumPericenterRadii( numberOfLegs );
    minimumPericenterRadii << 6778000.0, 6351800.0, 6351800.0, 6778000.0, 600000000.0, 600000000.0;

    // Create the trajectory problem.
    Trajectory Cassini1( numberOfLegs, legTypeVector, ephemerisVector,
                         gravitationalParameterVector, variableVector, sunGravitationalParameter,
                         minimumPericenterRadii, semiMajorAxes, eccentricities );

    // Vectors for the specific maneuvers and the total delta v
    std::vector< Eigen::Vector3d > positionVector;
    std::vector< double > timeVector;
    std::vector< double > deltaVVector;
    double resultingDeltaV;

    // Calculate the orbits
    Cassini1.calculateTrajectory( resultingDeltaV );
    Cassini1.maneuvers( positionVector, timeVector, deltaVVector );

    std::cout << " Cassini Mission: " << std::endl;
    std::cout << " Total Delta V needed: " << resultingDeltaV <<std::endl;
    std::cout << " Time of Earth departure: " << timeVector[ 0 ]
              << ". Delta V needed at Earth: " << deltaVVector[ 0 ] << std::endl;
    std::cout << " Time of Venus visit: " << timeVector[ 1 ]
              << ". Delta V needed at Venus: " << deltaVVector[ 1 ] << std::endl;
    std::cout << " Time of second Venus visit: " << timeVector[ 2 ]
              << ". Delta V needed at Venus: " << deltaVVector[ 2 ] << std::endl;
    std::cout << " Time of Earth visit: " << timeVector[ 3 ]
              << ". Delta V needed at Earth: " << deltaVVector[ 3 ] << std::endl;
    std::cout << " Time of Jupiter visit: " << timeVector[ 4 ]
              << ". Delta V needed at Jupiter: " << deltaVVector[ 4 ] << std::endl;
    std::cout << " Time of Saturn capture: " << timeVector[ 5 ]
              << ". Delta V needed at Saturn: " << deltaVVector[ 5 ] << std::endl;

    std::cout << std::endl << std::endl;

    //////////////////////////////////////////////////////////////////////////
    ////////////////////////// MESSENGER /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    // Specify required parameters
    // Specify the number of legs and type of legs.
    numberOfLegs = 5;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[ 0 ] = mga1DsmVelocity_Departure;
    legTypeVector[ 1 ] = mga1DsmVelocity_Swingby;
    legTypeVector[ 2 ] = mga1DsmVelocity_Swingby;
    legTypeVector[ 3 ] = mga1DsmVelocity_Swingby;
    legTypeVector[ 4 ] = capture;

    // Create the ephemeris vector.
    ephemerisVector.resize( numberOfLegs );
    ephemerisVector[ 0 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 1 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 2 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 3 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 4 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury );

    // Create gravitational parameter vector
    gravitationalParameterVector.resize( numberOfLegs );
    gravitationalParameterVector << 3.9860119e14, 3.9860119e14, 3.24860e14, 3.24860e14, 2.2321e13;

    // Create variable vector.
    variableVector.resize( numberOfLegs /*time of flight*/ + 1 /*start epoch*/ +
                           4 * ( numberOfLegs - 1 ) /*additional variables for model, except the final capture leg*/ );

    // Add the time of flight and start epoch, which are in JD.
    variableVector << 1171.64503236 * physical_constants::JULIAN_DAY,
            399.999999715 * physical_constants::JULIAN_DAY,
            178.372255301 * physical_constants::JULIAN_DAY,
            299.223139512 * physical_constants::JULIAN_DAY,
            180.510754824 * physical_constants::JULIAN_DAY,
            1, // The capture time is irrelevant for the final leg.
            // Add the additional variables.
            0.234594654679, 1408.99421278, 0.37992647165 * 2 * 3.14159265358979,
            std::acos(  2 * 0.498004040298 - 1. ) - 3.14159265358979 / 2, // 1st leg.
            0.0964769387134, 1.35077257078, 1.80629232251 * 6.378e6, 0.0, // 2nd leg.
            0.829948744508, 1.09554368115, 3.04129845698 * 6.052e6, 0.0, // 3rd leg.
            0.317174785637, 1.34317576594, 1.10000000891 * 6.052e6, 0.0; // 4th leg.

    // Create minimum pericenter radii vector
    minimumPericenterRadii.resize( numberOfLegs );
    minimumPericenterRadii << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN;

    // Create departure and capture variables.
    semiMajorAxes << std::numeric_limits< double >::infinity( ),
            std::numeric_limits< double >::infinity( );
    eccentricities << 0.0, 0.0;

    // Create the trajectory problem.
    Trajectory Messenger( numberOfLegs, legTypeVector, ephemerisVector,
                          gravitationalParameterVector, variableVector, sunGravitationalParameter,
                          minimumPericenterRadii, semiMajorAxes, eccentricities );

    // Vectors for the specific maneuvers and the total delta v
    std::vector< Eigen::Vector3d > positionVectorMessenger;
    std::vector< double > timeVectorMessenger;
    std::vector< double > deltaVVectorMessenger;
    double resultingDeltaVMessenger;

    // Calculate the orbits
    Messenger.calculateTrajectory( resultingDeltaVMessenger );
    Messenger.maneuvers( positionVectorMessenger, timeVectorMessenger, deltaVVectorMessenger );

    std::cout << " Messenger Mission: " << std::endl;
    std::cout << " Total Delta V: " << resultingDeltaVMessenger <<std::endl;
    std::cout << " Time of Earth departure: " << timeVectorMessenger[ 0 ]
              << ". Delta V needed at Earth: " << deltaVVectorMessenger[ 0 ] << std::endl;
    std::cout << " Time of 1st DSM: " << timeVectorMessenger[ 1 ]
              << ". Delta V needed for 1st DSM: " << deltaVVectorMessenger[ 1 ] << std::endl;
    std::cout << " Time of second Earth visit: " << timeVectorMessenger[ 2 ]
              << ". Delta V needed at Earth: " << deltaVVectorMessenger[ 2 ] << std::endl;
    std::cout << " Time of 2nd DSM: " << timeVectorMessenger[ 3 ]
              << ". Delta V needed for 2nd DSM: " << deltaVVectorMessenger[ 3 ] << std::endl;
    std::cout << " Time of Venus visit: " << timeVectorMessenger[ 4 ]
              << ". Delta V needed at Venus: " << deltaVVectorMessenger[ 4 ] << std::endl;
    std::cout << " Time of 3d DSM: " << timeVectorMessenger[ 5 ]
              << ". Delta V needed for 3d DSM: " << deltaVVectorMessenger[ 5 ] << std::endl;
    std::cout << " Time of second Venus visit: " << timeVectorMessenger[ 6 ]
              << ". Delta V needed at Venus: " << deltaVVectorMessenger[ 6 ] << std::endl;
    std::cout << " Time of 4th DSM: " << timeVectorMessenger[ 7 ]
              << ". Delta V needed for 4th DSM: " << deltaVVectorMessenger[ 7 ] << std::endl;
    std::cout << " Time of Mercury capture: " << timeVectorMessenger[ 8 ]
              << ". Delta V needed at Mercury: " << deltaVVectorMessenger[ 8 ] << std::endl;

    // Define vectors to calculate intermediate points
    std::vector< Eigen::Vector3d > interPositionVectorMessenger;
    std::vector< double > interTimeVectorMessenger;

    // Calculate intermediate points and write to file
    std::string outputFileTraj = tudat_applications::getOutputPath( ) + "messengerTrajectory.dat";
    Messenger.intermediatePoints( 1000.0 , interPositionVectorMessenger, interTimeVectorMessenger );
    writeTrajectoryToFile( interPositionVectorMessenger, interTimeVectorMessenger, outputFileTraj );

    // Define vectors to calculate intermediate points
    std::vector< Eigen::Vector3d > manPositionVectorMessenger;
    std::vector< double > manTimeVectorMessenger;
    std::vector< double > manDeltaVVectorMessenger;

    // Calculate maneuvers and write to file
    std::string outputFileMan = tudat_applications::getOutputPath( ) + "messengerManeuvers.dat";
    Messenger.maneuvers( manPositionVectorMessenger, manTimeVectorMessenger, manDeltaVVectorMessenger );
    writeTrajectoryToFile( manPositionVectorMessenger, manTimeVectorMessenger, outputFileMan );

    // Calculate trajectories of the planets and output to file
    std::vector< Eigen::Vector3d > positionVectorEarth;
    std::vector< double > timeVectorEarth;
    std::vector< Eigen::Vector3d > positionVectorVenus;
    std::vector< double > timeVectorVenus;
    std::vector< Eigen::Vector3d > positionVectorMercury;
    std::vector< double > timeVectorMercury;

    // Earth
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorEarth,
                timeVectorEarth );

    // Venus
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorVenus,
                timeVectorVenus );

    // Mercury
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorMercury,
                timeVectorMercury );

    std::string outputFilePlanetE = tudat_applications::getOutputPath(  ) + "earthTrajectory.dat";
    writeTrajectoryToFile( positionVectorEarth, timeVectorEarth, outputFilePlanetE );

    std::string outputFilePlanetV = tudat_applications::getOutputPath(  ) + "venusTrajectory.dat";
    writeTrajectoryToFile( positionVectorVenus, timeVectorVenus, outputFilePlanetV );

    std::string outputFilePlanetM = tudat_applications::getOutputPath(  ) + "mercuryTrajectory.dat";
    writeTrajectoryToFile( positionVectorMercury, timeVectorMercury, outputFilePlanetM );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
