/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include "Tudat/SimulationSetup/PropagationSetup/fullPropagationLambertTargeter.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "Tudat/Basics/testMacros.h"


// required for the MGA
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/planetTrajectory.h"


namespace tudat
{

namespace unit_tests
{

using namespace tudat;

//! Test if the difference between the Lambert targeter solution and the full dynamics problem is computed correctly.
BOOST_AUTO_TEST_CASE( testFullPropagationLambertTargeter )
{

    std::cout.precision(20);

    double initialTime = 0.0;
    double fixedStepSize = 1.0;

    Eigen::Vector3d cartesianPositionAtDeparture ( 2.0 * 6.378136e6, 0.0, 0.0 );
    Eigen::Vector3d cartesianPositionAtArrival ( 2.0 * 6.378136e6, 2.0 * std::sqrt( 3.0 ) * 6.378136e6, 0.0 );

    double timeOfFlight = 806.78 * 5.0;

    std::vector< std::string > bodiesToPropagate; bodiesToPropagate.push_back("spacecraft");
    std::vector< std::string > centralBodies; centralBodies.push_back("Earth");

    // Define integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
                ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize);

    // Define the body map.
    simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapLambertTargeter("Earth", "spacecraft");
    basic_astrodynamics::AccelerationMap accelerationModelMap = propagators::setupAccelerationMapLambertTargeter(
                "Earth", "spacecraft", bodyMap);

    std::vector< std::string > departureAndArrivalBodies;
    departureAndArrivalBodies.push_back("Earth");
    departureAndArrivalBodies.push_back("Mars");

   // Compute the difference in state between the full problem and the Lambert targeter solution at departure and at arrival
    std::pair< Eigen::Vector6d, Eigen::Vector6d > differenceState =
            propagators::getDifferenceFullPropagationWrtLambertTargeterAtDepartureAndArrival(cartesianPositionAtDeparture,
             cartesianPositionAtArrival, timeOfFlight, bodyMap, accelerationModelMap, bodiesToPropagate,
             centralBodies, integratorSettings, departureAndArrivalBodies);

    Eigen::Vector6d differenceStateAtDeparture = differenceState.first;
    Eigen::Vector6d differenceStateAtArrival = differenceState.second;

    std::cout << "differenceStateAtDeparture: " << differenceStateAtDeparture << "\n\n";
    std::cout << "differenceStateAtArrival: " << differenceStateAtArrival << "\n\n";


    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( differenceStateAtDeparture( i ) ), 1.0 );
        BOOST_CHECK_SMALL( std::fabs( differenceStateAtDeparture( i + 3 ) ), 1.0E-6 );
        BOOST_CHECK_SMALL( std::fabs( differenceStateAtArrival( i ) ), 1.0 );
        BOOST_CHECK_SMALL( std::fabs( differenceStateAtArrival( i + 3 ) ), 1.0E-6 );
    }




    //// TEST MGA TRAJECTORY


    int test = 2;

    if (test == 1){

    // Specify required parameters
    // Specify the number of legs and type of legs.
    int numberOfLegs = 6;
    std::vector< int > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[ 0 ] = transfer_trajectories::mga_Departure;
    legTypeVector[ 1 ] = transfer_trajectories::mga_Swingby;
    legTypeVector[ 2 ] = transfer_trajectories::mga_Swingby;
    legTypeVector[ 3 ] = transfer_trajectories::mga_Swingby;
    legTypeVector[ 4 ] = transfer_trajectories::mga_Swingby;
    legTypeVector[ 5 ] = transfer_trajectories::capture;

    // Name of the bodies involved in the trajectory
    std::vector< std::string > nameBodiesTrajectory;
    nameBodiesTrajectory.push_back("Earth");
    nameBodiesTrajectory.push_back("Venus");
    nameBodiesTrajectory.push_back("Venus");
    nameBodiesTrajectory.push_back("Earth");
    nameBodiesTrajectory.push_back("Jupiter");
    nameBodiesTrajectory.push_back("Saturn");

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


    std::vector< std::string > centralBody;
    centralBody.push_back( "Sun" );
    std::vector< std::string > bodyToPropagate;
    bodyToPropagate.push_back( "spacecraft" );
    std::map< double, Eigen::Vector6d > outputTest = propagators::fullPropagationMGA(numberOfLegs, nameBodiesTrajectory,
                                                                                     centralBody, bodyToPropagate,
                                                                                   legTypeVector, ephemerisVector,
                                                                                   gravitationalParameterVector,
                                                                                   variableVector, sunGravitationalParameter,
                                                                                   minimumPericenterRadii, semiMajorAxes,
                                                                                   eccentricities, integratorSettings);

    std::cout << "approximate position Earth: " << ephemerisVector[1]->getCartesianState( (-789.8117 + 158.302027105278) * physical_constants::JULIAN_DAY) << "\n\n";


    }

    else if (test == 2){

    // Specify required parameters
    // Specify the number of legs and type of legs.
    int numberOfLegs = 5;
    std::vector< int > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[ 0 ] = transfer_trajectories::mga1DsmVelocity_Departure;
    legTypeVector[ 1 ] = transfer_trajectories::mga1DsmVelocity_Swingby;
    legTypeVector[ 2 ] = transfer_trajectories::mga1DsmVelocity_Swingby;
    legTypeVector[ 3 ] = transfer_trajectories::mga1DsmVelocity_Swingby;
    legTypeVector[ 4 ] = transfer_trajectories::capture;

    // Name of the bodies involved in the trajectory
    std::vector< std::string > nameBodiesTrajectory;
    nameBodiesTrajectory.push_back("Earth");
    nameBodiesTrajectory.push_back("DSM1");
    nameBodiesTrajectory.push_back("Earth");
    nameBodiesTrajectory.push_back("DSM2");
    nameBodiesTrajectory.push_back("Venus");
    nameBodiesTrajectory.push_back("DSM3");
    nameBodiesTrajectory.push_back("Venus");
    nameBodiesTrajectory.push_back("DSM4");
    nameBodiesTrajectory.push_back("Mercury");

    // Create the ephemeris vector.
    std::vector< ephemerides::EphemerisPointer > ephemerisVector( numberOfLegs );
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
    Eigen::VectorXd gravitationalParameterVector( numberOfLegs );
    gravitationalParameterVector << 3.9860119e14, 3.9860119e14, 3.24860e14, 3.24860e14, 2.2321e13;

    // Create variable vector.
    Eigen::VectorXd variableVector;
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

    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;

    // Create minimum pericenter radii vector
    Eigen::VectorXd minimumPericenterRadii( numberOfLegs );
    minimumPericenterRadii << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN;

    // Create departure and capture variables.
    Eigen::VectorXd semiMajorAxes( 2 ), eccentricities( 2 );
    semiMajorAxes << std::numeric_limits< double >::infinity( ),
            std::numeric_limits< double >::infinity( );
    eccentricities << 0.0, 0.0;


    std::vector< std::string > centralBody;
    centralBody.push_back( "Sun" );
    std::vector< std::string > bodyToPropagate;
    bodyToPropagate.push_back( "spacecraft" );
    std::map< double, Eigen::Vector6d > outputTest = propagators::fullPropagationMGA(numberOfLegs, nameBodiesTrajectory,
                                                                                     centralBody, bodyToPropagate,
                                                                                   legTypeVector, ephemerisVector,
                                                                                   gravitationalParameterVector,
                                                                                   variableVector, sunGravitationalParameter,
                                                                                   minimumPericenterRadii, semiMajorAxes,
                                                                                   eccentricities, integratorSettings);




    std::cout << "approximate position Earth: " << ephemerisVector[4]->getCartesianState( 1.9265e+08) << "\n\n";

    }



}

}

}


