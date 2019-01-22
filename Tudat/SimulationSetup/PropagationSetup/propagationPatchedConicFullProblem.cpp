/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/propagationPatchedConicFullProblem.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationLambertTargeterFullProblem.h"

namespace tudat
{

namespace propagators
{


transfer_trajectories::Trajectory createTransferTrajectoryObject(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::vector< transfer_trajectories::TransferLegType >& transferLegTypes,
        const std::vector< double > trajectoryIndependentVariables,
        const std::vector< double > minimumPericenterRadii,
        const bool includeDepartureDeltaV,
        const double departureSemiMajorAxis,
        const double departureEccentricity,
        const bool includeArrivalDeltaV,
        const double arrivalSemiMajorAxis,
        const double arrivalEccentricity )
{
    int numberOfLegs = transferBodyOrder.size( );

    std::vector< ephemerides::EphemerisPointer > ephemerisVector;
    Eigen::VectorXd gravitationalParameterVector = Eigen::VectorXd::Zero( transferBodyOrder.size( ) );

    for( int i = 0; i < transferBodyOrder.size( ); i++ )
    {
        if( bodyMap.count( transferBodyOrder.at( i ) ) != 0 )
        {
            ephemerisVector.push_back( bodyMap.at( transferBodyOrder.at( i ) )->getEphemeris( ) );
            gravitationalParameterVector( i ) =
                    bodyMap.at( transferBodyOrder.at( i ) )->getGravityFieldModel( )->getGravitationalParameter( );
        }
    }

    double centralBodyGravitationalParameter;
    if( bodyMap.count( centralBody ) != 0 )
    {
        centralBodyGravitationalParameter = bodyMap.at( centralBody )->getGravityFieldModel( )->getGravitationalParameter( );
    }

    Eigen::VectorXd semiMajorAxesVector =
            ( Eigen::VectorXd( 2 ) << departureSemiMajorAxis, arrivalSemiMajorAxis ).finished( );
    Eigen::VectorXd eccentricityVector =
            ( Eigen::VectorXd( 2 ) << departureEccentricity, arrivalEccentricity ).finished( );

    return transfer_trajectories::Trajectory(
                numberOfLegs, transferLegTypes, ephemerisVector, gravitationalParameterVector,
                utilities::convertStlVectorToEigenVector( trajectoryIndependentVariables ), centralBodyGravitationalParameter,
                 utilities::convertStlVectorToEigenVector( minimumPericenterRadii ),
                semiMajorAxesVector, eccentricityVector, includeDepartureDeltaV, includeArrivalDeltaV );

}



void fullPropagationMGA(
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const std::vector< double >& trajectoryVariableVector,
        const std::vector< double >& minimumPericenterRadiiVector,
        const std::vector< double >& semiMajorAxesVector,
        const std::vector< double >& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg){


    // calculate the patched conic trajectory from the body map
    transfer_trajectories::Trajectory trajectory = propagators::createTransferTrajectoryObject(
            bodyMap, transferBodyOrder, centralBody, legTypeVector, trajectoryVariableVector, minimumPericenterRadiiVector, true,
            semiMajorAxesVector[0], eccentricitiesVector[0], true, semiMajorAxesVector[1], eccentricitiesVector[1]);


    int numberOfLegs = legTypeVector.size();
    int numberLegsIncludingDSM = ((trajectoryVariableVector.size() - 1 - numberOfLegs) / 4.0) + numberOfLegs ;

    std::vector< Eigen::Vector3d > positionVector;
    std::vector< double > timeVector;
    std::vector< double > deltaVVector;
    double totalDeltaV;


    std::map< double, std::pair<Eigen::Vector6d, Eigen::Vector6d> > stateDifferenceAtDepartureAndArrival;

    // Calculate the trajectory
    trajectory.calculateTrajectory( totalDeltaV );
    trajectory.maneuvers( positionVector, timeVector, deltaVVector );

    std::map< int, Eigen::Vector3d > cartesianPositionAtDepartureLambertTargeter;
    std::map< int, Eigen::Vector3d > cartesianPositionAtArrivalLambertTargeter;


    double timeOfFlight;
    int counterLegTotal = 0;
    int counterLegWithDSM = 0;

    double initialTime = timeVector[0];
    double totalTime = initialTime;
    std::vector< double > absoluteTimeVector; absoluteTimeVector.push_back(initialTime);
    std::vector< double > timeOfFlightVector;



    // Include manoeuvres between transfer bodies when required
    std::vector< std::string > bodiesAndManoeuvresOrder;

    int counterDSMs = 1;
    for (int i = 0 ; i < numberOfLegs ; i++){

        bodiesAndManoeuvresOrder.push_back(transferBodyOrder[i]);

        if (legTypeVector[i] != transfer_trajectories::mga_Departure &&
                legTypeVector[i] != transfer_trajectories::mga_Swingby){

            bodiesAndManoeuvresOrder.push_back("DSM" + std::to_string(counterDSMs));
            counterDSMs++;

        }

    }





    // Calculate the time of flight for each leg (considering that a deep space manoeuvre divides a leg into two smaller ones)
    for (int i = 0 ; i < numberOfLegs - 1 ; i ++){

        if (legTypeVector[i] == transfer_trajectories::mga_Departure ||
                legTypeVector[i] == transfer_trajectories::mga_Swingby ){

            timeOfFlight = trajectoryVariableVector[1 + counterLegTotal];
            timeOfFlightVector.push_back( timeOfFlight );

            totalTime += timeOfFlight;
            absoluteTimeVector.push_back(totalTime);

            counterLegTotal++;

        }

        else {

            // first part of the leg (from departure to DSM)
            timeOfFlight = trajectoryVariableVector[numberOfLegs + 1 + counterLegWithDSM * 4]
                    * trajectoryVariableVector[counterLegTotal + 1];
            timeOfFlightVector.push_back( timeOfFlight );

            totalTime += timeOfFlight;
            absoluteTimeVector.push_back(totalTime);

            // second part of the leg (from DSM to arrival)
            timeOfFlight = (1 - trajectoryVariableVector[numberOfLegs + 1 + counterLegWithDSM * 4])
                    * trajectoryVariableVector[counterLegTotal + 1];
            timeOfFlightVector.push_back( timeOfFlight );

            totalTime += timeOfFlight;
            absoluteTimeVector.push_back(totalTime);



            counterLegWithDSM++;
            counterLegTotal++;
        }

    }

    std::string testString = "DSM" + std::to_string(1);
    std::cout << "testString: " << testString << "\n\n";


    for (int i = 0; i<numberLegsIncludingDSM-1 ; i++)
    {

        cartesianPositionAtDepartureLambertTargeter[ i ] = positionVector[i];
        cartesianPositionAtArrivalLambertTargeter[ i ] = positionVector[i+1];

        std::vector< std::string > departureAndArrivalBodies;
        departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[i] );
        departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[1 + i]);
        std::cout << "departure and arrival bodies: " << departureAndArrivalBodies[0] << "\n\n";
        std::cout << "departure and arrival bodies: " << departureAndArrivalBodies[1] << "\n\n";


        Eigen::Vector3d cartesianPositionAtDeparture = cartesianPositionAtDepartureLambertTargeter[i];
        Eigen::Vector3d cartesianPositionAtArrival = cartesianPositionAtArrivalLambertTargeter[i];



       // Compute the difference in state between the full problem and the Lambert targeter solution at departure and at arrival
        std::map< double, Eigen::Vector6d > lambertTargeterResultForOneLeg;
        std::map< double, Eigen::Vector6d > fullProblemResultForOneLeg;

        integratorSettings->initialTime_ = initialTime;

        propagateLambertTargeterAndFullProblem( cartesianPositionAtDeparture, cartesianPositionAtArrival,
                timeOfFlightVector[i], initialTime, bodyMap, accelerationMap, bodyToPropagate, centralBody,
                integratorSettings, lambertTargeterResultForOneLeg, fullProblemResultForOneLeg,
                departureAndArrivalBodies, false, false);


        lambertTargeterResultForEachLeg[i] = lambertTargeterResultForOneLeg;
        fullProblemResultForEachLeg[i] = fullProblemResultForOneLeg;


    }


}




std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullPropagationWrtLambertTargeterMGA(simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const std::vector<double>& trajectoryVariableVector,
        const std::vector<double>& minimumPericenterRadiiVector,
        const std::vector<double>& semiMajorAxesVector,
        const std::vector<double>& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings)
{

    int numberOfLegs = legTypeVector.size();
    int numberLegsIncludingDSM = ((trajectoryVariableVector.size()-1-numberOfLegs)/4.0) + numberOfLegs ;


    std::map< int, std::map< double, Eigen::Vector6d > > lambertTargeterResultForEachLeg;
    std::map< int, std::map< double, Eigen::Vector6d > > fullProblemResultForEachLeg;



    fullPropagationMGA(bodyMap, accelerationMap, transferBodyOrder, centralBody, bodyToPropagate,
                         legTypeVector, trajectoryVariableVector, minimumPericenterRadiiVector, semiMajorAxesVector, eccentricitiesVector,
                         integratorSettings, lambertTargeterResultForEachLeg, fullProblemResultForEachLeg);


    std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > stateDifferenceAtArrivalAndDepartureForEachLeg;

    for (int i = 0 ; i< numberLegsIncludingDSM-1 ; i++){

        std::map< double, Eigen::Vector6d > lambertTargeterResultCurrentLeg = lambertTargeterResultForEachLeg[i];
        std::map< double, Eigen::Vector6d > fullProblemResultCurrentLeg = fullProblemResultForEachLeg[i];

        Eigen::Vector6d stateLambertTargeterAtDepartureForOneLeg = lambertTargeterResultCurrentLeg.begin( )->second;
        Eigen::Vector6d stateFullProblemAtDepartureForOneLeg = fullProblemResultCurrentLeg.begin( )->second;
        Eigen::Vector6d stateLambertTargeterAtArrivalForOneLeg = lambertTargeterResultCurrentLeg.rbegin( )->second;
        Eigen::Vector6d stateFullProblemAtArrivalForOneLeg = fullProblemResultCurrentLeg.rbegin( )->second;

        stateDifferenceAtArrivalAndDepartureForEachLeg[i] = std::make_pair( stateLambertTargeterAtDepartureForOneLeg -
                                                                            stateFullProblemAtDepartureForOneLeg,
                                                                            stateLambertTargeterAtArrivalForOneLeg -
                                                                            stateFullProblemAtArrivalForOneLeg);



    }


    return stateDifferenceAtArrivalAndDepartureForEachLeg;

}




}

}


