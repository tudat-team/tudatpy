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
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"

#include <Tudat/InputOutput/basicInputOutput.h>
#include <tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/applicationOutput.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

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

namespace propagators
{


//! Function to setup a body map corresponding to the assumptions of a patched conics trajectory,
//! using default ephemerides for the central and transfer bodies.
simulation_setup::NamedBodyMap setupBodyMapFromEphemeridesForPatchedConicsTrajectory(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& nameTransferBodies)
{
    spice_interface::loadStandardSpiceKernels( );

    // Create central and transfer bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( nameCentralBody );
    for ( unsigned int i = 0 ; i < nameTransferBodies.size() ; i++ ){
        bodiesToCreate.push_back( nameTransferBodies[i] );
    }


    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ nameCentralBody ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ nameCentralBody ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ nameCentralBody ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );


    // Define body to propagate.
    bodyMap[ nameBodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ nameBodyToPropagate ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                      std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                      < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );


    return bodyMap;
}



//! Function to setup a body map corresponding to the assumptions of the patched conics trajectory,
//! the ephemerides of the transfer bodies being provided as inputs.
simulation_setup::NamedBodyMap setupBodyMapFromUserDefinedEphemeridesForPatchedConicsTrajectory(const std::string& nameCentralBody,
                                                                                                const std::string& nameBodyToPropagate,
                                                                                                const std::vector< std::string >& nameTransferBodies,
                                                                                                const std::vector< ephemerides::EphemerisPointer >& ephemerisVectorTransferBodies,
                                                                                                const std::vector< double >& gravitationalParametersTransferBodies)
{

    spice_interface::loadStandardSpiceKernels( );


    // Create central body object.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( nameCentralBody );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "J2000";

    // Define central body ephemeris settings.
    bodySettings[ nameCentralBody ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ nameCentralBody ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ nameCentralBody ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ nameBodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ nameBodyToPropagate ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                      std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                      < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    // Define ephemeris and gravity field for the transfer bodies.
    for ( unsigned int i = 0 ; i < nameTransferBodies.size() ; i++){

        bodyMap[ nameTransferBodies[i] ] = std::make_shared< simulation_setup::Body >( );
        bodyMap[ nameTransferBodies[i] ]->setEphemeris( ephemerisVectorTransferBodies[i] );
        bodyMap[ nameTransferBodies[i] ]->setGravityFieldModel( simulation_setup::createGravityFieldModel(
                                                                    std::make_shared< simulation_setup::CentralGravityFieldSettings >( gravitationalParametersTransferBodies[i] ),
                                                                    nameTransferBodies[i] ) );
    }

    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    return bodyMap;

}



//! Function to directly setup a vector of acceleration maps for a patched conics trajectory.
std::vector < basic_astrodynamics::AccelerationMap > setupAccelerationMapPatchedConicsTrajectory(
        const double numberOfLegs,
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const simulation_setup::NamedBodyMap& bodyMap )
{

    std::vector< basic_astrodynamics::AccelerationMap > accelerationMapsVector;

    for (int i = 0 ; i < numberOfLegs ; i++){

        accelerationMapsVector.push_back( setupAccelerationMapLambertTargeter(nameCentralBody, nameBodyToPropagate, bodyMap) );

    }

    return accelerationMapsVector;

}




//! Function to create the trajectory from the body map.
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

    for( unsigned int i = 0; i < transferBodyOrder.size( ); i++ )
    {
        if( bodyMap.count( transferBodyOrder.at( i ) ) != 0 )
        {
            ephemerisVector.push_back( bodyMap.at( transferBodyOrder.at( i ) )->getEphemeris( ) );
            gravitationalParameterVector( i ) =
                    bodyMap.at( transferBodyOrder.at( i ) )->getGravityFieldModel( )->getGravitationalParameter( );

        }
    }

    double centralBodyGravitationalParameter = TUDAT_NAN;
    if( bodyMap.count( centralBody ) != 0 )
    {
        centralBodyGravitationalParameter = bodyMap.at( centralBody )->getGravityFieldModel( )->getGravitationalParameter( );
    }
    else
    {
        throw std::runtime_error( "Error, central body " + centralBody + " not found when creating transfer trajectory object" );
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



//! Function to calculate the patched conics trajectory and to propagate the corresponding full problem.
void fullPropagationPatchedConicsTrajectory(
        simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const std::vector< double >& trajectoryVariableVector,
        const std::vector< double >& minimumPericenterRadiiVector,
        const std::vector< double >& semiMajorAxesVector,
        const std::vector< double >& eccentricitiesVector,
        const std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
        std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg)
{
    int numberOfLegs = legTypeVector.size();
    int numberLegsIncludingDSM = ((trajectoryVariableVector.size() - 1 - numberOfLegs) / 4.0) + numberOfLegs ;

    // Define the patched conic trajectory from the body map.
    transfer_trajectories::Trajectory trajectory = propagators::createTransferTrajectoryObject(
                bodyMap, transferBodyOrder, centralBody, legTypeVector, trajectoryVariableVector, minimumPericenterRadiiVector, true,
                semiMajorAxesVector[0], eccentricitiesVector[0], true, semiMajorAxesVector[1], eccentricitiesVector[1]);

    // Clear output maps.
    lambertTargeterResultForEachLeg.clear();
    fullProblemResultForEachLeg.clear();

    // Calculate the trajectory.
    std::vector< Eigen::Vector3d > positionVector;
    std::vector< double > timeVector;
    std::vector< double > deltaVVector;
    double totalDeltaV;

    trajectory.calculateTrajectory( totalDeltaV );
    trajectory.maneuvers( positionVector, timeVector, deltaVVector );

    // Include manoeuvres between transfer bodies when required and associate acceleration maps
    //            (considering that a deep space manoeuvre divides a leg into two smaller ones).
    std::vector< std::string > bodiesAndManoeuvresOrder;
    std::vector< basic_astrodynamics::AccelerationMap > accelerationMapForEachLeg;

    int counterDSMs = 1;
    for (int i = 0 ; i < numberOfLegs ; i++){

        bodiesAndManoeuvresOrder.push_back(transferBodyOrder[i]);
        accelerationMapForEachLeg.push_back( accelerationMap[i] );

        if (legTypeVector[i] != transfer_trajectories::mga_Departure && legTypeVector[i] != transfer_trajectories::mga_Swingby){

            bodiesAndManoeuvresOrder.push_back("DSM" + std::to_string(counterDSMs));
            accelerationMapForEachLeg.push_back( accelerationMap[i] );

            counterDSMs++;
        }
    }



    // Compute swing-by (or DSM) time and corresponding time of flight for each leg.

    double timeOfFlight;
    int counterLegTotal = 0;
    int counterLegWithDSM = 0;
    double initialTime = timeVector[0];
    double totalTime = initialTime;

    std::vector< double > absoluteTimeVector; absoluteTimeVector.push_back(initialTime);
    std::vector< double > timeOfFlightVector;

    for (int i = 0 ; i < numberOfLegs - 1 ; i ++){

        if (legTypeVector[i] == transfer_trajectories::mga_Departure || legTypeVector[i] == transfer_trajectories::mga_Swingby ){

            timeOfFlight = trajectoryVariableVector[1 + counterLegTotal];
            timeOfFlightVector.push_back( timeOfFlight );

            totalTime += timeOfFlight;
            absoluteTimeVector.push_back(totalTime);

            counterLegTotal++;
        }


        else { // If a DSM is performed in the current leg

            // First part of the leg (from departure to DSM)
            timeOfFlight = trajectoryVariableVector[numberOfLegs + 1 + counterLegWithDSM * 4]
                    * trajectoryVariableVector[counterLegTotal + 1];
            timeOfFlightVector.push_back( timeOfFlight );

            totalTime += timeOfFlight;
            absoluteTimeVector.push_back(totalTime);

            // Second part of the leg (from DSM to arrival)
            timeOfFlight = (1 - trajectoryVariableVector[numberOfLegs + 1 + counterLegWithDSM * 4])
                    * trajectoryVariableVector[counterLegTotal + 1];
            timeOfFlightVector.push_back( timeOfFlight );

            totalTime += timeOfFlight;
            absoluteTimeVector.push_back(totalTime);


            counterLegWithDSM++;
            counterLegTotal++;
        }
    }




    for (int i = 0; i < numberLegsIncludingDSM - 1 ; i++)
    {
        std::vector< std::string > departureAndArrivalBodies;
        departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[i] );
        departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[1 + i]);

        // Define cartesian position at departure and arrival for the current leg.
        Eigen::Vector3d cartesianPositionAtDepartureCurrentLeg = positionVector[i]; //cartesianPositionAtDepartureForEachLeg[i];
        Eigen::Vector3d cartesianPositionAtArrivalCurrentLeg = positionVector[i+1]; //cartesianPositionAtArrivalForEachLeg[i];

        // Compute the difference in state between the full problem and the Lambert targeter solution for the current leg.
        std::map< double, Eigen::Vector6d > lambertTargeterResultCurrentLeg;
        std::map< double, Eigen::Vector6d > fullProblemResultCurrentLeg;

        integratorSettings->initialTime_ = absoluteTimeVector[i];

        propagateLambertTargeterAndFullProblem(
                    timeOfFlightVector[i], absoluteTimeVector[i], bodyMap, accelerationMapForEachLeg[i], bodyToPropagate, centralBody,
                    terminationSettings, integratorSettings, lambertTargeterResultCurrentLeg, fullProblemResultCurrentLeg,
                    departureAndArrivalBodies, bodyMap[ centralBody ]->getGravityFieldModel()->getGravitationalParameter(),
                    cartesianPositionAtDepartureCurrentLeg, cartesianPositionAtArrivalCurrentLeg);

        lambertTargeterResultForEachLeg[i] = lambertTargeterResultCurrentLeg;
        fullProblemResultForEachLeg[i] = fullProblemResultCurrentLeg;

    }
}



////! Function to calculate the patched conics trajectory and to propagate the corresponding full problem.
//void fullPropagationPatchedConicsTrajectory(
//        simulation_setup::NamedBodyMap& bodyMap,
//        const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
//        const std::vector< std::string >& transferBodyOrder,
//        const std::string& centralBody,
//        const std::string& bodyToPropagate,
//        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
//        const std::vector< double >& trajectoryVariableVector,
//        const std::vector< double >& minimumPericenterRadiiVector,
//        const std::vector< double >& semiMajorAxesVector,
//        const std::vector< double >& eccentricitiesVector,
//        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
//        const bool terminationSphereOfInfluence,
//        std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg,
//        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg)
//{

//    int numberOfLegs = legTypeVector.size();
//    int numberLegsIncludingDSM = ((trajectoryVariableVector.size() - 1 - numberOfLegs) / 4.0) + numberOfLegs ;

//    // Define the patched conic trajectory from the body map.
//    transfer_trajectories::Trajectory trajectory = propagators::createTransferTrajectoryObject(
//                bodyMap, transferBodyOrder, centralBody, legTypeVector, trajectoryVariableVector, minimumPericenterRadiiVector, true,
//                semiMajorAxesVector[0], eccentricitiesVector[0], true, semiMajorAxesVector[1], eccentricitiesVector[1]);



//    // Clear output maps.
//    lambertTargeterResultForEachLeg.clear();
//    fullProblemResultForEachLeg.clear();

//    // Calculate the trajectory.
//    std::vector< Eigen::Vector3d > positionVector;
//    std::vector< double > timeVector;
//    std::vector< double > deltaVVector;
//    double totalDeltaV;

//    trajectory.calculateTrajectory( totalDeltaV );
//    trajectory.maneuvers( positionVector, timeVector, deltaVVector );

//    std::cout << "resulting Delta V Messenger: " << totalDeltaV << "\n\n";
//    for (int i = 0 ; i < numberLegsIncludingDSM ; i++){
//        std::cout << "leg " << i << ": " << "\n\n";
//        std::cout << "positionVectorMessenger: " << positionVector[i] << "\n\n";
//        std::cout << "timeVectorMessenger: " << timeVector[i] << "\n\n";
//        std::cout << "deltaVVectorMessenger: " << deltaVVector[i] << "\n\n";
//    }


//    // Include manoeuvres between transfer bodies when required and associate acceleration maps
//    //            (considering that a deep space manoeuvre divides a leg into two smaller ones).
//    std::vector< std::string > bodiesAndManoeuvresOrder;
//    std::vector< basic_astrodynamics::AccelerationMap > accelerationMapForEachLeg;

//    int counterDSMs = 1;
//    for (int i = 0 ; i < numberOfLegs ; i++){

//        bodiesAndManoeuvresOrder.push_back(transferBodyOrder[i]);
//        accelerationMapForEachLeg.push_back( accelerationMap[i] );

//        if (legTypeVector[i] != transfer_trajectories::mga_Departure && legTypeVector[i] != transfer_trajectories::mga_Swingby){

//            bodiesAndManoeuvresOrder.push_back("DSM" + std::to_string(counterDSMs));
//            accelerationMapForEachLeg.push_back( accelerationMap[i] );

//            counterDSMs++;
//        }
//    }



//    // Compute swing-by (or DSM) time and corresponding time of flight for each leg.

//    double timeOfFlight;
//    int counterLegTotal = 0;
//    int counterLegWithDSM = 0;
//    double initialTime = timeVector[0];
//    double totalTime = initialTime;

//    std::vector< double > absoluteTimeVector; absoluteTimeVector.push_back(initialTime);
//    std::vector< double > timeOfFlightVector;

//    for (int i = 0 ; i < numberOfLegs - 1 ; i ++){

//        if (legTypeVector[i] == transfer_trajectories::mga_Departure || legTypeVector[i] == transfer_trajectories::mga_Swingby ){

//            timeOfFlight = trajectoryVariableVector[1 + counterLegTotal];
//            timeOfFlightVector.push_back( timeOfFlight );

//            totalTime += timeOfFlight;
//            absoluteTimeVector.push_back(totalTime);

//            counterLegTotal++;
//        }


//        else { // If a DSM is performed in the current leg

//            // First part of the leg (from departure to DSM)
//            timeOfFlight = trajectoryVariableVector[numberOfLegs + 1 + counterLegWithDSM * 4]
//                    * trajectoryVariableVector[counterLegTotal + 1];
//            timeOfFlightVector.push_back( timeOfFlight );

//            totalTime += timeOfFlight;
//            absoluteTimeVector.push_back(totalTime);

//            // Second part of the leg (from DSM to arrival)
//            timeOfFlight = (1 - trajectoryVariableVector[numberOfLegs + 1 + counterLegWithDSM * 4])
//                    * trajectoryVariableVector[counterLegTotal + 1];
//            timeOfFlightVector.push_back( timeOfFlight );

//            totalTime += timeOfFlight;
//            absoluteTimeVector.push_back(totalTime);


//            counterLegWithDSM++;
//            counterLegTotal++;
//        }
//    }




//    for (int i = 0; i < numberLegsIncludingDSM - 1 ; i++)
//    {
//        std::vector< std::string > departureAndArrivalBodies;
//        departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[i] );
//        departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[1 + i]);

//        // Define cartesian position at departure and arrival for the current leg.
//        Eigen::Vector3d cartesianPositionAtDepartureCurrentLeg = positionVector[i]; //cartesianPositionAtDepartureForEachLeg[i];
//        Eigen::Vector3d cartesianPositionAtArrivalCurrentLeg = positionVector[i+1]; //cartesianPositionAtArrivalForEachLeg[i];

//        // Compute the difference in state between the full problem and the Lambert targeter solution for the current leg.
//        std::map< double, Eigen::Vector6d > lambertTargeterResultCurrentLeg;
//        std::map< double, Eigen::Vector6d > fullProblemResultCurrentLeg;

//        integratorSettings->initialTime_ = absoluteTimeVector[i];

//        timeOfFlightVector[i] = timeVector[i+1] - timeVector[i];
//        absoluteTimeVector[i] = timeVector[i];

//        std::cout << "LEG " << i << ": " << "\n\n";
//        std::cout << "cartesian position at departure: " << cartesianPositionAtDepartureCurrentLeg << "\n\n";
//        std::cout << "cartesian position at arrival: " << cartesianPositionAtArrivalCurrentLeg << "\n\n";
//        std::cout << "time of flight: " << timeOfFlightVector[i] << "\n\n";
//        std::cout << "absolute time vector: " << absoluteTimeVector[i] << "\n\n";

//        propagateLambertTargeterAndFullProblem(
//            timeOfFlightVector[i], absoluteTimeVector[i], bodyMap, accelerationMapForEachLeg[i], bodyToPropagate, centralBody,
//            integratorSettings, lambertTargeterResultCurrentLeg, fullProblemResultCurrentLeg, departureAndArrivalBodies,
//            terminationSphereOfInfluence, cartesianPositionAtDepartureCurrentLeg, cartesianPositionAtArrivalCurrentLeg, TUDAT_NAN,
//            TUDAT_NAN, TUDAT_NAN);

//        lambertTargeterResultForEachLeg[i] = lambertTargeterResultCurrentLeg;
//        fullProblemResultForEachLeg[i] = fullProblemResultCurrentLeg;

//    }
//}



//! Function to calculate the patched conics trajectory and to propagate the corresponding full problem.
void fullPropagationPatchedConicsTrajectoryNewVersion(
        simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const std::vector< double >& trajectoryVariableVector,
        const std::vector< double >& minimumPericenterRadiiVector,
        const std::vector< double >& semiMajorAxesVector,
        const std::vector< double >& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        const bool terminationSphereOfInfluence,
        std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg)
{

    int numberOfLegs = legTypeVector.size();
    int numberLegsIncludingDSM = ((trajectoryVariableVector.size() - 1 - numberOfLegs) / 4.0) + numberOfLegs ;

    // Define the patched conic trajectory from the body map.
    transfer_trajectories::Trajectory trajectory = propagators::createTransferTrajectoryObject(
                bodyMap, transferBodyOrder, centralBody, legTypeVector, trajectoryVariableVector, minimumPericenterRadiiVector, true,
                semiMajorAxesVector[0], eccentricitiesVector[0], true, semiMajorAxesVector[1], eccentricitiesVector[1]);



    // Clear output maps.
    lambertTargeterResultForEachLeg.clear();
    fullProblemResultForEachLeg.clear();

    // Calculate the trajectory.
    std::vector< Eigen::Vector3d > positionVector;
    std::vector< double > timeVector;
    std::vector< double > deltaVVector;
    double totalDeltaV;

    trajectory.calculateTrajectory( totalDeltaV );
    trajectory.maneuvers( positionVector, timeVector, deltaVVector );



    // Include manoeuvres between transfer bodies when required and associate acceleration maps
    //            (considering that a deep space manoeuvre divides a leg into two smaller ones).
    std::vector< std::string > bodiesAndManoeuvresOrder;
    std::vector< basic_astrodynamics::AccelerationMap > accelerationMapForEachLeg;

    int counterDSMs = 1;
    for (int i = 0 ; i < numberOfLegs ; i++){

        bodiesAndManoeuvresOrder.push_back(transferBodyOrder[i]);
        accelerationMapForEachLeg.push_back( accelerationMap[i] );

        if (legTypeVector[i] != transfer_trajectories::mga_Departure && legTypeVector[i] != transfer_trajectories::mga_Swingby){

            bodiesAndManoeuvresOrder.push_back("DSM" + std::to_string(counterDSMs));
            accelerationMapForEachLeg.push_back( accelerationMap[i] );

            counterDSMs++;
        }
    }



    int counterLegs = 0;
    int counterLegWithDSM = 0;

    double deltaV;
    Eigen::Vector3d velocityBeforeArrival;

    Eigen::Vector3d departureBodyPosition;
    Eigen::Vector3d departureBodyVelocity;
    Eigen::Vector3d velocityAfterDeparture;

    for (int i = 0 ; i < numberOfLegs - 1 ; i++){

        std::cout << "LEG " << i << ": " << "\n\n";

        if (legTypeVector[i] == transfer_trajectories::mga_Departure || legTypeVector[i] == transfer_trajectories::mga_Swingby){

            std::vector< std::string > departureAndArrivalBodies;
            departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs] );
            departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs+1]);

            // Compute the difference in state between the full problem and the Lambert targeter solution for the current leg.
            std::map< double, Eigen::Vector6d > lambertTargeterResultCurrentLeg;
            std::map< double, Eigen::Vector6d > fullProblemResultCurrentLeg;

            propagators::propagateMgaWithoutDsmAndFullProblem(bodyMap, accelerationMapForEachLeg[i], departureAndArrivalBodies, centralBody,
                         bodyToPropagate, positionVector[counterLegs], positionVector[counterLegs+1], timeVector[counterLegs], timeVector[counterLegs+1]
                         - timeVector[counterLegs], integratorSettings, terminationSphereOfInfluence, lambertTargeterResultCurrentLeg,
                         fullProblemResultCurrentLeg);

            lambertTargeterResultForEachLeg[counterLegs] = lambertTargeterResultCurrentLeg;
            fullProblemResultForEachLeg[counterLegs] = fullProblemResultCurrentLeg;

            counterLegs++;

        }

        if ( legTypeVector[i] == transfer_trajectories::mga1DsmVelocity_Departure || legTypeVector[i] == transfer_trajectories::mga1DsmVelocity_Swingby ){

        std::vector< std::string > departureAndArrivalBodies;
        departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs] );
        departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs+2]);

        // Define cartesian position at departure and arrival for the current leg.
//        Eigen::Vector3d cartesianPositionAtDepartureCurrentLeg = positionVector[counterLegs];
//        Eigen::Vector3d cartesianPositionAtArrivalCurrentLeg = positionVector[counterLegs+1];

        std::vector< double > trajectoryVariableVectorLeg;
        trajectoryVariableVectorLeg.push_back( trajectoryVariableVector[ numberOfLegs + 1 + (counterLegWithDSM * 4) ] );
        trajectoryVariableVectorLeg.push_back( trajectoryVariableVector[ numberOfLegs + 2 + (counterLegWithDSM * 4) ] );
        trajectoryVariableVectorLeg.push_back( trajectoryVariableVector[ numberOfLegs + 3 + (counterLegWithDSM * 4) ] );
        trajectoryVariableVectorLeg.push_back( trajectoryVariableVector[ numberOfLegs + 4 + (counterLegWithDSM * 4) ] );

        std::map< double, Eigen::Vector6d > patchedConicsResultFromDepartureToDsm;
        std::map< double, Eigen::Vector6d > fullProblemResultFromDepartureToDsm;
        std::map< double, Eigen::Vector6d > patchedConicsResultFromDsmToArrival;
        std::map< double, Eigen::Vector6d > fullProblemResultFromDsmToArrival;

        propagators::propagateMga1DsmVelocityAndFullProblem( bodyMap, accelerationMapForEachLeg[i], departureAndArrivalBodies,
               bodiesAndManoeuvresOrder[counterLegs+1], centralBody, bodyToPropagate, positionVector[counterLegs], positionVector[counterLegs+1],
               positionVector[counterLegs+2], timeVector[counterLegs], timeVector[counterLegs+1], timeVector[counterLegs+2],
               legTypeVector[i], trajectoryVariableVectorLeg, integratorSettings, terminationSphereOfInfluence, patchedConicsResultFromDepartureToDsm,
               fullProblemResultFromDepartureToDsm, patchedConicsResultFromDsmToArrival, fullProblemResultFromDsmToArrival, velocityBeforeArrival,
               velocityAfterDeparture, semiMajorAxesVector[0], eccentricitiesVector[0]);

        // results of the first part of the leg (from departure body to DSM)
        lambertTargeterResultForEachLeg[counterLegs] = patchedConicsResultFromDepartureToDsm;
        fullProblemResultForEachLeg[counterLegs] = fullProblemResultFromDepartureToDsm;
        counterLegs++;

        // results of the second part of the leg (from DSM to arrival body)
        lambertTargeterResultForEachLeg[counterLegs] = patchedConicsResultFromDsmToArrival;
        fullProblemResultForEachLeg[counterLegs] = fullProblemResultFromDsmToArrival;
        counterLegs++;
        counterLegWithDSM++;


//        if (legTypeVector[i] == transfer_trajectories::mga1DsmVelocity_Departure){ // && legTypeVector[i] != transfer_trajectories::mga_Swingby){

//            std::shared_ptr< transfer_trajectories::MissionLeg > missionLeg;
//            std::shared_ptr< transfer_trajectories::DepartureLegMga1DsmVelocity > departureLegMga1DsmVelocity =
//                    std::make_shared< transfer_trajectories::DepartureLegMga1DsmVelocity >(
//                         positionVector[ counterLegs ], positionVector[ counterLegs + 2], trajectoryVariableVector[ counterLegs + 1 ],
//                         bodyMap[ departureAndArrivalBodies[0] ]->getEphemeris()->getCartesianState(timeVector[counterLegs]).segment(3,3),
//                         bodyMap[ centralBody ]->getGravityFieldModel()->getGravitationalParameter(),
//                         bodyMap[ departureAndArrivalBodies[0] ]->getGravityFieldModel()->getGravitationalParameter(),
//                         semiMajorAxesVector[0], eccentricitiesVector[0],
//                         trajectoryVariableVector[ numberOfLegs + 1 + (counterLegWithDSM * 4) ],
//                         trajectoryVariableVector[ numberOfLegs + 2 + (counterLegWithDSM * 4) ],
//                         trajectoryVariableVector[ numberOfLegs + 3 + (counterLegWithDSM * 4) ],
//                         trajectoryVariableVector[ numberOfLegs + 4 + (counterLegWithDSM * 4) ], true  );

//            departureLegMga1DsmVelocity->calculateLeg( velocityBeforeArrival, deltaV );

//            departureLegMga1DsmVelocity->returnDepartureVariables( departureBodyPosition, departureBodyVelocity, velocityAfterDeparture );

//        }

//        if (legTypeVector[i] == transfer_trajectories::mga1DsmVelocity_Swingby){


//            Eigen::Vector3d velocityBeforePlanet = velocityBeforeArrival;

//            std::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet = std::make_shared< Eigen::Vector3d > ( velocityBeforePlanet );


//            std::shared_ptr< transfer_trajectories::MissionLeg > missionLeg;
//            std::shared_ptr< transfer_trajectories::SwingbyLegMga1DsmVelocity > swingbyLegMga1DsmVelocity =
//                    std::make_shared< transfer_trajectories::SwingbyLegMga1DsmVelocity >(
//                         positionVector[ counterLegs ], positionVector[ counterLegs + 2], trajectoryVariableVector[ i + 1 ],
//                         bodyMap[ departureAndArrivalBodies[0] ]->getEphemeris()->getCartesianState(timeVector[counterLegs]).segment(3,3),
//                         bodyMap[ centralBody ]->getGravityFieldModel()->getGravitationalParameter(),
//                         bodyMap[ departureAndArrivalBodies[0] ]->getGravityFieldModel()->getGravitationalParameter(),
//                         pointerToVelocityBeforePlanet,
//                         trajectoryVariableVector[ numberOfLegs + 1 + (counterLegWithDSM * 4) ],
//                         trajectoryVariableVector[ numberOfLegs + 2 + (counterLegWithDSM * 4) ],
//                         trajectoryVariableVector[ numberOfLegs + 3 + (counterLegWithDSM * 4) ],
//                         trajectoryVariableVector[ numberOfLegs + 4 + (counterLegWithDSM * 4) ] );

//            swingbyLegMga1DsmVelocity->calculateLeg( velocityBeforeArrival, deltaV );

//            swingbyLegMga1DsmVelocity->returnDepartureVariables( departureBodyPosition, departureBodyVelocity, velocityAfterDeparture );


//        }

//            // FIRST PART OF THE LEG: FROM DEPARTURE BODY TO DSM

//            std::map< double, Eigen::Vector6d > lambertTargeterResultCurrentLeg;
//            std::map< double, Eigen::Vector6d > fullProblemResultCurrentLeg;

//            integratorSettings->initialTime_ = timeVector[counterLegs];

//            propagatePatchedConicsLegAndFullProblem(
//                timeVector[counterLegs+1] - timeVector[counterLegs], timeVector[counterLegs], bodyMap, accelerationMapForEachLeg[i],
//                bodyToPropagate, centralBody, integratorSettings, lambertTargeterResultCurrentLeg, fullProblemResultCurrentLeg,
//                departureAndArrivalBodies, terminationSphereOfInfluence, cartesianPositionAtDepartureCurrentLeg,
//                cartesianPositionAtArrivalCurrentLeg, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, velocityAfterDeparture);

//            lambertTargeterResultForEachLeg[counterLegs] = lambertTargeterResultCurrentLeg;
//            fullProblemResultForEachLeg[counterLegs] = fullProblemResultCurrentLeg;

//            counterLegs++;

//            // SECOND PART OF THE LEG: FROM DSM TO ARRIVAL BODY

//            // Compute the difference in state between the full problem and the Lambert targeter solution for the current leg.
//            departureAndArrivalBodies.clear();
//            departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs] );
//            departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs+1]);

//            cartesianPositionAtDepartureCurrentLeg = positionVector[counterLegs];
//            cartesianPositionAtArrivalCurrentLeg = positionVector[counterLegs+1];

//            lambertTargeterResultCurrentLeg.clear();
//            fullProblemResultCurrentLeg.clear();

//            integratorSettings->initialTime_ = timeVector[counterLegs];


//            propagateLambertTargeterAndFullProblem(
//                timeVector[counterLegs+1] - timeVector[counterLegs], timeVector[counterLegs], bodyMap, accelerationMapForEachLeg[i],
//                bodyToPropagate, centralBody, integratorSettings, lambertTargeterResultCurrentLeg, fullProblemResultCurrentLeg,
//                departureAndArrivalBodies, terminationSphereOfInfluence, cartesianPositionAtDepartureCurrentLeg,
//                cartesianPositionAtArrivalCurrentLeg, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN);

//            lambertTargeterResultForEachLeg[counterLegs] = lambertTargeterResultCurrentLeg;
//            fullProblemResultForEachLeg[counterLegs] = fullProblemResultCurrentLeg;


//            counterLegs++;
//            counterLegWithDSM++;

        }

        if (legTypeVector[i] == transfer_trajectories::mga1DsmPosition_Departure
                        || legTypeVector[i] == transfer_trajectories::mga1DsmPosition_Swingby) {


            std::vector< std::string > departureAndArrivalBodies;
            departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs] );
            departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs+1]);

            // Define cartesian position at departure and arrival for the current leg.
            Eigen::Vector3d cartesianPositionAtDepartureCurrentLeg = positionVector[counterLegs];
            Eigen::Vector3d cartesianPositionAtArrivalCurrentLeg = positionVector[counterLegs+1];


            if (legTypeVector[i] == transfer_trajectories::mga1DsmPosition_Departure){ // && legTypeVector[i] != transfer_trajectories::mga_Swingby){

                std::shared_ptr< transfer_trajectories::MissionLeg > missionLeg;
                std::shared_ptr< transfer_trajectories::DepartureLegMga1DsmPosition > departureLegMga1DsmPosition =
                        std::make_shared< transfer_trajectories::DepartureLegMga1DsmPosition >(
                             positionVector[ counterLegs ], positionVector[ counterLegs + 2], trajectoryVariableVector[ counterLegs + 1 ],
                             bodyMap[ departureAndArrivalBodies[0] ]->getEphemeris()->getCartesianState(timeVector[counterLegs]).segment(3,3),
                             bodyMap[ centralBody ]->getGravityFieldModel()->getGravitationalParameter(),
                             bodyMap[ departureAndArrivalBodies[0] ]->getGravityFieldModel()->getGravitationalParameter(),
                             semiMajorAxesVector[0], eccentricitiesVector[0],
                             trajectoryVariableVector[ numberOfLegs + 1 + (counterLegWithDSM * 4) ],
                             trajectoryVariableVector[ numberOfLegs + 2 + (counterLegWithDSM * 4) ],
                             trajectoryVariableVector[ numberOfLegs + 3 + (counterLegWithDSM * 4) ],
                             trajectoryVariableVector[ numberOfLegs + 4 + (counterLegWithDSM * 4) ], true  );

                departureLegMga1DsmPosition->calculateLeg( velocityBeforeArrival, deltaV );

                departureLegMga1DsmPosition->returnDepartureVariables( departureBodyPosition, departureBodyVelocity, velocityAfterDeparture );

            }

            if (legTypeVector[i] == transfer_trajectories::mga1DsmPosition_Swingby){


                Eigen::Vector3d velocityBeforePlanet = velocityBeforeArrival;

                std::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet = std::make_shared< Eigen::Vector3d > ( velocityBeforePlanet );


                std::shared_ptr< transfer_trajectories::MissionLeg > missionLeg;
                std::shared_ptr< transfer_trajectories::SwingbyLegMga1DsmPosition > swingbyLegMga1DsmPosition =
                        std::make_shared< transfer_trajectories::SwingbyLegMga1DsmPosition >(
                             positionVector[ counterLegs ], positionVector[ counterLegs + 2], trajectoryVariableVector[ i + 1 ],
                             bodyMap[ departureAndArrivalBodies[0] ]->getEphemeris()->getCartesianState(timeVector[counterLegs]).segment(3,3),
                             bodyMap[ centralBody ]->getGravityFieldModel()->getGravitationalParameter(),
                             bodyMap[ departureAndArrivalBodies[0] ]->getGravityFieldModel()->getGravitationalParameter(),
                             pointerToVelocityBeforePlanet, minimumPericenterRadiiVector[i],
                             trajectoryVariableVector[ numberOfLegs + 1 + (counterLegWithDSM * 4) ],
                             trajectoryVariableVector[ numberOfLegs + 2 + (counterLegWithDSM * 4) ],
                             trajectoryVariableVector[ numberOfLegs + 3 + (counterLegWithDSM * 4) ],
                             trajectoryVariableVector[ numberOfLegs + 4 + (counterLegWithDSM * 4) ] );

                swingbyLegMga1DsmPosition->calculateLeg( velocityBeforeArrival, deltaV );

                swingbyLegMga1DsmPosition->returnDepartureVariables( departureBodyPosition, departureBodyVelocity, velocityAfterDeparture );


            }



            // FIRST PART OF THE LEG: FROM DEPARTURE BODY TO DSM

            std::map< double, Eigen::Vector6d > lambertTargeterResultCurrentLeg;
            std::map< double, Eigen::Vector6d > fullProblemResultCurrentLeg;

            integratorSettings->initialTime_ = timeVector[counterLegs];

            propagateLambertTargeterAndFullProblem(
                timeVector[counterLegs+1] - timeVector[counterLegs], timeVector[counterLegs], bodyMap, accelerationMapForEachLeg[i],
                bodyToPropagate, centralBody, integratorSettings, lambertTargeterResultCurrentLeg, fullProblemResultCurrentLeg,
                departureAndArrivalBodies, terminationSphereOfInfluence, cartesianPositionAtDepartureCurrentLeg,
                cartesianPositionAtArrivalCurrentLeg, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN);

            lambertTargeterResultForEachLeg[counterLegs] = lambertTargeterResultCurrentLeg;
            fullProblemResultForEachLeg[counterLegs] = fullProblemResultCurrentLeg;

            counterLegs++;



            // SECOND PART OF THE LEG: FROM DSM TO ARRIVAL BODY

            // Compute the difference in state between the full problem and the Lambert targeter solution for the current leg.
            departureAndArrivalBodies.clear();
            departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs] );
            departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs+1]);

            cartesianPositionAtDepartureCurrentLeg = positionVector[counterLegs];
            cartesianPositionAtArrivalCurrentLeg = positionVector[counterLegs+1];

            lambertTargeterResultCurrentLeg.clear();
            fullProblemResultCurrentLeg.clear();

            integratorSettings->initialTime_ = timeVector[counterLegs];


            propagateLambertTargeterAndFullProblem(
                timeVector[counterLegs+1] - timeVector[counterLegs], timeVector[counterLegs], bodyMap, accelerationMapForEachLeg[i],
                bodyToPropagate, centralBody, integratorSettings, lambertTargeterResultCurrentLeg, fullProblemResultCurrentLeg,
                departureAndArrivalBodies, terminationSphereOfInfluence, cartesianPositionAtDepartureCurrentLeg,
                cartesianPositionAtArrivalCurrentLeg, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN);

            lambertTargeterResultForEachLeg[counterLegs] = lambertTargeterResultCurrentLeg;
            fullProblemResultForEachLeg[counterLegs] = fullProblemResultCurrentLeg;


            counterLegs++;
            counterLegWithDSM++;

        }

    }
}

void propagateMgaWithoutDsmAndFullProblem(
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const std::vector< std::string > departureAndArrivalBodies,
        const std::string centralBody,
        const std::string bodyToPropagate,
        const Eigen::Vector3d cartesianPositionAtDeparture,
        const Eigen::Vector3d cartesianPositionAtArrival,
        const double initialTime,
        const double timeOfFlight,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        const bool terminationSphereOfInfluence,
        std::map< double, Eigen::Vector6d >& patchedConicsResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult
        ){

    integratorSettings->initialTime_ = initialTime;

    // Compute the difference in state between the full problem and the Lambert targeter solution for the current leg.
    propagators::propagateLambertTargeterAndFullProblem(
        timeOfFlight, initialTime, bodyMap, accelerationMap, bodyToPropagate, centralBody, integratorSettings,
        patchedConicsResult, fullProblemResult,
        departureAndArrivalBodies, terminationSphereOfInfluence, cartesianPositionAtDeparture,
        cartesianPositionAtArrival, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN);

}

void propagateMga1DsmVelocityAndFullProblem(
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const std::vector< std::string > departureAndArrivalBodies,
        const std::string dsm,
        const std::string centralBody,
        const std::string bodyToPropagate,
        const Eigen::Vector3d cartesianPositionAtDeparture,
        const Eigen::Vector3d cartesianPositionDSM,
        const Eigen::Vector3d cartesianPositionAtArrival,
        const double initialTime,
        const double timeDsm,
        const double timeArrival,
        const transfer_trajectories::TransferLegType& legType,
        const std::vector< double >& trajectoryVariableVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        const bool terminationSphereOfInfluence,
        std::map< double, Eigen::Vector6d >& patchedConicsResultFromDepartureToDsm,
        std::map< double, Eigen::Vector6d >& fullProblemResultFromDepartureToDsm,
        std::map< double, Eigen::Vector6d >& patchedConicsResultFromDsmToArrival,
        std::map< double, Eigen::Vector6d >& fullProblemResultFromDsmToArrival,
        Eigen::Vector3d& velocityBeforeArrival,
        Eigen::Vector3d& velocityAfterDeparture,
        const double semiMajorAxis,
        const double eccentricity)
        {




//    simulation_setup::NamedBodyMap& bodyMap,
//    const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
//    const std::vector< std::string >& transferBodyOrder,
//    const std::string& centralBody,
//    const std::string& bodyToPropagate,
//    const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
//    const std::vector< double >& trajectoryVariableVector,
//    const std::vector< double >& minimumPericenterRadiiVector,
//    const std::vector< double >& semiMajorAxesVector,
//    const std::vector< double >& eccentricitiesVector,

//    )
//    {

//    std::vector< std::string > departureAndArrivalBodies;
//    departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs] );
//    departureAndArrivalBodies.push_back( bodiesAndManoeuvresOrder[counterLegs+1]);

    // Define cartesian position at departure and arrival for the current leg.
//    Eigen::Vector3d cartesianPositionAtDepartureCurrentLeg = positionVector[counterLegs];
//    Eigen::Vector3d cartesianPositionAtArrivalCurrentLeg = positionVector[counterLegs+1];



    if (legType == transfer_trajectories::mga1DsmVelocity_Departure){ // && legTypeVector[i] != transfer_trajectories::mga_Swingby){

        std::shared_ptr< transfer_trajectories::DepartureLegMga1DsmVelocity > departureLegMga1DsmVelocity =
                std::make_shared< transfer_trajectories::DepartureLegMga1DsmVelocity >(
                     cartesianPositionAtDeparture, cartesianPositionAtArrival, timeArrival - initialTime,
                     bodyMap[ departureAndArrivalBodies[0] ]->getEphemeris()->getCartesianState(initialTime).segment(3,3),
                     bodyMap[ centralBody ]->getGravityFieldModel()->getGravitationalParameter(),
                     bodyMap[ departureAndArrivalBodies[0] ]->getGravityFieldModel()->getGravitationalParameter(),
                     semiMajorAxis, eccentricity,
                     trajectoryVariableVector[ 0 ],
                     trajectoryVariableVector[ 1 ],
                     trajectoryVariableVector[ 2 ],
                     trajectoryVariableVector[ 3 ], true  );

        double deltaV;
        Eigen::Vector3d departureBodyPosition;
        Eigen::Vector3d departureBodyVelocity;

        departureLegMga1DsmVelocity->calculateLeg( velocityBeforeArrival, deltaV );

        departureLegMga1DsmVelocity->returnDepartureVariables( departureBodyPosition, departureBodyVelocity, velocityAfterDeparture );

    }

    if (legType == transfer_trajectories::mga1DsmVelocity_Swingby){


//        Eigen::Vector3d velocityBeforePlanet = velocityBeforeArrival;
        std::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforeArrival = std::make_shared< Eigen::Vector3d > ( velocityBeforeArrival );


        std::shared_ptr< transfer_trajectories::MissionLeg > missionLeg;
        std::shared_ptr< transfer_trajectories::SwingbyLegMga1DsmVelocity > swingbyLegMga1DsmVelocity =
                std::make_shared< transfer_trajectories::SwingbyLegMga1DsmVelocity >(
                     cartesianPositionAtDeparture, cartesianPositionAtArrival, timeArrival - initialTime,
                     bodyMap[ departureAndArrivalBodies[0] ]->getEphemeris()->getCartesianState(initialTime).segment(3,3),
                     bodyMap[ centralBody ]->getGravityFieldModel()->getGravitationalParameter(),
                     bodyMap[ departureAndArrivalBodies[0] ]->getGravityFieldModel()->getGravitationalParameter(),
                     pointerToVelocityBeforeArrival,
                     trajectoryVariableVector[ 0 ],
                     trajectoryVariableVector[ 1 ],
                     trajectoryVariableVector[ 2 ],
                     trajectoryVariableVector[ 3 ]);

        double deltaV;
        Eigen::Vector3d departureBodyPosition;
        Eigen::Vector3d departureBodyVelocity;

        swingbyLegMga1DsmVelocity->calculateLeg( velocityBeforeArrival, deltaV );

        swingbyLegMga1DsmVelocity->returnDepartureVariables( departureBodyPosition, departureBodyVelocity, velocityAfterDeparture );


    }

        // FIRST PART OF THE LEG: FROM DEPARTURE BODY TO DSM

//        std::map< double, Eigen::Vector6d > lambertTargeterResultCurrentLeg;
//        std::map< double, Eigen::Vector6d > fullProblemResultCurrentLeg;

        integratorSettings->initialTime_ = initialTime;

        std::vector< std::string >legDepartureAndArrival;
        legDepartureAndArrival.push_back( departureAndArrivalBodies[0] );
        legDepartureAndArrival.push_back( dsm );

        propagatePatchedConicsLegAndFullProblem(
            timeDsm - initialTime, initialTime, bodyMap, accelerationMap, bodyToPropagate, centralBody, integratorSettings,
            patchedConicsResultFromDepartureToDsm, fullProblemResultFromDepartureToDsm,
            legDepartureAndArrival, terminationSphereOfInfluence, cartesianPositionAtDeparture,
            cartesianPositionDSM, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, velocityAfterDeparture);

//        lambertTargeterResultForEachLeg[counterLegs] = lambertTargeterResultCurrentLeg;
//        fullProblemResultForEachLeg[counterLegs] = fullProblemResultCurrentLeg;

//        counterLegs++;

        // SECOND PART OF THE LEG: FROM DSM TO ARRIVAL BODY

        // Compute the difference in state between the full problem and the Lambert targeter solution for the current leg.
        legDepartureAndArrival.clear();
        legDepartureAndArrival.push_back( dsm );
        legDepartureAndArrival.push_back( departureAndArrivalBodies[1] );

//        cartesianPositionAtDepartureCurrentLeg = positionVector[counterLegs];
//        cartesianPositionAtArrivalCurrentLeg = positionVector[counterLegs+1];

//        lambertTargeterResultCurrentLeg.clear();
//        fullProblemResultCurrentLeg.clear();

        integratorSettings->initialTime_ = timeDsm;


        propagateLambertTargeterAndFullProblem(
            timeArrival - timeDsm, timeDsm, bodyMap, accelerationMap, bodyToPropagate, centralBody, integratorSettings,
            patchedConicsResultFromDsmToArrival, fullProblemResultFromDsmToArrival,
            legDepartureAndArrival, terminationSphereOfInfluence, cartesianPositionDSM,
            cartesianPositionAtArrival, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN);

}




void propagatePatchedConicsLegAndFullProblem(
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
        std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        const std::vector<std::string>& departureAndArrivalBodies,
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const Eigen::Vector3d& velocityAfterDeparture)
{
    // Clear output maps
    lambertTargeterResult.clear( );
    fullProblemResult.clear( );

    // Retrieve the gravitational parameter of the relevant bodies.
    double gravitationalParameterCentralBody = ( centralBodyGravitationalParameter == centralBodyGravitationalParameter ) ?
                centralBodyGravitationalParameter :
                bodyMap.at( centralBody )->getGravityFieldModel( )->getGravitationalParameter( );

    // Get halved value of the time of flight, later used as initial time for the propagation.
    double halvedTimeOfFlight = timeOfFlight / 2.0;

    // Time at the end of the transfer
    double finalTime = initialTime + timeOfFlight;

    // Retrieve positions of departure and arrival bodies from ephemerides
    Eigen::Vector3d cartesianPositionAtDepartureForLambertTargeter, cartesianPositionAtArrivalForLambertTargeter;
    if ( cartesianPositionAtDeparture != cartesianPositionAtDeparture )
    {
        // Cartesian state at departure
        if ( bodyMap.at( departureAndArrivalBodies.at( 0 ) )->getEphemeris( ) == nullptr)
        {
            throw std::runtime_error( "Ephemeris not defined for departure body." );
        }
        else
        {
            Eigen::Vector6d cartesianStateDepartureBody =
                    bodyMap.at( departureAndArrivalBodies.at( 0 ) )->getEphemeris( )->getCartesianState( initialTime);
            cartesianPositionAtDepartureForLambertTargeter = cartesianStateDepartureBody.segment(0,3);
        }
    }
    else
    {
        cartesianPositionAtDepartureForLambertTargeter = cartesianPositionAtDeparture;
    }

    if( cartesianPositionAtArrival != cartesianPositionAtArrival )
    {

        // Cartesian state at arrival
        if ( bodyMap.at( departureAndArrivalBodies.at( 1 ) )->getEphemeris( ) == nullptr){
            throw std::runtime_error( "Ephemeris not defined for arrival body." );
        }
        else{
            Eigen::Vector6d cartesianStateArrivalBody =
                    bodyMap.at( departureAndArrivalBodies.at( 1 ) )->getEphemeris( )->getCartesianState(finalTime);
            cartesianPositionAtArrivalForLambertTargeter =  cartesianStateArrivalBody.segment(0,3);
        }
    }
    else
    {
        cartesianPositionAtArrivalForLambertTargeter = cartesianPositionAtArrival;
    }

    Eigen::Vector6d cartesianStateAtDeparture;
    cartesianStateAtDeparture.segment(0,3) = cartesianPositionAtDepartureForLambertTargeter;
    cartesianStateAtDeparture.segment(3,3) = velocityAfterDeparture;


    // Compute cartesian state at halved time of flight.
    Eigen::Vector6d initialStatePropagationCartesianElements = computeCartesianStateHalfTimeOfFlightLambertTargeter(
                cartesianStateAtDeparture, gravitationalParameterCentralBody, timeOfFlight);



    // Initialise variables for propagatation.
    std::vector< std::string > centralBodiesPropagation;
    centralBodiesPropagation.push_back( "SSB" );
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back(bodyToPropagate);

    Eigen::Vector6d cartesianStatePatchedConics = initialStatePropagationCartesianElements;


    // Define forward propagator settings variables.
    integratorSettings->initialTime_ = initialTime + halvedTimeOfFlight;

    // Define forward propagation settings
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsForwardPropagation;
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsBackwardPropagation;

    propagatorSettingsForwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > > (
                centralBodiesPropagation, accelerationModelMap, bodiesToPropagate, initialStatePropagationCartesianElements,
                terminationSettings.second );
    propagatorSettingsBackwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > > (
                centralBodiesPropagation, accelerationModelMap, bodiesToPropagate, initialStatePropagationCartesianElements,
                terminationSettings.first );


    // Perform forward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards(
                bodyMap, integratorSettings, propagatorSettingsForwardPropagation );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.
            getEquationsOfMotionNumericalSolution( );

    // Calculate the difference between the full problem and the Lambert targeter solution along the forward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {

        cartesianStatePatchedConics = orbital_element_conversions::convertKeplerianToCartesianElements(
             orbital_element_conversions::propagateKeplerOrbit(
                        orbital_element_conversions::convertCartesianToKeplerianElements( initialStatePropagationCartesianElements,
                                                                                          gravitationalParameterCentralBody ),
                        itr->first - ( initialTime + halvedTimeOfFlight ),
                        gravitationalParameterCentralBody), gravitationalParameterCentralBody);

        lambertTargeterResult[ itr->first ] = cartesianStatePatchedConics;
        fullProblemResult[ itr->first ] = itr->second;
    }


    cartesianStatePatchedConics = initialStatePropagationCartesianElements;

    // Define backward propagator settings variables.
    integratorSettings->initialTimeStep_ = -integratorSettings->initialTimeStep_;
    integratorSettings->initialTime_ = initialTime + halvedTimeOfFlight;

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards(bodyMap, integratorSettings, propagatorSettingsBackwardPropagation );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation =
            dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );

    // Calculate the difference between the full problem and the Lambert targeter solution along the backward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
    {

        cartesianStatePatchedConics = orbital_element_conversions::convertKeplerianToCartesianElements(
             orbital_element_conversions::propagateKeplerOrbit(
                        orbital_element_conversions::convertCartesianToKeplerianElements( initialStatePropagationCartesianElements,
                                                                                          gravitationalParameterCentralBody ),
                        - (initialTime + halvedTimeOfFlight) + itr->first,
                        gravitationalParameterCentralBody), gravitationalParameterCentralBody);


        lambertTargeterResult[ itr->first ] = cartesianStatePatchedConics;
        fullProblemResult[ itr->first ] = itr->second;

    }

    integratorSettings->initialTimeStep_ = -integratorSettings->initialTimeStep_;

}


//! Function to propagate the full dynamics problem and the Lambert targeter solution.
void propagatePatchedConicsLegAndFullProblem(
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        const std::vector<std::string>& departureAndArrivalBodies,
        const bool terminationSphereOfInfluence,
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double departureBodyGravitationalParameter,
        const double arrivalBodyGravitationalParameter,
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& velocityAfterDeparture)
{

    // Retrieve the gravitational parameter of the relevant bodies.
    double gravitationalParameterCentralBody = ( centralBodyGravitationalParameter == centralBodyGravitationalParameter ) ?
                centralBodyGravitationalParameter :
                bodyMap.at( centralBody )->getGravityFieldModel( )->getGravitationalParameter( );


    // Retrieve positions of departure and arrival bodies from ephemerides if not defined as inputs.
    Eigen::Vector3d cartesianPositionAtDepartureForLambertTargeter, cartesianPositionAtArrivalForLambertTargeter;
    if ( cartesianPositionAtDeparture != cartesianPositionAtDeparture )
    {
        // Cartesian state at departure
        if ( bodyMap.at( departureAndArrivalBodies.at( 0 ) )->getEphemeris( ) == nullptr)
        {
            throw std::runtime_error( "Ephemeris not defined for departure body." );
        }
        else
        {
            Eigen::Vector6d cartesianStateDepartureBody =
                    bodyMap.at( departureAndArrivalBodies.at( 0 ) )->getEphemeris( )->getCartesianState( initialTime);
            cartesianPositionAtDepartureForLambertTargeter = cartesianStateDepartureBody.segment(0,3);
        }
    }
    else
    {
        cartesianPositionAtDepartureForLambertTargeter = cartesianPositionAtDeparture;
    }

    if( cartesianPositionAtArrival != cartesianPositionAtArrival )
    {

        // Cartesian state at arrival
        if ( bodyMap.at( departureAndArrivalBodies.at( 1 ) )->getEphemeris( ) == nullptr){
            throw std::runtime_error( "Ephemeris not defined for arrival body." );
        }
        else{
            Eigen::Vector6d cartesianStateArrivalBody =
                    bodyMap.at( departureAndArrivalBodies.at( 1 ) )->getEphemeris( )->getCartesianState( initialTime + timeOfFlight );
            cartesianPositionAtArrivalForLambertTargeter =  cartesianStateArrivalBody.segment(0,3);
        }
    }
    else
    {
        cartesianPositionAtArrivalForLambertTargeter = cartesianPositionAtArrival;
    }



    // Calculate radii sphere of influence about departure and arrival bodies
    double radiusSphereOfInfluenceDeparture;
    double radiusSphereOfInfluenceArrival;

    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettings;

    if (terminationSphereOfInfluence == true)
    {

        double gravitationalParameterDepartureBody = ( departureBodyGravitationalParameter == departureBodyGravitationalParameter ) ?
                    departureBodyGravitationalParameter :
                    bodyMap.at( departureAndArrivalBodies[0] )->getGravityFieldModel( )->getGravitationalParameter( );
        double gravitationalParameterArrivalBody = ( arrivalBodyGravitationalParameter == arrivalBodyGravitationalParameter ) ?
                    arrivalBodyGravitationalParameter :
                    bodyMap.at( departureAndArrivalBodies[1] )->getGravityFieldModel( )->getGravitationalParameter( );

        double distanceDepartureToCentralBodies =
                bodyMap.at( centralBody )->getEphemeris( )->getCartesianState(
                    initialTime ).segment( 0, 3 ).norm( ) - cartesianPositionAtDepartureForLambertTargeter.segment( 0, 3 ).norm( );
        double distanceArrivalToCentralBodies =
                bodyMap.at( centralBody )->getEphemeris( )->getCartesianState(
                    initialTime + timeOfFlight ).segment( 0, 3 ).norm( ) - cartesianPositionAtArrivalForLambertTargeter.segment( 0, 3 ).norm( );


        // Calculate radius sphere of influence for departure body.
        radiusSphereOfInfluenceDeparture = tudat::mission_geometry::computeSphereOfInfluence(
                    distanceDepartureToCentralBodies, gravitationalParameterDepartureBody, gravitationalParameterCentralBody);

        // Calculate radius sphere of influence for arrival body.
        radiusSphereOfInfluenceArrival = tudat::mission_geometry::computeSphereOfInfluence(
                    distanceArrivalToCentralBodies, gravitationalParameterArrivalBody, gravitationalParameterCentralBody);

        // Calculate the synodic period.
        double orbitalPeriodDepartureBody = basic_astrodynamics::computeKeplerOrbitalPeriod(
              orbital_element_conversions::convertCartesianToKeplerianElements( bodyMap.at( departureAndArrivalBodies[0] )->
              getEphemeris()->getCartesianState(initialTime), gravitationalParameterCentralBody)[ orbital_element_conversions::semiMajorAxisIndex ],
              gravitationalParameterCentralBody, gravitationalParameterDepartureBody);

        double orbitalPeriodArrivalBody = basic_astrodynamics::computeKeplerOrbitalPeriod(
              orbital_element_conversions::convertCartesianToKeplerianElements( bodyMap.at( departureAndArrivalBodies[1] )->
              getEphemeris()->getCartesianState(initialTime), gravitationalParameterCentralBody)[ orbital_element_conversions::semiMajorAxisIndex ],
              gravitationalParameterCentralBody, gravitationalParameterArrivalBody);

        double synodicPeriod;
        if (orbitalPeriodDepartureBody < orbitalPeriodArrivalBody){
            synodicPeriod = basic_astrodynamics::computeSynodicPeriod(orbitalPeriodDepartureBody, orbitalPeriodArrivalBody);
        }
        else {
            synodicPeriod = basic_astrodynamics::computeSynodicPeriod(orbitalPeriodArrivalBody, orbitalPeriodDepartureBody);
        }


        // Create total propagator termination settings.
        std::vector< std::shared_ptr< PropagationTerminationSettings > >  forwardPropagationTerminationSettingsList;
        forwardPropagationTerminationSettingsList.push_back(
                    std::make_shared< PropagationDependentVariableTerminationSettings >(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            relative_distance_dependent_variable, bodyToPropagate, departureAndArrivalBodies[ 1 ] ), radiusSphereOfInfluenceArrival, false ) );
        forwardPropagationTerminationSettingsList.push_back(
                    std::make_shared< PropagationTimeTerminationSettings >( 2 * synodicPeriod ) );


        std::shared_ptr< PropagationTerminationSettings > forwardPropagationTerminationSettings =
                std::make_shared< PropagationHybridTerminationSettings >( forwardPropagationTerminationSettingsList, true );


        std::vector< std::shared_ptr< PropagationTerminationSettings > >  backwardPropagationTerminationSettingsList;
        backwardPropagationTerminationSettingsList.push_back(
            std::make_shared< PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    relative_distance_dependent_variable, bodyToPropagate, departureAndArrivalBodies[ 0 ] ), radiusSphereOfInfluenceDeparture, false ) );
        backwardPropagationTerminationSettingsList.push_back(
                    std::make_shared< PropagationTimeTerminationSettings >( 2 * synodicPeriod ) );

        \
        std::shared_ptr< PropagationTerminationSettings > backwardPropagationTerminationSettings =
                std::make_shared< PropagationHybridTerminationSettings >( backwardPropagationTerminationSettingsList, true );

        terminationSettings = std::make_pair( backwardPropagationTerminationSettings, forwardPropagationTerminationSettings );


    }
    else
    {


        terminationSettings = std::make_pair(
                    std::make_shared< propagators::PropagationTimeTerminationSettings >( initialTime ),
                    std::make_shared< propagators::PropagationTimeTerminationSettings >( initialTime + timeOfFlight ) );
    }

    propagatePatchedConicsLegAndFullProblem(
            timeOfFlight, initialTime, bodyMap, accelerationModelMap, bodyToPropagate, centralBody, terminationSettings,
            integratorSettings, lambertTargeterResult, fullProblemResult, departureAndArrivalBodies,
            centralBodyGravitationalParameter, cartesianPositionAtDeparture, cartesianPositionAtArrival, velocityAfterDeparture );
}








//! Function to calculate the patched conics trajectory and to propagate the corresponding full problem
//! with the same acceleration map for every leg.
void fullPropagationPatchedConicsTrajectory(
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
        const bool terminationSphereOfInfluence,
        std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg)
{
    std::vector< basic_astrodynamics::AccelerationMap > accelerationMapForEachLeg;
    int numberOfLegs = legTypeVector.size(  );

    // Create vector with identical acceleration maps.
    for (int i = 0 ; i < numberOfLegs; i++){
        accelerationMapForEachLeg.push_back( accelerationMap );
    }

    // Compute difference between patched conics trajectory and full problem.
    fullPropagationPatchedConicsTrajectoryNewVersion(
                bodyMap, accelerationMapForEachLeg, transferBodyOrder, centralBody, bodyToPropagate, legTypeVector,
                trajectoryVariableVector, minimumPericenterRadiiVector, semiMajorAxesVector, eccentricitiesVector,
                integratorSettings, terminationSphereOfInfluence, lambertTargeterResultForEachLeg, fullProblemResultForEachLeg);

}




//! Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem,
//! at both departure and arrival positions for each leg.
std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullProblemWrtPatchedConicsTrajectory(
        simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< basic_astrodynamics::AccelerationMap >& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const std::vector<double>& trajectoryVariableVector,
        const std::vector<double>& minimumPericenterRadiiVector,
        const std::vector<double>& semiMajorAxesVector,
        const std::vector<double>& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        const bool terminationSphereOfInfluence )
{

    int numberOfLegs = legTypeVector.size( );
    int numberLegsIncludingDSM = ((trajectoryVariableVector.size() - 1 - numberOfLegs) / 4.0) + numberOfLegs ;


    // Compute difference between patched conics and full problem along the trajectory.
    std::map< int, std::map< double, Eigen::Vector6d > > lambertTargeterResultForEachLeg;
    std::map< int, std::map< double, Eigen::Vector6d > > fullProblemResultForEachLeg;

    fullPropagationPatchedConicsTrajectoryNewVersion(
                bodyMap, accelerationMap, transferBodyOrder, centralBody, bodyToPropagate, legTypeVector,
                trajectoryVariableVector, minimumPericenterRadiiVector, semiMajorAxesVector, eccentricitiesVector,
                integratorSettings, terminationSphereOfInfluence, lambertTargeterResultForEachLeg, fullProblemResultForEachLeg );


    // Compute difference at departure and at arrival for each leg (considering that a leg including a DSM consists of two sub-legs).
    std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > stateDifferenceAtArrivalAndDepartureForEachLeg;

    for (int i = 0 ; i < numberLegsIncludingDSM - 1 ; i++){

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




//! Function to compute the difference in cartesian state between patched conics trajectory and full dynamics problem,
//! at both departure and arrival positions for each leg when using the same accelerations for each leg.
std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullProblemWrtPatchedConicsTrajectory(
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const std::vector< std::string >& transferBodyOrder,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const std::vector<double>& trajectoryVariableVector,
        const std::vector<double>& minimumPericenterRadiiVector,
        const std::vector<double>& semiMajorAxesVector,
        const std::vector<double>& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        const bool terminationSphereOfInfluence )
{

    int numberOfLegs = legTypeVector.size(  );

    // Create vector with identical acceleration maps.
    std::vector< basic_astrodynamics::AccelerationMap > accelerationMapForEachLeg;
    for (int i = 0 ; i < numberOfLegs; i++){
        accelerationMapForEachLeg.push_back( accelerationMap );
    }


    // Compute difference at departure and at arrival for each leg (considering that a leg including a DSM consists of two sub-legs).
    std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > stateDifferenceAtArrivalAndDepartureForEachLeg;

    stateDifferenceAtArrivalAndDepartureForEachLeg = getDifferenceFullProblemWrtPatchedConicsTrajectory(
                bodyMap, accelerationMapForEachLeg,
                transferBodyOrder, centralBody, bodyToPropagate,
                legTypeVector, trajectoryVariableVector,
                minimumPericenterRadiiVector, semiMajorAxesVector,
                eccentricitiesVector, integratorSettings, terminationSphereOfInfluence );
    return stateDifferenceAtArrivalAndDepartureForEachLeg;

}



}

}


