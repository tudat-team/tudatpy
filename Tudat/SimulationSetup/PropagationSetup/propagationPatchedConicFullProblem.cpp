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

std::vector< double > getDefaultMinimumPericenterRadii( const std::vector< std::string >& bodyNames )
{
   std::vector< double > pericenterRadii;
   for( unsigned int i = 0; i < bodyNames.size( ); i++ )
   {
       if( bodyNames.at( i ) == "Mercury" )
       {
            pericenterRadii.push_back( 2639.7E3 );
       }
       else if( bodyNames.at( i ) == "Venus" )
       {
           pericenterRadii.push_back( 6251.8E3 );
       }
       else if( bodyNames.at( i ) == "Earth" )
       {
           pericenterRadii.push_back( 6578.1E3 );
       }
       else if( bodyNames.at( i ) == "Mars" )
       {
           pericenterRadii.push_back( 3596.2E3 );
       }
       else if( bodyNames.at( i ) == "Jupiter" )
       {
           pericenterRadii.push_back( 72000.0E3 );
       }
       else if( bodyNames.at( i ) == "Saturn" )
       {
           pericenterRadii.push_back( 61000.0E3 );
       }
       else if( bodyNames.at( i ) == "Uranus" )
       {
           pericenterRadii.push_back( 26000.0E3 );
       }
       else if( bodyNames.at( i ) == "Neptune" )
       {
           pericenterRadii.push_back( 25000.0E3 );
       }
       else if( bodyNames.at( i ) == "Pluto" )
       {
           pericenterRadii.push_back( 1395.0E3 );
       }
       else
       {
           throw std::runtime_error( "Error, could not recognize body " + bodyNames.at( i ) + " when getting minimum periapsis radius" );
       }
   }
   return pericenterRadii;
}


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
    for (int i = 0 ; i < numberOfLegs ; i++)
    {
        bodiesAndManoeuvresOrder.push_back(transferBodyOrder[i]);
        accelerationMapForEachLeg.push_back( accelerationMap[i] );

        if (legTypeVector[i] != transfer_trajectories::mga_Departure && legTypeVector[i] != transfer_trajectories::mga_Swingby)
        {
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

    for (int i = 0 ; i < numberOfLegs - 1 ; i ++)
    {
        if (legTypeVector[i] == transfer_trajectories::mga_Departure || legTypeVector[i] == transfer_trajectories::mga_Swingby )
        {
            timeOfFlight = trajectoryVariableVector[1 + counterLegTotal];
            timeOfFlightVector.push_back( timeOfFlight );

            totalTime += timeOfFlight;
            absoluteTimeVector.push_back(totalTime);

            counterLegTotal++;
        }
        // If a DSM is performed in the current leg
        else
        {
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
                    cartesianPositionAtDepartureCurrentLeg, cartesianPositionAtArrivalCurrentLeg,
                    timeOfFlightVector[i], absoluteTimeVector[i], bodyMap, accelerationMapForEachLeg[i], bodyToPropagate, centralBody,
                    integratorSettings, lambertTargeterResultCurrentLeg, fullProblemResultCurrentLeg, departureAndArrivalBodies, false,
                    terminationSphereOfInfluence );

        lambertTargeterResultForEachLeg[i] = lambertTargeterResultCurrentLeg;
        fullProblemResultForEachLeg[i] = fullProblemResultCurrentLeg;

    }
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
    fullPropagationPatchedConicsTrajectory(
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

    fullPropagationPatchedConicsTrajectory(
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


