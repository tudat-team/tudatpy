/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/createStateDerivativeModel.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterIzzo.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>

#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"

#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/planetTrajectory.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"

namespace tudat
{

namespace propagators
{

//! Function to setup a body map corresponding to the assumptions of the Lambert targeter,
//! using default ephemerides for the central, departure and arrival bodies.
simulation_setup::NamedBodyMap setupBodyMapFromEphemeridesForLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& departureAndArrivalBodies )
{

    spice_interface::loadStandardSpiceKernels( );

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( nameCentralBody );
    bodiesToCreate.push_back( departureAndArrivalBodies[0] );
    bodiesToCreate.push_back( departureAndArrivalBodies[1] );


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

    bodyMap[ nameBodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( nameBodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );


    return bodyMap;
}



//! Function to setup a body map corresponding to the assumptions of the patched conics trajectory,
//! the ephemerides of the transfer bodies being provided as inputs.
simulation_setup::NamedBodyMap setupBodyMapFromUserDefinedEphemeridesForLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& departureAndArrivalBodies,
        const std::vector< ephemerides::EphemerisPointer >& ephemerisVectorDepartureAndArrivalBodies)
{

    spice_interface::loadStandardSpiceKernels( );

    // Create central body object.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( nameCentralBody );

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

    bodyMap[ nameBodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( nameBodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );



    // Define ephemeris for the departure and arrival bodies.
    for ( unsigned int i = 0 ; i < departureAndArrivalBodies.size() ; i++){

        bodyMap[ departureAndArrivalBodies[i] ] = std::make_shared< simulation_setup::Body >( );
        bodyMap[ departureAndArrivalBodies[i] ]->setEphemeris( ephemerisVectorDepartureAndArrivalBodies[i] );

    }


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    return bodyMap;

}




//! Function to setup a body map corresponding to the assumptions of the Lambert targeter,
//! using default ephemerides for the central body only, while the positions of departure and arrival bodies are provided as inputs.
simulation_setup::NamedBodyMap setupBodyMapFromUserDefinedStatesForLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& departureAndArrivalBodies,
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival )
{

    spice_interface::loadStandardSpiceKernels( );


    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Create ephemeris vector for departure and arrival bodies.
    std::vector< ephemerides::EphemerisPointer > ephemerisVectorDepartureAndArrivalBodies(2);

    ephemerisVectorDepartureAndArrivalBodies[ 0 ] = std::make_shared< ephemerides::ConstantEphemeris > (
                ( Eigen::Vector6d( ) << cartesianPositionAtDeparture[0], cartesianPositionAtDeparture[1], cartesianPositionAtDeparture[2],
                   0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    ephemerisVectorDepartureAndArrivalBodies[ 1 ] = std::make_shared< ephemerides::ConstantEphemeris >(
                ( Eigen::Vector6d( ) << cartesianPositionAtArrival[0], cartesianPositionAtArrival[1], cartesianPositionAtArrival[2],
                  0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );



    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = setupBodyMapFromUserDefinedEphemeridesForLambertTargeter(nameCentralBody, nameBodyToPropagate,
                                                                                                      departureAndArrivalBodies,
                                                                                                      ephemerisVectorDepartureAndArrivalBodies);

    return bodyMap;

}



//! Function to directly setup an acceleration map for the Lambert targeter.
basic_astrodynamics::AccelerationMap setupAccelerationMapLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const simulation_setup::NamedBodyMap& bodyMap )
{

    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( nameBodyToPropagate );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( nameCentralBody );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[nameCentralBody].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ nameBodyToPropagate ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    return accelerationModelMap;

}





//! Function to determine the cartesian state at a given time for a keplerian orbit, based on the initial state.
Eigen::Vector6d computeCartesianStateFromKeplerianOrbit(
        const Eigen::Vector6d& initialState,
        const double finalPropagationTime,
        const double gravitationalParameter)
{
    Eigen::Vector6d keplerianInitialState = orbital_element_conversions::convertCartesianToKeplerianElements(
                initialState, gravitationalParameter);
    // Retrieve the semi-major axis and eccentricty of the keplerian orbit.
    double semiMajorAxis = keplerianInitialState[orbital_element_conversions::semiMajorAxisIndex];
    double eccentricity = keplerianInitialState[orbital_element_conversions::eccentricityIndex];

    // Calculate the initial mean anomaly.
    double initialTrueAnomaly = keplerianInitialState[orbital_element_conversions::trueAnomalyIndex];
    double initialMeanAnomaly;

    if( eccentricity < 1.0 )
    {
        initialMeanAnomaly = orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
                    orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly(initialTrueAnomaly, eccentricity), eccentricity);

    }
    else
    {
        initialMeanAnomaly = orbital_element_conversions::convertHyperbolicEccentricAnomalyToMeanAnomaly(
                    orbital_element_conversions::convertTrueAnomalyToHyperbolicEccentricAnomaly(
                        initialTrueAnomaly, eccentricity), eccentricity);
    }

    // Calculate the mean anomaly at the final time.
    double meanAnomalyEndPropagation = initialMeanAnomaly + orbital_element_conversions::convertElapsedTimeToMeanAnomalyChange(
                finalPropagationTime, gravitationalParameter, semiMajorAxis );

    double trueAnomalyEndPropagation;

    if( eccentricity < 1.0 )
    {
        trueAnomalyEndPropagation = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                    orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly(
                        eccentricity, meanAnomalyEndPropagation), eccentricity );
    }
    else
    {
        trueAnomalyEndPropagation = orbital_element_conversions::convertHyperbolicEccentricAnomalyToTrueAnomaly(
                    orbital_element_conversions::convertMeanAnomalyToHyperbolicEccentricAnomaly(
                        eccentricity, meanAnomalyEndPropagation), eccentricity );
    }

    // Determine the keplerian state at final time.
    Eigen::Vector6d finalKeplerianState = keplerianInitialState;
    finalKeplerianState[orbital_element_conversions::trueAnomalyIndex] = trueAnomalyEndPropagation;

    // Convert keplerian to cartesian state at final time.
    Eigen::Vector6d cartesianStateLambertSolution = orbital_element_conversions::convertKeplerianToCartesianElements(
                finalKeplerianState, gravitationalParameter);

    return cartesianStateLambertSolution;
}



//! Function to propagate the full dynamics problem and the Lambert targeter solution.
void propagateLambertTargeterAndFullProblem(
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string& centralBody,
        std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        std::map< double, Eigen::VectorXd >& dependentVariableResult,
        const std::vector<std::string>& departureAndArrivalBodies,
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival )
{
    // Clear output maps
    lambertTargeterResult.clear( );
    fullProblemResult.clear( );
    dependentVariableResult.clear( );

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


    // Run the Lambert targeter.
    mission_segments::LambertTargeterIzzo lambertTargeter(
                cartesianPositionAtDepartureForLambertTargeter, cartesianPositionAtArrivalForLambertTargeter,
                timeOfFlight, gravitationalParameterCentralBody );

    // Retrieve cartesian state at departure.
    Eigen::Vector3d cartesianVelocityAtDeparture = lambertTargeter.getInertialVelocityAtDeparture( );
    Eigen::Vector6d cartesianStateAtDeparture;
    cartesianStateAtDeparture.segment(0,3) = cartesianPositionAtDepartureForLambertTargeter;
    cartesianStateAtDeparture.segment(3,3) = cartesianVelocityAtDeparture;


    // Convert into keplerian elements
    Eigen::Vector6d keplerianStateAtDeparture = orbital_element_conversions::convertCartesianToKeplerianElements(
                cartesianStateAtDeparture, gravitationalParameterCentralBody );

    // Propagate the keplerian elements until half of the time of flight.
    Eigen::Vector6d keplerianStateAtHalvedTimeOfFlight = orbital_element_conversions::propagateKeplerOrbit( keplerianStateAtDeparture,
                halvedTimeOfFlight, gravitationalParameterCentralBody );

    // Convert the keplerian elements back into Cartesian elements.
    Eigen::Vector6d initialStatePropagationCartesianElements = orbital_element_conversions::convertKeplerianToCartesianElements(
                keplerianStateAtHalvedTimeOfFlight, gravitationalParameterCentralBody );


    Eigen::Vector6d cartesianStateLambertSolution;


    // Define forward propagator settings variables.
    integratorSettings->initialTime_ = initialTime + halvedTimeOfFlight;

    // Define forward propagation settings
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsForwardPropagation;
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsBackwardPropagation;


    propagatorSettingsForwardPropagation = propagatorSettings.second;
    propagatorSettingsForwardPropagation->resetInitialStates( initialStatePropagationCartesianElements );

    propagatorSettingsBackwardPropagation = propagatorSettings.first;
    propagatorSettingsBackwardPropagation->resetInitialStates( initialStatePropagationCartesianElements );

    // Perform forward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards(
                bodyMap, integratorSettings, propagatorSettingsForwardPropagation );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.
            getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation =
            dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );

    // Calculate the difference between the full problem and the Lambert targeter solution along the forward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {
        cartesianStateLambertSolution = computeCartesianStateFromKeplerianOrbit(
                    initialStatePropagationCartesianElements, itr->first - ( initialTime + halvedTimeOfFlight ),
                    gravitationalParameterCentralBody );
        lambertTargeterResult[ itr->first ] = cartesianStateLambertSolution;
        fullProblemResult[ itr->first ] = itr->second;
        dependentVariableResult[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
    }


    // Define backward propagator settings variables.
    integratorSettings->initialTimeStep_ = -integratorSettings->initialTimeStep_;
    integratorSettings->initialTime_ = initialTime + halvedTimeOfFlight;

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards(bodyMap, integratorSettings, propagatorSettingsBackwardPropagation );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation =
            dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryBackwardPropagation =
            dynamicsSimulatorIntegrationBackwards.getDependentVariableHistory( );

    // Calculate the difference between the full problem and the Lambert targeter solution along the backward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
    {
        cartesianStateLambertSolution = computeCartesianStateFromKeplerianOrbit(
                    initialStatePropagationCartesianElements, - (initialTime + halvedTimeOfFlight) + itr->first,
                    gravitationalParameterCentralBody);

        lambertTargeterResult[ itr->first ] = cartesianStateLambertSolution;
        fullProblemResult[ itr->first ] = itr->second;
        dependentVariableResult[ itr->first ] = dependentVariableHistoryBackwardPropagation[ itr->first ];

    }

//    std::cout<<"In loop: "<<std::endl<<std::setprecision( 16 )<<
//               lambertTargeterResult.begin( )->second.segment( 0, 3 ).transpose( )<<std::endl<<
//               cartesianPositionAtDeparture.transpose( )<<std::endl<<
//               lambertTargeterResult.rbegin( )->second.segment( 0, 3 ).transpose( )<<std::endl<<
//               cartesianPositionAtArrival.transpose( )<<std::endl<<std::endl;



    integratorSettings->initialTimeStep_ = -integratorSettings->initialTimeStep_;

}


//! Function to propagate the full dynamics problem and the Lambert targeter solution.
void propagateLambertTargeterAndFullProblem(
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        std::map< double, Eigen::VectorXd >& dependentVariableResult,
        const std::vector<std::string>& departureAndArrivalBodies,
        const bool terminationSphereOfInfluence,
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double departureBodyGravitationalParameter,
        const double arrivalBodyGravitationalParameter,
        const double centralBodyGravitationalParameter,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = std::shared_ptr< DependentVariableSaveSettings >( ),
        const TranslationalPropagatorType propagator = cowell)
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

    Eigen::Vector6d initialState;

    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings;

    std::vector< std::string > centralBodyPropagation; centralBodyPropagation.push_back( centralBody );
    std::vector< std::string > bodyToPropagatePropagation; bodyToPropagatePropagation.push_back( bodyToPropagate );

    propagatorSettings.first = std::make_shared< TranslationalStatePropagatorSettings< double > >(
        centralBodyPropagation, accelerationModelMap, bodyToPropagatePropagation, initialState,
        terminationSettings.first, propagator, dependentVariablesToSave );

    propagatorSettings.second = std::make_shared< TranslationalStatePropagatorSettings< double > >(
        centralBodyPropagation, accelerationModelMap, bodyToPropagatePropagation, initialState,
        terminationSettings.second, propagator, dependentVariablesToSave );

    propagateLambertTargeterAndFullProblem(
            timeOfFlight, initialTime, bodyMap, centralBody, propagatorSettings,
            integratorSettings, lambertTargeterResult, fullProblemResult, dependentVariableResult, departureAndArrivalBodies,
            centralBodyGravitationalParameter, cartesianPositionAtDeparture, cartesianPositionAtArrival );
}


//! Function to compute the difference in cartesian state between Lambert targeter solution and full dynamics problem, both at departure
//! and at arrival.
std::pair< Eigen::Vector6d, Eigen::Vector6d > getDifferenceFullPropagationWrtLambertTargeterAtDepartureAndArrival(
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double timeOfFlight,
        const double initialTime,
        const simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::vector< std::string >& departureAndArrivalBodies,
        const bool terminationSphereOfInfluence,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = std::shared_ptr< DependentVariableSaveSettings >( ),
        const TranslationalPropagatorType propagator = cowell)

{
    std::map< double, Eigen::Vector6d > lambertTargeterResult;
    std::map< double, Eigen::Vector6d > fullProblemResult;
    std::map< double, Eigen::VectorXd > dependentVariableResult;



    // Compute full problem and Lambert targeter solution at both departure and arrival.
    propagateLambertTargeterAndFullProblem(
                timeOfFlight, initialTime, bodyMap, accelerationModelMap, bodyToPropagate, centralBody, integratorSettings,
                lambertTargeterResult, fullProblemResult, dependentVariableResult, departureAndArrivalBodies,
                terminationSphereOfInfluence, cartesianPositionAtDeparture, cartesianPositionAtArrival, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN,
                dependentVariablesToSave, propagator );

    Eigen::Vector6d stateLambertTargeterAtDeparture = lambertTargeterResult.begin( )->second;
    Eigen::Vector6d propagatedStateFullProblemAtDeparture = fullProblemResult.begin( )->second;
    Eigen::Vector6d stateLambertTargeterAtArrival = lambertTargeterResult.rbegin( )->second;
    Eigen::Vector6d propagatedStateFullProblemAtArrival = fullProblemResult.rbegin( )->second;


    // Difference between the Lambert targeter and full problem results at departure and arrival.
    return std::make_pair( stateLambertTargeterAtDeparture - propagatedStateFullProblemAtDeparture,
                           stateLambertTargeterAtArrival - propagatedStateFullProblemAtArrival);
}




}

}
