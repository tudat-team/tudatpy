///*    Copyright (c) 2010-2019, Delft University of Technology
// *    All rigths reserved
// *
// *    This file is part of the Tudat. Redistribution and use in source and
// *    binary forms, with or without modification, are permitted exclusively
// *    under the terms of the Modified BSD license. You should have received
// *    a copy of the license with this file. If not, please or visit:
// *    http://tudat.tudelft.nl/LICENSE.
// */

//#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"
//#include "tudat/astro/mission_segments/lambertTargeter.h"
//#include "tudat/astro/mission_segments/lambertTargeterIzzo.h"
//#include "tudat/astro/mission_segments/lambertRoutines.h"
//#include "tudat/simulation/environment_setup/defaultBodies.h"
//#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
//#include "tudat/simulation/propagation_setup/propagationLambertTargeterFullProblem.h"
//#include <tudat/simulation/estimation.h>

//#include "tudat/astro/gravitation/librationPoint.h"
//#include "tudat/astro/basic_astro/celestialBodyConstants.h"
//#include "tudat/astro/basic_astro/missionGeometry.h"

//#include "tudat/astro/trajectory_design/trajectory.h"
//#include "tudat/astro/trajectory_design/exportTrajectory.h"
//#include "tudat/astro/trajectory_design/planetTrajectory.h"

//#include "tudat/interface/spice/spiceInterface.h"

//#include "tudat/io/basicInputOutput.h"

//namespace tudat
//{

//namespace propagators
//{

////! Function to setup a system of bodies corresponding to the assumptions of the Lambert targeter,
////! using default ephemerides for the central, departure and arrival bodies.
//simulation_setup::SystemOfBodies setupBodyMapFromEphemeridesForLambertTargeter(
//        const std::string& nameCentralBody,
//        const std::string& nameBodyToPropagate,
//        const std::pair< std::string, std::string >& departureAndArrivalBodies )
//{
    
//    spice_interface::loadStandardSpiceKernels( );
    
//    // Create central, departure and arrival bodies.
//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( nameCentralBody );
//    bodiesToCreate.push_back( departureAndArrivalBodies.first );
//    bodiesToCreate.push_back( departureAndArrivalBodies.second );
    
//    std::string frameOrigin = "SSB";
//    std::string frameOrientation = "ECLIPJ2000";
//    simulation_setup::BodyListSettings bodySettings =
//            simulation_setup::getDefaultBodySettings( bodiesToCreate, frameOrigin, frameOrientation );
    
//    // Define central body ephemeris settings.
//    bodySettings.at( nameCentralBody )->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
//                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );
    
//    bodySettings.at( nameCentralBody )->ephemerisSettings->resetFrameOrientation( frameOrientation );
//    bodySettings.at( nameCentralBody )->rotationModelSettings->resetOriginalFrame( frameOrientation );
//<<<<<<< HEAD


//    // Create system of bodies.
//    simulation_setup::SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//    bodies.createEmptyBody( nameBodyToPropagate );
//    bodies.at( nameBodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
//                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


//    return bodies;
//=======
    
    
//    // Create body map.
//    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );
    
//    bodyMap.addNewBody( nameBodyToPropagate );
//    bodyMap.at( nameBodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
//                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );
    
    
//    return bodyMap;
//>>>>>>> dominic-origin/features/mission_segments_refactor
//}



////! Function to setup a system of bodies corresponding to the assumptions of the patched conics trajectory,
////! the ephemerides of the transfer bodies being provided as inputs.
//simulation_setup::SystemOfBodies setupBodyMapFromUserDefinedEphemeridesForLambertTargeter(
//        const std::string& nameCentralBody,
//        const std::string& nameBodyToPropagate,
//        const std::pair< std::string, std::string >& departureAndArrivalBodies,
//        const std::vector< ephemerides::EphemerisPointer >& ephemerisVectorDepartureAndArrivalBodies)
//{
//    spice_interface::loadStandardSpiceKernels( );
    
//    // Create central body object.
//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( nameCentralBody );
    
//    std::string frameOrigin = "SSB";
//    std::string frameOrientation = "ECLIPJ2000";
//    simulation_setup::BodyListSettings bodySettings =
//            simulation_setup::getDefaultBodySettings( bodiesToCreate, frameOrigin, frameOrientation );
    
//    // Define central body ephemeris settings.
//    bodySettings.at( nameCentralBody )->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
//                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );
//<<<<<<< HEAD

//    bodySettings.at( nameCentralBody )->ephemerisSettings->resetFrameOrientation( frameOrientation );
//    bodySettings.at( nameCentralBody )->rotationModelSettings->resetOriginalFrame( frameOrientation );



//    // Create system of bodies.
//    simulation_setup::SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//    bodies.createEmptyBody( nameBodyToPropagate );
//    bodies.at( nameBodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
//=======
    
//    // Create body map.
//    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );
    
//    bodyMap.addNewBody( nameBodyToPropagate );
//    bodyMap.at( nameBodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
//>>>>>>> dominic-origin/features/mission_segments_refactor
//                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );

//    bodyMap.addNewBody( departureAndArrivalBodies.first );
//    bodyMap.at( departureAndArrivalBodies.first )->setEphemeris( ephemerisVectorDepartureAndArrivalBodies.at( 0 ) );

//<<<<<<< HEAD

//    // Define ephemeris for the departure and arrival bodies.
//    for ( unsigned int i = 0 ; i < departureAndArrivalBodies.size() ; i++){

//        bodies.createEmptyBody( departureAndArrivalBodies[i] );
//        bodies.at( departureAndArrivalBodies[i] )->setEphemeris( ephemerisVectorDepartureAndArrivalBodies[i] );

//    }


//    return bodies;

//}




////! Function to setup a system of bodies corresponding to the assumptions of the Lambert targeter,
////! using default ephemerides for the central body only, while the positions of departure and arrival bodies are provided as inputs.
//simulation_setup::SystemOfBodies setupBodyMapFromUserDefinedStatesForLambertTargeter(
//        const std::string& nameCentralBody,
//        const std::string& nameBodyToPropagate,
//        const std::vector< std::string >& departureAndArrivalBodies,
//        const Eigen::Vector3d& cartesianPositionAtDeparture,
//        const Eigen::Vector3d& cartesianPositionAtArrival )
//{

//    spice_interface::loadStandardSpiceKernels( );


//    std::string frameOrigin = "SSB";
//    std::string frameOrientation = "ECLIPJ2000";


//    // Create ephemeris vector for departure and arrival bodies.
//    std::vector< ephemerides::EphemerisPointer > ephemerisVectorDepartureAndArrivalBodies(2);

//    ephemerisVectorDepartureAndArrivalBodies[ 0 ] = std::make_shared< ephemerides::ConstantEphemeris > (
//                ( Eigen::Vector6d( ) << cartesianPositionAtDeparture[0], cartesianPositionAtDeparture[1], cartesianPositionAtDeparture[2],
//                   0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

//    ephemerisVectorDepartureAndArrivalBodies[ 1 ] = std::make_shared< ephemerides::ConstantEphemeris >(
//                ( Eigen::Vector6d( ) << cartesianPositionAtArrival[0], cartesianPositionAtArrival[1], cartesianPositionAtArrival[2],
//                  0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );



//    // Create system of bodies.
//    simulation_setup::SystemOfBodies bodies = setupBodyMapFromUserDefinedEphemeridesForLambertTargeter(nameCentralBody, nameBodyToPropagate,
//                                                                                                      departureAndArrivalBodies,
//                                                                                                      ephemerisVectorDepartureAndArrivalBodies);

//    return bodies;

//=======
//    bodyMap.addNewBody( departureAndArrivalBodies.second );
//    bodyMap.at( departureAndArrivalBodies.second )->setEphemeris( ephemerisVectorDepartureAndArrivalBodies.at( 1 ) );
    
    
//    return bodyMap;
    
//>>>>>>> dominic-origin/features/mission_segments_refactor
//}


////! Function to directly setup an acceleration map for the Lambert targeter.
//basic_astrodynamics::AccelerationMap setupAccelerationMapLambertTargeter(
//        const std::string& nameCentralBody,
//        const std::string& nameBodyToPropagate,
//        const simulation_setup::SystemOfBodies& bodies )
//{
    
//    std::vector< std::string > bodiesToPropagate;
//    bodiesToPropagate.push_back( nameBodyToPropagate );
//    std::vector< std::string > centralBodies;
//    centralBodies.push_back( nameCentralBody );
    
//    // Acceleration from the central body.
//    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
//    bodyToPropagateAccelerations[nameCentralBody].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
//                                                                basic_astrodynamics::central_gravity ) );
    
//    simulation_setup::SelectedAccelerationMap accelerationMap;
//    accelerationMap[ nameBodyToPropagate ] = bodyToPropagateAccelerations;
    
//    // Create the acceleration map.
//    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
//<<<<<<< HEAD
//                bodies, accelerationMap, bodiesToPropagate, centralBodies );


//=======
//                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );
    
    
//>>>>>>> dominic-origin/features/mission_segments_refactor
//    return accelerationModelMap;
    
//}

//std::shared_ptr< mission_segments::LambertTargeterIzzo > createLambertTargeterIzzo(
//        const double timeOfFlight,
//        const double initialTime,
//        const simulation_setup::NamedBodyMap& bodyMap,
//        const std::string& centralBody,
//        const std::pair< std::string, std::string >& departureAndArrivalBodies )
//{

//    // Cartesian state at departure
//    Eigen::Vector3d cartesianPositionAtDepartureForLambertTargeter;
//    if ( bodyMap.at( departureAndArrivalBodies.first )->getEphemeris( ) == nullptr)
//    {
//        throw std::runtime_error( "Ephemeris not defined for departure body." );
//    }
//    else
//    {
//        cartesianPositionAtDepartureForLambertTargeter =
//                bodyMap.at( departureAndArrivalBodies.first )->getEphemeris( )->getCartesianState( initialTime).segment( 0, 3 );
//    }

//    // Cartesian state at arrival
//    Eigen::Vector3d cartesianPositionAtArrivalForLambertTargeter;
//    if ( bodyMap.at( departureAndArrivalBodies.second )->getEphemeris( ) == nullptr)
//    {
//        throw std::runtime_error( "Ephemeris not defined for arrival body." );
//    }
//    else
//    {
//        cartesianPositionAtArrivalForLambertTargeter =
//                bodyMap.at( departureAndArrivalBodies.second )->getEphemeris( )->getCartesianState(
//                    initialTime + timeOfFlight ).segment( 0, 3 );
//    }

//    // Run the Lambert targeter.
//    double centralBodyGravitationalParameter =
//            bodyMap.at( centralBody )->getGravityFieldModel( )->getGravitationalParameter( );
//    return std::make_shared< mission_segments::LambertTargeterIzzo >(
//                cartesianPositionAtDepartureForLambertTargeter, cartesianPositionAtArrivalForLambertTargeter,
//                timeOfFlight, centralBodyGravitationalParameter );
//}

//void propagateLambertTargeterAndFullProblem(
//        const double timeOfFlight,
//        const double initialTime,
//        const simulation_setup::NamedBodyMap& bodyMap,
//        const std::string& centralBody,
//        const std::pair< std::string, std::string >& departureAndArrivalBodies ,
//        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
//        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
//        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
//        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
//        std::map< double, Eigen::Vector6d >& fullProblemResult )
//{
//    std::map< double, Eigen::VectorXd > dependentVariableResult;
//    propagateLambertTargeterAndFullProblem(
//                timeOfFlight, initialTime, bodyMap, centralBody, departureAndArrivalBodies, propagatorSettings,
//                integratorSettings, lambertTargeterResult,fullProblemResult );
//}

////! Function to propagate the full dynamics problem and the Lambert targeter solution.
//void propagateLambertTargeterAndFullProblem(
//        const double timeOfFlight,
//        const double initialTime,
//        const simulation_setup::SystemOfBodies& bodies,
//        const std::string& centralBody,
//        const std::pair< std::string, std::string >& departureAndArrivalBodies,
//        const std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
//        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
//        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
//        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
//        std::map< double, Eigen::Vector6d >& fullProblemResult,
//        std::map< double, Eigen::VectorXd >& dependentVariableResult )
//{

//    // Clear output maps
//    lambertTargeterResult.clear( );
//    fullProblemResult.clear( );
//    dependentVariableResult.clear( );

//<<<<<<< HEAD
//    // Retrieve the gravitational parameter of the relevant bodies.
//    double gravitationalParameterCentralBody = ( centralBodyGravitationalParameter == centralBodyGravitationalParameter ) ?
//                centralBodyGravitationalParameter :
//                bodies.at( centralBody )->getGravityFieldModel( )->getGravitationalParameter( );

//    // Get halved value of the time of flight, later used as initial time for the propagation.
//    double halvedTimeOfFlight = timeOfFlight / 2.0;

//    // Time at the end of the transfer
//    double finalTime = initialTime + timeOfFlight;

//    // Retrieve positions of departure and arrival bodies from ephemerides
//    Eigen::Vector3d cartesianPositionAtDepartureForLambertTargeter, cartesianPositionAtArrivalForLambertTargeter;
//    if ( cartesianPositionAtDeparture != cartesianPositionAtDeparture )
//    {
//        // Cartesian state at departure
//        if ( bodies.at( departureAndArrivalBodies.at( 0 ) )->getEphemeris( ) == nullptr)
//        {
//            throw std::runtime_error( "Ephemeris not defined for departure body." );
//        }
//        else
//        {
//            Eigen::Vector6d cartesianStateDepartureBody =
//                    bodies.at( departureAndArrivalBodies.at( 0 ) )->getEphemeris( )->getCartesianState( initialTime);
//            cartesianPositionAtDepartureForLambertTargeter = cartesianStateDepartureBody.segment(0,3);
//        }
//    }
//    else
//    {
//        cartesianPositionAtDepartureForLambertTargeter = cartesianPositionAtDeparture;
//    }

//    if( cartesianPositionAtArrival != cartesianPositionAtArrival )
//    {

//        // Cartesian state at arrival
//        if ( bodies.at( departureAndArrivalBodies.at( 1 ) )->getEphemeris( ) == nullptr){
//            throw std::runtime_error( "Ephemeris not defined for arrival body." );
//        }
//        else{
//            Eigen::Vector6d cartesianStateArrivalBody =
//                    bodies.at( departureAndArrivalBodies.at( 1 ) )->getEphemeris( )->getCartesianState(finalTime);
//            cartesianPositionAtArrivalForLambertTargeter =  cartesianStateArrivalBody.segment(0,3);
//        }
//    }
//    else
//    {
//        cartesianPositionAtArrivalForLambertTargeter = cartesianPositionAtArrival;
//    }


//    // Run the Lambert targeter.
//    mission_segments::LambertTargeterIzzo lambertTargeter(
//                cartesianPositionAtDepartureForLambertTargeter, cartesianPositionAtArrivalForLambertTargeter,
//                timeOfFlight, gravitationalParameterCentralBody );

//    // Retrieve cartesian state at departure.
//    Eigen::Vector3d cartesianVelocityAtDeparture = lambertTargeter.getInertialVelocityAtDeparture( );
//    Eigen::Vector6d cartesianStateAtDeparture;
//    cartesianStateAtDeparture.segment(0,3) = cartesianPositionAtDepartureForLambertTargeter;
//    cartesianStateAtDeparture.segment(3,3) = cartesianVelocityAtDeparture;


//    // Convert into keplerian elements
//    Eigen::Vector6d keplerianStateAtDeparture = orbital_element_conversions::convertCartesianToKeplerianElements(
//                cartesianStateAtDeparture, gravitationalParameterCentralBody );

//    // Propagate the keplerian elements until half of the time of flight.
//    Eigen::Vector6d keplerianStateAtHalvedTimeOfFlight = orbital_element_conversions::propagateKeplerOrbit( keplerianStateAtDeparture,
//                halvedTimeOfFlight, gravitationalParameterCentralBody );
//=======
//    // Get halved value of the time of flight, later used as initial time for the propagation.
//    double halvedTimeOfFlight = timeOfFlight / 2.0;

//    std::shared_ptr< mission_segments::LambertTargeterIzzo > lambertTargeter = createLambertTargeterIzzo(
//                timeOfFlight, initialTime, bodyMap, centralBody, departureAndArrivalBodies );
//    double centralBodyGravitationalParameter = bodyMap.at( centralBody )->getGravityFieldModel( )->getGravitationalParameter( );
//>>>>>>> dominic-origin/features/mission_segments_refactor

//    // Convert the keplerian elements back into Cartesian elements.
//    Eigen::Vector6d initialStatePropagationCartesianElements =
//            mission_segments::getLambertTargeterCartesianStateDuringTransfer(
//                *lambertTargeter.get( ), halvedTimeOfFlight );


//    // Define forward propagator settings variables.
//    integratorSettings->initialTime_ = initialTime + halvedTimeOfFlight;
    
//    // Define forward propagation settings
//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsForwardPropagation;
//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsBackwardPropagation;
    
//    propagatorSettingsForwardPropagation = propagatorSettings.second;
//    propagatorSettingsForwardPropagation->resetInitialStates( initialStatePropagationCartesianElements );
    
//    propagatorSettingsBackwardPropagation = propagatorSettings.first;
//    propagatorSettingsBackwardPropagation->resetInitialStates( initialStatePropagationCartesianElements );
    
//    // Perform forward propagation.
//    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards(
//                bodies, integratorSettings, propagatorSettingsForwardPropagation );
//    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.
//            getEquationsOfMotionNumericalSolution( );
//    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation =
//            dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );


//    // Calculate the difference between the full problem and the Lambert targeter solution along the forward propagation direction.
//    Eigen::Vector6d cartesianStateLambertSolution;
//    for( auto stateIterator : stateHistoryFullProblemForwardPropagation )
//    {
//        cartesianStateLambertSolution = orbital_element_conversions::propagateCartesianStateAlongKeplerOrbit(
//                    initialStatePropagationCartesianElements, stateIterator.first - ( initialTime + halvedTimeOfFlight ),
//                    centralBodyGravitationalParameter );
//        lambertTargeterResult[ stateIterator.first ] = cartesianStateLambertSolution;
//        fullProblemResult[ stateIterator.first ] = stateIterator.second;
//        dependentVariableResult[ stateIterator.first ] = dependentVariableHistoryForwardPropagation[ stateIterator.first ];
//    }


//    // Define backward propagator settings variables.
//    integratorSettings->initialTimeStep_ = -integratorSettings->initialTimeStep_;
//    integratorSettings->initialTime_ = initialTime + halvedTimeOfFlight;
    
//    // Perform the backward propagation.
//<<<<<<< HEAD
//    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards(bodies, integratorSettings, propagatorSettingsBackwardPropagation );
//=======
//    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards(
//                bodyMap, integratorSettings, propagatorSettingsBackwardPropagation );
//>>>>>>> dominic-origin/features/mission_segments_refactor
//    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation =
//            dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );
//    std::map< double, Eigen::VectorXd > dependentVariableHistoryBackwardPropagation =
//            dynamicsSimulatorIntegrationBackwards.getDependentVariableHistory( );
    
//    // Calculate the difference between the full problem and the Lambert targeter solution along the backward propagation direction.
//    for( auto stateIterator : stateHistoryFullProblemBackwardPropagation )
//    {
//        cartesianStateLambertSolution = orbital_element_conversions::propagateCartesianStateAlongKeplerOrbit(
//                    initialStatePropagationCartesianElements, - (initialTime + halvedTimeOfFlight) + stateIterator.first,
//                    centralBodyGravitationalParameter);
        
//        lambertTargeterResult[ stateIterator.first ] = cartesianStateLambertSolution;
//        fullProblemResult[ stateIterator.first ] = stateIterator.second;
//        dependentVariableResult[ stateIterator.first ] = dependentVariableHistoryBackwardPropagation[ stateIterator.first ];
        
//    }


//    integratorSettings->initialTimeStep_ = -integratorSettings->initialTimeStep_;
    
//}


//std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
//std::shared_ptr< propagators::PropagationTerminationSettings > > getLambertTargeterTerminationSettings(
//        const double timeOfFlight,
//        const double initialTime,
//<<<<<<< HEAD
//        const simulation_setup::SystemOfBodies& bodies,
//        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
//=======
//        const simulation_setup::NamedBodyMap& bodyMap,
//>>>>>>> dominic-origin/features/mission_segments_refactor
//        const std::string& bodyToPropagate,
//        const std::string& centralBody,
//        const std::pair< std::string, std::string >& departureAndArrivalBodies,
//        const bool setSphereOfInfluenceTermination,
//        const std::pair< double, double > distanceTerminationAsSoiFraction,
//        const double timeTerminationInSynodicPeriods )
//{
//    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
//            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettings;

//<<<<<<< HEAD
//    // Retrieve the gravitational parameter of the relevant bodies.
//    double gravitationalParameterCentralBody = ( centralBodyGravitationalParameter == centralBodyGravitationalParameter ) ?
//                centralBodyGravitationalParameter :
//                bodies.at( centralBody )->getGravityFieldModel( )->getGravitationalParameter( );


//    // Retrieve positions of departure and arrival bodies from ephemerides if not defined as inputs.
//    Eigen::Vector3d cartesianPositionAtDepartureForLambertTargeter, cartesianPositionAtArrivalForLambertTargeter;
//    if ( cartesianPositionAtDeparture != cartesianPositionAtDeparture )
//    {
//        // Cartesian state at departure
//        if ( bodies.at( departureAndArrivalBodies.at( 0 ) )->getEphemeris( ) == nullptr)
//=======
//    if (setSphereOfInfluenceTermination == true)
//    {
//        if( bodyMap.count( centralBody ) == 0 )
//>>>>>>> dominic-origin/features/mission_segments_refactor
//        {
//            throw std::runtime_error( "Error when getting Lambert targeter termination conditions, no central body " +
//                                      centralBody + " found" );
//        }
//        else if( bodyMap.at( centralBody )->getGravityFieldModel( ) == nullptr )
//        {
//<<<<<<< HEAD
//            Eigen::Vector6d cartesianStateDepartureBody =
//                    bodies.at( departureAndArrivalBodies.at( 0 ) )->getEphemeris( )->getCartesianState( initialTime);
//            cartesianPositionAtDepartureForLambertTargeter = cartesianStateDepartureBody.segment(0,3);
//=======
//            throw std::runtime_error( "Error when getting Lambert targeter termination conditions, central body " +
//                                      centralBody + " has no gravity field" );
//        }
//        else if( bodyMap.at( centralBody )->getEphemeris( ) == nullptr )
//        {
//            throw std::runtime_error( "Error when getting Lambert targeter termination conditions, central body " +
//                                      centralBody + " has no ephemeris" );
//>>>>>>> dominic-origin/features/mission_segments_refactor
//        }

//<<<<<<< HEAD
//        // Cartesian state at arrival
//        if ( bodies.at( departureAndArrivalBodies.at( 1 ) )->getEphemeris( ) == nullptr){
//            throw std::runtime_error( "Ephemeris not defined for arrival body." );
//        }
//        else{
//            Eigen::Vector6d cartesianStateArrivalBody =
//                    bodies.at( departureAndArrivalBodies.at( 1 ) )->getEphemeris( )->getCartesianState( initialTime + timeOfFlight );
//            cartesianPositionAtArrivalForLambertTargeter =  cartesianStateArrivalBody.segment(0,3);
//        }
//    }
//    else
//    {
//        cartesianPositionAtArrivalForLambertTargeter = cartesianPositionAtArrival;
//    }



//    // Calculate radii sphere of influence about departure and arrival bodies
//    double radiusSphereOfInfluenceDeparture;
//    double radiusSphereOfInfluenceArrival;

//    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
//            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettings;

//    if (terminationSphereOfInfluence == true)
//    {

//        double gravitationalParameterDepartureBody = ( departureBodyGravitationalParameter == departureBodyGravitationalParameter ) ?
//                    departureBodyGravitationalParameter :
//                    bodies.at( departureAndArrivalBodies[0] )->getGravityFieldModel( )->getGravitationalParameter( );
//        double gravitationalParameterArrivalBody = ( arrivalBodyGravitationalParameter == arrivalBodyGravitationalParameter ) ?
//                    arrivalBodyGravitationalParameter :
//                    bodies.at( departureAndArrivalBodies[1] )->getGravityFieldModel( )->getGravitationalParameter( );

//        double distanceDepartureToCentralBodies =
//                bodies.at( centralBody )->getEphemeris( )->getCartesianState(
//                    initialTime ).segment( 0, 3 ).norm( ) - cartesianPositionAtDepartureForLambertTargeter.segment( 0, 3 ).norm( );
//        double distanceArrivalToCentralBodies =
//                bodies.at( centralBody )->getEphemeris( )->getCartesianState(
//                    initialTime + timeOfFlight ).segment( 0, 3 ).norm( ) - cartesianPositionAtArrivalForLambertTargeter.segment( 0, 3 ).norm( );
//=======
//        if( bodyMap.count( departureAndArrivalBodies.first ) == 0 )
//        {
//            throw std::runtime_error( "Error when getting Lambert targeter termination conditions, no departure body " +
//                                      departureAndArrivalBodies.first + " found" );
//        }
//        else if( bodyMap.at( departureAndArrivalBodies.first )->getGravityFieldModel( ) == nullptr )
//        {
//            throw std::runtime_error( "Error when getting Lambert targeter termination conditions, departure body " +
//                                      departureAndArrivalBodies.first + " has no gravity field" );
//        }
//        else if( bodyMap.at( departureAndArrivalBodies.first )->getEphemeris( ) == nullptr )
//        {
//            throw std::runtime_error( "Error when getting Lambert targeter termination conditions, departure body " +
//                                      departureAndArrivalBodies.first + " has no ephemeris" );
//        }

//        if( bodyMap.count( departureAndArrivalBodies.second ) == 0 )
//        {
//            throw std::runtime_error( "Error when getting Lambert targeter termination conditions, no arrival body " +
//                                      departureAndArrivalBodies.second + " found" );
//        }
//        else if( bodyMap.at( departureAndArrivalBodies.second )->getGravityFieldModel( ) == nullptr )
//        {
//            throw std::runtime_error( "Error when getting Lambert targeter termination conditions, arrival body " +
//                                      departureAndArrivalBodies.second + " has no gravity field" );
//        }
//        else if( bodyMap.at( departureAndArrivalBodies.second )->getEphemeris( ) == nullptr )
//        {
//            throw std::runtime_error( "Error when getting Lambert targeter termination conditions, arrival body " +
//                                      departureAndArrivalBodies.second + " has no ephemeris" );
//        }
//>>>>>>> dominic-origin/features/mission_segments_refactor

//        double gravitationalParameterCentralBody =
//                bodyMap.at( centralBody )->getGravityFieldModel( )->getGravitationalParameter( );
//        double gravitationalParameterDepartureBody =
//                bodyMap.at( departureAndArrivalBodies.first )->getGravityFieldModel( )->getGravitationalParameter( );
//        double gravitationalParameterArrivalBody =
//                bodyMap.at( departureAndArrivalBodies.second )->getGravityFieldModel( )->getGravitationalParameter( );

//        double distanceDepartureToCentralBodY =
//                ( bodyMap.at( centralBody )->getEphemeris( )->getCartesianState( initialTime ).segment( 0, 3 ) -
//                  bodyMap.at( departureAndArrivalBodies.first )->getEphemeris( )->getCartesianState( initialTime ).segment( 0, 3 ) ).norm( );
//        double distanceArrivalToCentralBodY =
//                ( bodyMap.at( centralBody )->getEphemeris( )->getCartesianState( initialTime + timeOfFlight ).segment( 0, 3 ) -
//                  bodyMap.at( departureAndArrivalBodies.second )->getEphemeris( )->getCartesianState( initialTime + timeOfFlight ).segment( 0, 3 ) ).norm( );

//        // Calculate radius sphere of influence for departure body.
//        double radiusSphereOfInfluenceDeparture = tudat::mission_geometry::computeSphereOfInfluence(
//                    distanceDepartureToCentralBodY, gravitationalParameterDepartureBody, gravitationalParameterCentralBody);

//        // Calculate radius sphere of influence for arrival body.
//        double radiusSphereOfInfluenceArrival = tudat::mission_geometry::computeSphereOfInfluence(
//                    distanceArrivalToCentralBodY, gravitationalParameterArrivalBody, gravitationalParameterCentralBody);

//        // Calculate the synodic period.
//        double orbitalPeriodDepartureBody = basic_astrodynamics::computeKeplerOrbitalPeriod(
//<<<<<<< HEAD
//              orbital_element_conversions::convertCartesianToKeplerianElements( bodies.at( departureAndArrivalBodies[0] )->
//              getEphemeris()->getCartesianState(initialTime), gravitationalParameterCentralBody)[ orbital_element_conversions::semiMajorAxisIndex ],
//              gravitationalParameterCentralBody, gravitationalParameterDepartureBody);

//        double orbitalPeriodArrivalBody = basic_astrodynamics::computeKeplerOrbitalPeriod(
//              orbital_element_conversions::convertCartesianToKeplerianElements( bodies.at( departureAndArrivalBodies[1] )->
//              getEphemeris()->getCartesianState(initialTime), gravitationalParameterCentralBody)[ orbital_element_conversions::semiMajorAxisIndex ],
//              gravitationalParameterCentralBody, gravitationalParameterArrivalBody);
//=======
//                    orbital_element_conversions::convertCartesianToKeplerianElements(
//                        bodyMap.at( departureAndArrivalBodies.first )->getEphemeris()->getCartesianState(initialTime),
//                        gravitationalParameterCentralBody)[ orbital_element_conversions::semiMajorAxisIndex ],
//                gravitationalParameterCentralBody, gravitationalParameterDepartureBody);

//        double orbitalPeriodArrivalBody = basic_astrodynamics::computeKeplerOrbitalPeriod(
//                    orbital_element_conversions::convertCartesianToKeplerianElements(
//                        bodyMap.at( departureAndArrivalBodies.second )-> getEphemeris()->getCartesianState(initialTime),
//                        gravitationalParameterCentralBody)[ orbital_element_conversions::semiMajorAxisIndex ],
//                gravitationalParameterCentralBody, gravitationalParameterArrivalBody);
//>>>>>>> dominic-origin/features/mission_segments_refactor

//        double synodicPeriod;
//        if (orbitalPeriodDepartureBody < orbitalPeriodArrivalBody)
//        {
//            synodicPeriod = basic_astrodynamics::computeSynodicPeriod(orbitalPeriodDepartureBody, orbitalPeriodArrivalBody);
//        }
//        else
//        {
//            synodicPeriod = basic_astrodynamics::computeSynodicPeriod(orbitalPeriodArrivalBody, orbitalPeriodDepartureBody);
//        }


//        // Create total propagator termination settings.
//        std::vector< std::shared_ptr< PropagationTerminationSettings > >  forwardPropagationTerminationSettingsList;
//        forwardPropagationTerminationSettingsList.push_back(
//                    std::make_shared< PropagationDependentVariableTerminationSettings >(
//                        std::make_shared< SingleDependentVariableSaveSettings >(
//                            relative_distance_dependent_variable, bodyToPropagate, departureAndArrivalBodies.second ),
//                        distanceTerminationAsSoiFraction.second * radiusSphereOfInfluenceArrival, false ) );
//        forwardPropagationTerminationSettingsList.push_back(
//                    std::make_shared< PropagationTimeTerminationSettings >( timeTerminationInSynodicPeriods * synodicPeriod ) );


//        std::shared_ptr< PropagationTerminationSettings > forwardPropagationTerminationSettings =
//                std::make_shared< PropagationHybridTerminationSettings >( forwardPropagationTerminationSettingsList, true );


//        std::vector< std::shared_ptr< PropagationTerminationSettings > >  backwardPropagationTerminationSettingsList;
//        backwardPropagationTerminationSettingsList.push_back(
//                    std::make_shared< PropagationDependentVariableTerminationSettings >(
//                        std::make_shared< SingleDependentVariableSaveSettings >(
//                            relative_distance_dependent_variable, bodyToPropagate, departureAndArrivalBodies.first ),
//                        distanceTerminationAsSoiFraction.first * radiusSphereOfInfluenceDeparture, false ) );
//        backwardPropagationTerminationSettingsList.push_back(
//                    std::make_shared< PropagationTimeTerminationSettings >( timeTerminationInSynodicPeriods * synodicPeriod ) );

//        \
//        std::shared_ptr< PropagationTerminationSettings > backwardPropagationTerminationSettings =
//                std::make_shared< PropagationHybridTerminationSettings >( backwardPropagationTerminationSettingsList, true );

//        terminationSettings = std::make_pair( backwardPropagationTerminationSettings, forwardPropagationTerminationSettings );


//    }
//    else
//    {
//        terminationSettings = std::make_pair(
//                    std::make_shared< propagators::PropagationTimeTerminationSettings >( initialTime ),
//                    std::make_shared< propagators::PropagationTimeTerminationSettings >( initialTime + timeOfFlight ) );
//    }
//    return terminationSettings;

//<<<<<<< HEAD
//    Eigen::Vector6d initialState;

//    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings;

//    std::vector< std::string > centralBodyPropagation; centralBodyPropagation.push_back( centralBody );
//    std::vector< std::string > bodyToPropagatePropagation; bodyToPropagatePropagation.push_back( bodyToPropagate );

//    propagatorSettings.first = std::make_shared< TranslationalStatePropagatorSettings< double > >(
//        centralBodyPropagation, accelerationModelMap, bodyToPropagatePropagation, initialState,
//        terminationSettings.first, propagator, dependentVariablesToSave );

//    propagatorSettings.second = std::make_shared< TranslationalStatePropagatorSettings< double > >(
//        centralBodyPropagation, accelerationModelMap, bodyToPropagatePropagation, initialState,
//        terminationSettings.second, propagator, dependentVariablesToSave );

//    propagateLambertTargeterAndFullProblem(
//            timeOfFlight, initialTime, bodies, centralBody, propagatorSettings,
//            integratorSettings, lambertTargeterResult, fullProblemResult, dependentVariableResult, departureAndArrivalBodies,
//            centralBodyGravitationalParameter, cartesianPositionAtDeparture, cartesianPositionAtArrival );
//=======
//>>>>>>> dominic-origin/features/mission_segments_refactor
//}
////! Function to propagate the full dynamics problem and the Lambert targeter solution.
//void propagateLambertTargeterAndFullProblem(
//        const double timeOfFlight,
//        const double initialTime,
//        const simulation_setup::SystemOfBodies& bodies,
//        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
//        const std::string& bodyToPropagate,
//        const std::string& centralBody,
//        const std::pair< std::string, std::string >& departureAndArrivalBodies,
//        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
//        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
//        std::map< double, Eigen::Vector6d >& fullProblemResult,
//        std::map< double, Eigen::VectorXd >& dependentVariableResult,
//        const bool terminationSphereOfInfluence,
//        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave ,
//        const TranslationalPropagatorType propagator )
//{

//    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
//            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettings =
//            getLambertTargeterTerminationSettings( timeOfFlight, initialTime,
//                bodyMap, bodyToPropagate, centralBody, departureAndArrivalBodies, terminationSphereOfInfluence );

//    Eigen::Vector6d initialState;
    
//    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
//            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings;
    
//    std::vector< std::string > centralBodyPropagation; centralBodyPropagation.push_back( centralBody );
//    std::vector< std::string > bodyToPropagatePropagation; bodyToPropagatePropagation.push_back( bodyToPropagate );
    
//    propagatorSettings.first = std::make_shared< TranslationalStatePropagatorSettings< double > >(
//                centralBodyPropagation, accelerationModelMap, bodyToPropagatePropagation, initialState,
//                terminationSettings.first, propagator, dependentVariablesToSave );
    
//    propagatorSettings.second = std::make_shared< TranslationalStatePropagatorSettings< double > >(
//                centralBodyPropagation, accelerationModelMap, bodyToPropagatePropagation, initialState,
//                terminationSettings.second, propagator, dependentVariablesToSave );
    
//    propagateLambertTargeterAndFullProblem(
//<<<<<<< HEAD
//                timeOfFlight, initialTime, bodies, accelerationModelMap, bodyToPropagate, centralBody, integratorSettings,
//                lambertTargeterResult, fullProblemResult, dependentVariableResult, departureAndArrivalBodies,
//                terminationSphereOfInfluence, cartesianPositionAtDeparture, cartesianPositionAtArrival, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN,
//                dependentVariablesToSave, propagator );

//    Eigen::Vector6d stateLambertTargeterAtDeparture = lambertTargeterResult.begin( )->second;
//    Eigen::Vector6d propagatedStateFullProblemAtDeparture = fullProblemResult.begin( )->second;
//    Eigen::Vector6d stateLambertTargeterAtArrival = lambertTargeterResult.rbegin( )->second;
//    Eigen::Vector6d propagatedStateFullProblemAtArrival = fullProblemResult.rbegin( )->second;


//    // Difference between the Lambert targeter and full problem results at departure and arrival.
//    return std::make_pair( stateLambertTargeterAtDeparture - propagatedStateFullProblemAtDeparture,
//                           stateLambertTargeterAtArrival - propagatedStateFullProblemAtArrival);
//=======
//                timeOfFlight, initialTime, bodyMap, centralBody, departureAndArrivalBodies, propagatorSettings,
//                integratorSettings, lambertTargeterResult, fullProblemResult, dependentVariableResult );
//>>>>>>> dominic-origin/features/mission_segments_refactor
//}


//}

//}
