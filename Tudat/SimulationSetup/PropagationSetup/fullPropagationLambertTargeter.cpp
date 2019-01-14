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

namespace tudat
{

namespace propagators
{

//! Function to directly setup a body map corresponding to the assumptions of the Lambert targeter.
simulation_setup::NamedBodyMap setupBodyMapLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate)
{

    spice_interface::loadStandardSpiceKernels( );


    // Create central body object.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( nameCentralBody );
    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
                    simulation_setup::getDefaultBodySettings( bodiesToCreate );


    // Define central body ephemeris settings.
    std::string frameOrigin = "SSB";
    std::string frameOrientation = "J2000";
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

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );



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

    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[nameCentralBody].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                          basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ nameBodyToPropagate ] = bodyToPropagateAccelerations;

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    return accelerationModelMap;

}


//! Function to propagate the Lambert targeter solution from the associated environment




//! Function to determine the cartesian state at a given time for a keplerian orbit, based on the initial state.
Eigen::Vector6d propagateLambertTargeterSolution(
        const Eigen::Vector6d& initialState,
        const double finalPropagationTime,
        const double gravitationalParameter)
{

    Eigen::Vector6d keplerianInitialState = orbital_element_conversions::convertCartesianToKeplerianElements(initialState,
                                                                                                             gravitationalParameter);
    // Retrieve the semi-major axis and eccentricty of the keplerian orbit.
    double semiMajorAxis = keplerianInitialState[orbital_element_conversions::semiMajorAxisIndex];
    double eccentricity = keplerianInitialState[orbital_element_conversions::eccentricityIndex];

    // Calculate the initial mean anomaly.
    double initialTrueAnomaly = keplerianInitialState[orbital_element_conversions::trueAnomalyIndex];
    double initialMeanAnomaly = orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
                orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly(initialTrueAnomaly, eccentricity), eccentricity);

    // Calculate the mean anomaly at the final time.
    double meanAnomalyEndPropagation = initialMeanAnomaly + orbital_element_conversions::convertElapsedTimeToMeanAnomalyChange(
                finalPropagationTime, gravitationalParameter, semiMajorAxis);

    // Determine the final
    Eigen::Vector6d finalKeplerianState = keplerianInitialState;
    finalKeplerianState[orbital_element_conversions::trueAnomalyIndex] =
            orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly(eccentricity, meanAnomalyEndPropagation), eccentricity);

    Eigen::Vector6d cartesianStateLambertSolution = orbital_element_conversions::convertKeplerianToCartesianElements(
                finalKeplerianState, gravitationalParameter);

    return cartesianStateLambertSolution;
}

//! Function to propagate the full dynamics problem and the Lambert targeter solution.
void propagateLambertTargeterAndFullProblem( Eigen::Vector3d cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector<std::string>& bodiesToPropagate,
        const std::vector<std::string>& centralBody,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        std::map< double, Eigen::Vector6d >& fullProblemResult,
        const std::vector<std::string>& arrivalAndDepartureBodies,
        const bool arrivalAndDepartureInitialisationFromEphemerides = false,
        const bool terminationSettings = false)
{

    lambertTargeterResult.clear( );
    fullProblemResult.clear( );

    // Retrieve the gravitational parameter of the main body.
    double gravitationalParameter = bodyMap[centralBody[0]]->getGravityFieldModel()->getGravitationalParameter();

    // Get halved value of the time of flight, used as initial time for the propagation.
    double halvedTimeOfFlight = timeOfFlight / 2.0;



    // Retrieve positions from ephemerides
    double radiusSphereOfInfluenceDeparture = 1000;
    double radiusSphereOfInfluenceArrival = 1000;

    if (arrivalAndDepartureInitialisationFromEphemerides == true)
    {
//        std::vector< std::string > bodiesWithEphemerisData;
//        bodiesWithEphemerisData.push_back("mercury");
//        bodiesWithEphemerisData.push_back("venus");
//        bodiesWithEphemerisData.push_back("earthMoonBarycenter");
//        bodiesWithEphemerisData.push_back("mars");
//        bodiesWithEphemerisData.push_back("jupiter");
//        bodiesWithEphemerisData.push_back("saturn");
//        bodiesWithEphemerisData.push_back("uranus");
//        bodiesWithEphemerisData.push_back("neptune");
//        bodiesWithEphemerisData.push_back("pluto");


        std::map< std::string, std::shared_ptr< tudat::simulation_setup::BodySettings > > bodySettings =
                  tudat::simulation_setup::getDefaultBodySettings( arrivalAndDepartureBodies );
        bodySettings[ arrivalAndDepartureBodies[0] ]->ephemerisSettings = std::make_shared<
                simulation_setup::DirectSpiceEphemerisSettings >("SSB", "ECLIPJ2000");
        bodySettings[ arrivalAndDepartureBodies[1] ]->ephemerisSettings = std::make_shared<
                simulation_setup::DirectSpiceEphemerisSettings> ("SSB", "ECLIPJ2000");

//        bodySettings[arrivalAndDepartureBodies[0]]->ephemerisSettings->
//                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), "SSB", "J2000" );

        simulation_setup::NamedBodyMap bodyMapArrivalAndDepartureBodies = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMapArrivalAndDepartureBodies, "SSB", "ECLIPJ2000" );

        Eigen::Vector6d terminationStateArrivalBody =
                bodyMapArrivalAndDepartureBodies[ arrivalAndDepartureBodies[0] ]->getEphemeris()->getCartesianState(timeOfFlight);
        Eigen::Vector6d terminationStateDepartureBody =
                bodyMapArrivalAndDepartureBodies[ arrivalAndDepartureBodies[1] ]->getEphemeris()->getCartesianState(timeOfFlight);

        std::cout << "terminationStateArrivalBody: " << terminationStateArrivalBody << "\n\n";
        std::cout << "terminationStateArrivalBody: " << terminationStateDepartureBody << "\n\n";


        // If we stop the propagation at the sphere of influence

        if (terminationSettings == true) {

            double distanceToCentralBody = ( bodyMap[ centralBody[0] ]->getState() - terminationStateArrivalBody ).norm();
            double massDepartureBody = bodyMapArrivalAndDepartureBodies[ arrivalAndDepartureBodies[0] ]
                    ->getGravityFieldModel()->getGravitationalParameter();
            double massArrivalBody = bodyMapArrivalAndDepartureBodies[ arrivalAndDepartureBodies[1] ]
                    ->getGravityFieldModel()->getGravitationalParameter();
            double massCentralBody = bodyMap[ centralBody[0] ]->getGravityFieldModel()->getGravitationalParameter();
            radiusSphereOfInfluenceDeparture = tudat::mission_geometry::computeSphereOfInfluence(
                        distanceToCentralBody, massDepartureBody, massCentralBody);
            radiusSphereOfInfluenceArrival = tudat::mission_geometry::computeSphereOfInfluence(
                        distanceToCentralBody, massArrivalBody, massCentralBody);


        }

    }



    // Run the Lambert targeter.
    mission_segments::LambertTargeterIzzo LambertTargeter(
                cartesianPositionAtDeparture, cartesianPositionAtArrival, timeOfFlight, gravitationalParameter );

    // Retrieve cartesian state at departure.
    Eigen::Vector3d cartesianVelocityAtDeparture = LambertTargeter.getInertialVelocityAtDeparture();
    Eigen::Vector6d cartesianStateAtDeparture;
    cartesianStateAtDeparture.segment(0,3) = cartesianPositionAtDeparture;
    cartesianStateAtDeparture.segment(3,3) = cartesianVelocityAtDeparture;

    // Keplerian state at departure.
    Eigen::Vector6d keplerianElementsAtDeparture = tudat::orbital_element_conversions::convertCartesianToKeplerianElements(
                cartesianStateAtDeparture, gravitationalParameter);

    double semiMajorAxis = LambertTargeter.getSemiMajorAxis();
    double eccentricity = keplerianElementsAtDeparture( orbital_element_conversions::eccentricityIndex );

    double trueAnomalyAtDeparture = keplerianElementsAtDeparture(orbital_element_conversions::trueAnomalyIndex);
    double meanAnomalyAtDeparture = orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
                orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly(trueAnomalyAtDeparture, eccentricity),
                eccentricity);


    // Calculate the true anomaly at half the time of flight.
    double meanAnomalyChangeHalfTimeOfFlight = orbital_element_conversions::convertElapsedTimeToMeanAnomalyChange(halvedTimeOfFlight,
                                                                                      gravitationalParameter, semiMajorAxis);

    double meanAnomalyHalfTimeOfFlight = meanAnomalyChangeHalfTimeOfFlight + meanAnomalyAtDeparture;
    double trueAnomalyHalfTimeOfFlight = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( eccentricity, meanAnomalyHalfTimeOfFlight ), eccentricity);


    // Define the state at half of the time of flight (initial state for the propagation).
    Eigen::Vector6d initialStatePropagationKeplerianElements;
    initialStatePropagationKeplerianElements.segment(0,5) = keplerianElementsAtDeparture.segment(0,5);
    initialStatePropagationKeplerianElements[orbital_element_conversions::trueAnomalyIndex] = trueAnomalyHalfTimeOfFlight;

    Eigen::Vector6d initialStatePropagationCartesianElements = orbital_element_conversions::convertKeplerianToCartesianElements(
                initialStatePropagationKeplerianElements, gravitationalParameter);



    // Initialise variables for propagatation.
    std::vector< std::string > centralBodiesPropagation;
    centralBodiesPropagation.push_back( "SSB" );

    Eigen::Vector6d cartesianStateLambertSolution;




    // Define forward propagator settings variables.
    integratorSettings->initialTime_ = 0.0;

    // Termination settings if the forward propagation stops at the sphere of influence of the arrival body
    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariableAtArrival =
          std::make_shared< SingleDependentVariableSaveSettings >(
              relative_distance_dependent_variable, bodiesToPropagate[0], arrivalAndDepartureBodies[1] );
    std::shared_ptr< PropagationTerminationSettings > forwardPropagationTerminationSettings =
            std::make_shared< PropagationDependentVariableTerminationSettings >( terminationDependentVariableAtArrival,
                                                                                 radiusSphereOfInfluenceArrival, true);

    // Define forward propagation settings
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsForwardPropagation;

//    if (terminationSettings == false){
        propagatorSettingsForwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                ( centralBodiesPropagation, accelerationModelMap, bodiesToPropagate, initialStatePropagationCartesianElements,
                                                                                                         halvedTimeOfFlight );
//    }
//    else {
//        propagatorSettingsForwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
//                ( centralBodiesPropagation, accelerationModelMap, bodiesToPropagate, initialStatePropagationCartesianElements,
//                  forwardPropagationTerminationSettings, propagators::cowell );

//    }


    // Perform forward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards(bodyMap, integratorSettings, propagatorSettingsForwardPropagation );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.getEquationsOfMotionNumericalSolution( );

    // Calculate the difference between the full problem and the Lambert targeter solution along the forward propagation branch of the orbit.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {

        cartesianStateLambertSolution = propagateLambertTargeterSolution(initialStatePropagationCartesianElements, itr->first,
                                                                         gravitationalParameter);

        lambertTargeterResult[ itr->first ] = cartesianStateLambertSolution;
        fullProblemResult[ itr->first ] = itr->second;

    }





    // Define backward propagator settings variables.
    integratorSettings->initialTimeStep_ = -1 * integratorSettings->initialTimeStep_;
    integratorSettings->initialTime_ = 0.0;

    // Termination settings if the backward propagation stops at the sphere of influence of the departure body
    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariableAtDeparture =
          std::make_shared< SingleDependentVariableSaveSettings >(
              relative_distance_dependent_variable, bodiesToPropagate[0], arrivalAndDepartureBodies[0] );
    std::shared_ptr< PropagationTerminationSettings > backwardPropagationTerminationSettings =
            std::make_shared< PropagationDependentVariableTerminationSettings >( terminationDependentVariableAtDeparture,
                                                                                 radiusSphereOfInfluenceDeparture, true);

    // Define backward propagation settings.
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsBackwardPropagation;
    if (terminationSettings == false){
        propagatorSettingsBackwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                ( centralBodiesPropagation, accelerationModelMap, bodiesToPropagate, initialStatePropagationCartesianElements,
                                                                                                  -1.0 * halvedTimeOfFlight );
    }
    else {
        propagatorSettingsBackwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                ( centralBodiesPropagation, accelerationModelMap, bodiesToPropagate, initialStatePropagationCartesianElements,
                  backwardPropagationTerminationSettings, propagators::cowell );
    }


    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards(bodyMap, integratorSettings,
                                                                                     propagatorSettingsBackwardPropagation );

    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );

    // Calculate the difference between the full problem and the Lambert targeter solution along the forward propagation branch of the orbit.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
    {

        cartesianStateLambertSolution = propagateLambertTargeterSolution(initialStatePropagationCartesianElements, itr->first,
                                                                         gravitationalParameter);

        lambertTargeterResult[ itr->first ] = cartesianStateLambertSolution;
        fullProblemResult[ itr->first ] = itr->second;

    }

}


//! Function to compute the difference in cartesian state between Lambert targeter solution and full dynamics problem, both at departure
//! and at arrival.
std::pair< Eigen::Vector6d, Eigen::Vector6d > getDifferenceFullPropagationWrtLambertTargeterAtDepartureAndArrival(
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::vector< std::string >& departureAndArrivalBodies)

{
    std::map< double, Eigen::Vector6d > lambertTargeterResult;
    std::map< double, Eigen::Vector6d > fullProblemResult;

    // compute full problem and Lambert targeter solution both at departure and arrival.
    propagateLambertTargeterAndFullProblem(cartesianPositionAtDeparture, cartesianPositionAtArrival, timeOfFlight,
                                           bodyMap, accelerationModelMap, bodiesToPropagate, centralBodies, integratorSettings,
                                           lambertTargeterResult, fullProblemResult, departureAndArrivalBodies, true, true);

    Eigen::Vector6d stateLambertTargeterAtDeparture = lambertTargeterResult.begin( )->second;
    Eigen::Vector6d propagatedStateFullProblemAtDeparture = fullProblemResult.begin( )->second;
    Eigen::Vector6d stateLambertTargeterAtArrival = lambertTargeterResult.rbegin( )->second;
    Eigen::Vector6d propagatedStateFullProblemAtArrival = fullProblemResult.rbegin( )->second;

    // Difference between the two propagated states at departure and arrival.
    return std::make_pair( stateLambertTargeterAtDeparture - propagatedStateFullProblemAtDeparture,
                           stateLambertTargeterAtArrival - propagatedStateFullProblemAtArrival);
}



}

}
