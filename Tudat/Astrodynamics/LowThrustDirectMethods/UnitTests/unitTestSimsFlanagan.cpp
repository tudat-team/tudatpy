/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wakker, K. F. (2007), Lecture Notes Astrodynamics II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <math.h>
#include <iostream>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
//#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanagan.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganLeg.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganOptimisationSetup.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "pagmo/algorithms/de1220.hpp"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"

namespace tudat
{
namespace unit_tests
{

//! Test Sims Flanagan implementation.
BOOST_AUTO_TEST_SUITE( test_Sims_Flanagan )


//! Limit case: if the maximum thrust is set to 0, the Sims Flanagan trajectory should be a keplerian one.
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_limit_case )
{
    using namespace low_thrust_direct_methods;

    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.0;
    double specificImpulse = 3000.0;
    double mass = 1800.0;
    int numberSegments = 10;

    // Define (constant) specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return specificImpulse;
    };

    double julianDate = 1000.0 * physical_constants::JULIAN_DAY; //2458849.5;
    double timeOfFlight = 700.0 * physical_constants::JULIAN_DAY;

    std::string bodyToPropagate = "Vehicle";
    std::string centralBody = "Sun";

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ centralBody ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ centralBody ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ centralBody ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ bodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( bodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );


    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Define state at departure and arrival.
    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    // Define integrator settings.
    double stepSize = timeOfFlight / 500.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize );


    // Define thrust throttles.
    std::vector< Eigen::Vector3d > throttles;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        throttles.push_back( ( Eigen::Vector3d( ) << 0.3, 0.3, 0.3 ).finished( ) );
    }


    // Re-initialise mass of the spacecraft in the body map.
    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );
    double centralBodyGravitationalParameter = bodyMap[ centralBody ]->getGravityFieldModel()->getGravitationalParameter();

    SimsFlanaganLeg simsFlanaganLeg = SimsFlanaganLeg( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction,
                                                       timeOfFlight, bodyMap, throttles, bodyToPropagate, centralBody );

    simsFlanaganLeg.propagateForwardFromDepartureToMatchPoint( );
    std::cout << "state at match point from forward propagation: " << simsFlanaganLeg.getStateAtMatchPointForwardPropagation( ).transpose() << "\n\n";
    std::cout << "state forward propagation keplerian orbit: " << orbital_element_conversions::convertKeplerianToCartesianElements(
                     orbital_element_conversions::propagateKeplerOrbit(
                    orbital_element_conversions::convertCartesianToKeplerianElements( stateAtDeparture, centralBodyGravitationalParameter),
                     timeOfFlight / 2.0, centralBodyGravitationalParameter ), centralBodyGravitationalParameter ).transpose() << "\n\n";
    simsFlanaganLeg.propagateBackwardFromArrivalToMatchPoint( );
    std::cout << "state at match point from backward propagation: " << simsFlanaganLeg.getStateAtMatchPointBackwardPropagation( ).transpose() << "\n\n";
    std::cout << "state backward propagation keplerian orbit: " << orbital_element_conversions::convertKeplerianToCartesianElements(
                     orbital_element_conversions::propagateKeplerOrbit(
                    orbital_element_conversions::convertCartesianToKeplerianElements( stateAtArrival, centralBodyGravitationalParameter),
                     - timeOfFlight / 2.0, centralBodyGravitationalParameter ), centralBodyGravitationalParameter ).transpose() << "\n\n";


    /// TEST PROPAGATE SIMS FLANAGAN SOLUTION TO GIVEN TIME.
    maximumThrust = 0.8;
    simsFlanaganLeg = SimsFlanaganLeg( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction,
                                                           timeOfFlight, bodyMap, throttles, bodyToPropagate, centralBody );

    simsFlanaganLeg.propagateForwardFromDepartureToMatchPoint( );
    simsFlanaganLeg.propagateBackwardFromArrivalToMatchPoint( );

    int numberSegmentsForwardPropagation = ( numberSegments + 1 ) / 2;
    int numberSegmentsBackwardPropagation = numberSegments / 2;
    double segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    double segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );
    std::vector< double > timesAtNodes;

    Eigen::Vector6d currentState = stateAtDeparture;
    for ( int i = 0 ; i <= numberSegmentsForwardPropagation ; i++ )
    {
        timesAtNodes.push_back( i * segmentDurationForwardPropagation );
    }
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        currentState = simsFlanaganLeg.propagateInsideForwardSegment( timesAtNodes[ i ], timesAtNodes[ i + 1 ], segmentDurationForwardPropagation,
                currentState );
    }
    std::cout << "state propagated segment-wise (forward): " << currentState.transpose( ) << "\n\n";
    std::cout << "comparison (forward propagation): " << simsFlanaganLeg.getStateAtMatchPointForwardPropagation( ).transpose() << "\n\n";

    currentState = stateAtArrival;
    timesAtNodes.clear( );
    for ( int i = 0 ; i <= numberSegmentsBackwardPropagation ; i++ )
    {
        timesAtNodes.push_back( timeOfFlight / 2.0 + i * segmentDurationBackwardPropagation );
    }


    for ( int i = timesAtNodes.size( ) - 1 ; i > 0 ; i-- )
    {
        currentState = simsFlanaganLeg.propagateInsideBackwardSegment( timesAtNodes[ i ], timesAtNodes[ i - 1 ], segmentDurationBackwardPropagation,
                currentState );
    }
    std::cout << "state propagated segment-wise (backward): " << currentState.transpose( ) << "\n\n";
    std::cout << "comparison (backward propagation): " << simsFlanaganLeg.getStateAtMatchPointBackwardPropagation( ).transpose() << "\n\n";

}


//! Test Sims-Flanagan implementation by comparing it with trajectory obtained with successive impulsive shots applied at times corresponding to
//! half of each of the Sims-Flanagan segments.
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_implementation )
{
    using namespace low_thrust_direct_methods;

    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.8;
    double specificImpulse = 3000.0;
    double mass = 1800.0;
    int numberSegments = 10;

    // Define (constant) specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return specificImpulse;
    };

    double julianDate = 1000.0 * physical_constants::JULIAN_DAY; //2458849.5;
    double timeOfFlight = 700.0 * physical_constants::JULIAN_DAY;

    std::string bodyToPropagate = "Vehicle";
    std::string centralBody = "Sun";


    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ centralBody ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ centralBody ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ centralBody ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ bodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( bodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );


    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Define state at departure and arrival.
    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    // Define integrator settings.
    double stepSize = ( timeOfFlight ) / static_cast< double >( 50 );
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize / 10.0 );



    // Define optimisation algorithm.
    algorithm optimisationAlgorithm{ pagmo::de1220() };

    SimsFlanagan simsFlanagan = SimsFlanagan( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction, numberSegments,
                                              timeOfFlight, bodyMap, bodyToPropagate, centralBody, optimisationAlgorithm, 1, 10 );

    std::vector< double > fitnessVector = simsFlanagan.getBestIndividualFitness( );
    std::vector< double > bestIndividual = simsFlanagan.getBestIndividual( );

    std::vector< Eigen::Vector3d > bestThrottles;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        bestThrottles.push_back( ( Eigen::Vector3d( ) << bestIndividual[ i * 3 ], bestIndividual[ i * 3 + 1 ],
                bestIndividual[ i * 3 + 2 ] ).finished( ) );
    }

    ///! TEST SIMS FLANAGAN FULL PROPAGATION

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, std::vector< std::string >{ bodyToPropagate },
                std::vector< std::string >{ centralBody } );

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight, true );


    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Define pair of propagatorSettings for backward and forward propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings;

    propagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate }, stateAtDeparture,
                          terminationConditions.first, propagators::cowell, dependentVariablesToSave );

    propagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate }, stateAtArrival,
                          terminationConditions.second, propagators::cowell, dependentVariablesToSave );

    // Define empty maps to store the propagation results.
    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::Vector6d > simsFlanaganResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    // Compute full propagation.
    simsFlanagan.computeSimsFlanaganTrajectoryAndFullPropagation( integratorSettings, propagatorSettings, fullPropagationResults,
                                                                  simsFlanaganResults, dependentVariablesHistory );


    //! Test full propagation w.r.t. impulsive shots as thrust acceleration.

    // Calculate number of segments for both the forward propagation (from departure to match point)
    // and the backward propagation (from arrival to match point).
    int numberSegmentsForwardPropagation = ( numberSegments + 1 ) / 2;
    int numberSegmentsBackwardPropagation = numberSegments / 2;
    int segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    int segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );

    // Compute times at half of each segment.
    std::vector< double > thrustMidTimes;
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        thrustMidTimes.push_back( segmentDurationForwardPropagation / 2.0 + i * segmentDurationForwardPropagation );
    }
    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        thrustMidTimes.push_back( segmentDurationBackwardPropagation / 2.0 + timeOfFlight / 2.0 + i * segmentDurationBackwardPropagation );
    }


    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );
    SimsFlanaganLeg simsFlanaganLegTest = SimsFlanaganLeg( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction,
                                                           timeOfFlight, bodyMap, bestThrottles, bodyToPropagate, centralBody );

    // Compute state at half of the time of flight.
    simsFlanaganLegTest.propagateForwardFromDepartureToMatchPoint( );
    simsFlanaganLegTest.propagateBackwardFromArrivalToMatchPoint( );
    Eigen::Vector6d stateAtHalfTimeOfFlight = simsFlanaganLegTest.getStateAtMatchPointForwardPropagation( );
    Eigen::Vector6d stateAtHalfTimeOfFlightBackwardPropagation = simsFlanaganLegTest.getStateAtMatchPointBackwardPropagation( );
    std::cout << "state half TOF test unit test: " << stateAtHalfTimeOfFlight.transpose() << "\n\n";
    std::cout << "state half TOF test unit test: " << stateAtHalfTimeOfFlightBackwardPropagation.transpose() << "\n\n";

    // Compute mass at half of the time of flight.
    double massAtHalfTimeOfFlight = simsFlanaganLegTest.getMassAtMatchPointForwardPropagation( );

    // Compute deltaVs.
    double totalManeuverTime = 90.0;
    double maneuverRiseTime = 15.0;
    double currentMass = mass;
    std::vector< Eigen::Vector3d > deltaVs;

    // Compute deltaVs for the forward propagation half.
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        Eigen::Vector3d currentDeltaVvector = maximumThrust * bestThrottles[ i ] * segmentDurationForwardPropagation / currentMass;
        deltaVs.push_back( currentDeltaVvector );

        // Update mass.
        currentMass *= std::exp( - currentDeltaVvector.norm() /
                                 ( specificImpulseFunction( 0.0 ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );
    }
    // Compute deltaVs for the backward propagation half.
    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        Eigen::Vector3d currentDeltaVvector = maximumThrust * bestThrottles[ i + numberSegmentsForwardPropagation ] *
                segmentDurationBackwardPropagation / currentMass;
        deltaVs.push_back( currentDeltaVvector );

        // Update mass.
        currentMass *= std::exp( - currentDeltaVvector.norm() /
                                 ( specificImpulseFunction( 0.0 ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );
    }
    std::cout << "current mass: " << currentMass << "\n\n";


    bodyToPropagateAccelerations.clear();
    bodyToPropagateAccelerations[ centralBody ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ bodyToPropagate ].push_back( std::make_shared< simulation_setup::MomentumWheelDesaturationAccelerationSettings >(
                                                                   thrustMidTimes, deltaVs, totalManeuverTime, maneuverRiseTime ) );

    accelerationMap.clear();
    accelerationMap[ bodyToPropagate ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    accelerationModelMap.clear();
    accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, std::vector< std::string >{ bodyToPropagate },
                std::vector< std::string >{ centralBody } );

    // BACKWARD PROPAGATION.

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( massAtHalfTimeOfFlight );

    // Create termination conditions settings.
    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings
            = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );

    //  Create propagator settings.
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsImpulsiveDeltaV =
            std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                std::vector< std::string >{ centralBody }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate },
                stateAtHalfTimeOfFlight, terminationSettings, propagators::cowell );

    integratorSettings->initialTimeStep_ = - std::fabs( integratorSettings->initialTimeStep_ / 10000.0 );
    integratorSettings->initialTime_ = timeOfFlight / 2.0;

    Eigen::Vector6d currentState = stateAtHalfTimeOfFlight;

    // Backward propagation (corresponding to forward propagation in Sims-Flanagan method).
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        integratorSettings->initialTime_ = timeOfFlight / 2.0 - i * segmentDurationForwardPropagation;
        terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >(
                    ( timeOfFlight / 2.0 ) - ( i + 1 ) * segmentDurationForwardPropagation, true );
        propagatorSettingsImpulsiveDeltaV->resetTerminationSettings( terminationSettings );
        propagatorSettingsImpulsiveDeltaV->resetInitialStates( currentState );

        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettingsImpulsiveDeltaV );

        std::map< double, Eigen::VectorXd > impulsiveDeltaVsResults = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
        currentState = impulsiveDeltaVsResults.begin( )->second;

        std::cout << "impulsive deltaV segment " << std::to_string( numberSegmentsForwardPropagation - 1 - i ) << " : " <<
                     impulsiveDeltaVsResults[ ( timeOfFlight / 2.0 ) - ( i + 1 ) * segmentDurationForwardPropagation ].transpose() << "\n\n"; //.rbegin()->second.transpose() << "\n\n";
        std::cout << "Sims-Flanagan " << std::to_string( numberSegmentsForwardPropagation - 1 - i ) << " : " <<
                     simsFlanaganResults[ ( timeOfFlight / 2.0 ) - ( i + 1 ) * segmentDurationForwardPropagation ].transpose() << "\n\n";
    }


    // FORWARD PROPAGATION.

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( massAtHalfTimeOfFlight );

    // Create termination conditions settings.
    terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight, true );

    //  Create propagator settings.
    propagatorSettingsImpulsiveDeltaV = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                std::vector< std::string >{ centralBody }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate },
                stateAtHalfTimeOfFlightBackwardPropagation, terminationSettings, propagators::cowell );

    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );
    integratorSettings->initialTime_ = timeOfFlight / 2.0;

    currentState = stateAtHalfTimeOfFlightBackwardPropagation;

    // Forward propagation (corresponding to backward propagation in Sims-Flanagan method).
    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        integratorSettings->initialTime_ = ( timeOfFlight/ 2.0 ) +  i * segmentDurationBackwardPropagation;
        terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >(
                    ( timeOfFlight / 2.0 ) + ( i + 1 ) * segmentDurationBackwardPropagation, true );
        propagatorSettingsImpulsiveDeltaV->resetTerminationSettings( terminationSettings );
        propagatorSettingsImpulsiveDeltaV->resetInitialStates( currentState );

        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettingsImpulsiveDeltaV );

        currentState = dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin()->second;

        std::cout << "impulsive deltaV segment " << std::to_string( i + numberSegmentsForwardPropagation ) << " : " <<
                     currentState.transpose() << "\n\n";
        std::cout << "Sims-Flanagan " << std::to_string( i + numberSegmentsForwardPropagation ) << " : " <<
                     simsFlanaganResults[ ( timeOfFlight / 2.0 ) + ( i + 1 ) * segmentDurationBackwardPropagation ].transpose() << "\n\n";
    }

}

//! Test full propagation for Sims Flanagan (assuming constant thrust over each segment).
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_full_propagation )
{
    using namespace low_thrust_direct_methods;

    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.8;
    double specificImpulse = 3000.0;
    double mass = 1800.0;
    int numberSegments = 10;

    // Define (constant) specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return specificImpulse;
    };

    double julianDate = 1000.0 * physical_constants::JULIAN_DAY; //2458849.5;
    double timeOfFlight = 700.0 * physical_constants::JULIAN_DAY;

    std::string bodyToPropagate = "Vehicle";
    std::string centralBody = "Sun";

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ centralBody ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ centralBody ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ centralBody ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ bodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( bodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );


    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Define state at departure and arrival.
    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    // Define integrator settings.
    double stepSize = ( timeOfFlight ) / 500.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize / 10.0 );



    // Define optimisation algorithm.
    algorithm optimisationAlgorithm{ pagmo::de1220() };

    SimsFlanagan simsFlanagan = SimsFlanagan( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction, numberSegments,
                                              timeOfFlight, bodyMap, bodyToPropagate, centralBody, optimisationAlgorithm, 1, 10 );

    std::vector< double > fitnessVector = simsFlanagan.getBestIndividualFitness( );
    std::vector< double > bestIndividual = simsFlanagan.getBestIndividual( );

    std::vector< Eigen::Vector3d > bestThrottles;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        bestThrottles.push_back( ( Eigen::Vector3d( ) << bestIndividual[ i * 3 ], bestIndividual[ i * 3 + 1 ],
                bestIndividual[ i * 3 + 2 ] ).finished( ) );
    }



    ///! TEST SIMS FLANAGAN FULL PROPAGATION

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, std::vector< std::string >{ bodyToPropagate }, std::vector< std::string >{ centralBody } );

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight, true );


    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Define pair of propagatorSettings for backward and forward propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings;

    propagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate }, stateAtDeparture,
                          terminationConditions.first, propagators::cowell, dependentVariablesToSave );

    propagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate }, stateAtArrival,
                          terminationConditions.second, propagators::cowell, dependentVariablesToSave );

    // Define empty maps to store the propagation results.
    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::Vector6d > simsFlanaganResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    // Compute full propagation.
    simsFlanagan.computeSimsFlanaganTrajectoryAndFullPropagation( integratorSettings, propagatorSettings, fullPropagationResults,
                                                                  simsFlanaganResults, dependentVariablesHistory );


    std::cout << "state at departure Sims Flanagan: " << simsFlanaganResults.begin()->second << "\n\n";
    std::cout << "state at departure full propagation: " << fullPropagationResults.begin()->second << "\n\n";
    std::cout << "state at arrival Sims Flanagan: " << simsFlanaganResults.rbegin()->second << "\n\n";
    std::cout << "state at arrival full propagation: " << fullPropagationResults.rbegin()->second << "\n\n";

    input_output::writeDataMapToTextFile( simsFlanaganResults,
                                          "simsFlanaganHFResults.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( fullPropagationResults,
                                          "fullPropagationSFResults.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( dependentVariablesHistory,
                                          "dependentVariablesHistorySimsFlanagan.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    std::cout << "time of flight: " << timeOfFlight << "\n\n";


}


//! Test Sims-Flanagan implementation by comparing it with trajectory obtained with successive impulsive shots applied at times corresponding to
//! half of each of the Sims-Flanagan segments.
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_Shape_Based )
{
    using namespace low_thrust_direct_methods;
    using namespace shape_based_methods;


    /// Shape-based Earth-Mars transfer.

    double julianDate = 9264.5 * physical_constants::JULIAN_DAY;

    double timeOfFlight = 1000.0 * physical_constants::JULIAN_DAY;

    int numberOfRevolutions = 2;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Initialize free coefficients vector for radial velocity function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );


    // Retrieve cartesian state at departure and arrival.
    Eigen::Vector6d cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;

    double scaleFactor = 1.0 / timeOfFlight;

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );


    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    shape_based_methods::HodographicShaping VelocityShapingMethod(
                cartesianStateDepartureBody, cartesianStateArrivalBody,
                timeOfFlight, numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );


    int numberOfSteps = 10000;
    double stepSize = timeOfFlight / static_cast< double >( numberOfSteps );
    double peakAcceleration = 0.0;

    for ( int currentStep = 0 ; currentStep <= numberOfSteps ; currentStep++ ){

        double currentTime = currentStep * stepSize;

        double currentAccelerationMagnitude = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();

        if ( currentAccelerationMagnitude > peakAcceleration )
        {
            peakAcceleration = currentAccelerationMagnitude;
        }

    }


    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 5.0; //0.8;
    double specificImpulse = 3000.0;
    double mass = 2800.0;
    int numberSegments = 50;

    // Calculate number of segments for both the forward propagation (from departure to match point)
    // and the backward propagation (from arrival to match point).
    int numberSegmentsForwardPropagation = ( numberSegments + 1 ) / 2;
    int numberSegmentsBackwardPropagation = numberSegments / 2;
    int segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    int segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );

    std::vector< double > timesAtNodes;
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++)
    {
        timesAtNodes.push_back( i * segmentDurationForwardPropagation );
    }
    for ( int i = 0 ; i <= numberSegmentsBackwardPropagation ; i++ )
    {
        timesAtNodes.push_back( timeOfFlight / 2.0 + i * segmentDurationBackwardPropagation );
    }


    std::map< double, Eigen::Vector6d > thrustMap;
    std::map< double, Eigen::Vector6d > stateMap;
    int currentSegment = 0;
    double currentMass = mass;
    for ( int currentStep = 0 ; currentStep <= numberOfSteps ; currentStep++ ){

        double currentTime = currentStep * stepSize;
        if ( currentTime >= timesAtNodes[ currentSegment + 1 ] )
        {
            currentSegment++;
        }

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        double currentThrust = currentAcceleration * currentMass;
        currentMass *= std::exp( - currentThrust * stepSize /
                                 ( currentMass * specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );

        thrustMap[ currentTime ] = ( Eigen::Vector6d( ) << currentAcceleration, currentThrust, currentMass, 0.0, 0.0, 0.0 ).finished( );
        stateMap[ currentTime ] = VelocityShapingMethod.computeCurrentStateVector( currentTime );
    }

    std::vector< double > thrustMagnitudesPerSegment;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        thrustMagnitudesPerSegment.push_back( ( thrustMap[ timesAtNodes[ i ] ][ 1 ] + thrustMap[ timesAtNodes[ i + 1 ] ][ 1 ] ) / 2.0 );
    }


    thrustMap.clear( );
    currentSegment = 0;
    currentMass = mass;
    double currentMassApprox = mass;
    for ( int currentStep = 0 ; currentStep <= numberOfSteps ; currentStep++ ){

        double currentTime = currentStep * stepSize;
        if ( currentTime > timesAtNodes[ currentSegment + 1 ] )
        {
            currentSegment++;
        }

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        double currentThrust = currentAcceleration * currentMass;
        currentMass *= std::exp( - currentAcceleration * stepSize /
                                 ( specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );

        currentMassApprox += - thrustMagnitudesPerSegment[ currentSegment ] / ( specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION )
                * stepSize;
        double currentAccelerationApprox = thrustMagnitudesPerSegment[ currentSegment ] / currentMassApprox;

        thrustMap[ currentTime ] = ( Eigen::Vector6d( ) << currentAcceleration, currentThrust, currentMass,
                                     currentAccelerationApprox, thrustMagnitudesPerSegment[ currentSegment ], currentMassApprox ).finished( );
    }



    input_output::writeDataMapToTextFile( thrustMap,
                                          "thrustMapSimsFlanagan.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10 );

    input_output::writeDataMapToTextFile( stateMap,
                                          "stateMapSimsFlanagan.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10 );

    std::function< Eigen::Vector3d( const double ) > functionInitialGuess = [ = ]( const double currentTime )
    {
        Eigen::Vector3d currentThrustVector;

        std::function< Eigen::Vector3d( const double ) > computeCurrentThrustNormTest = std::bind(
                    &HodographicShaping::computeCurrentThrustAccelerationVector, VelocityShapingMethod, std::placeholders::_1 );

        int indexSegment;

        if ( currentTime <= timeOfFlight / 2.0 )
        {
            indexSegment = currentTime / segmentDurationForwardPropagation;
        }
        else if ( currentTime == timeOfFlight )
        {
            indexSegment = numberSegments - 1;
        }
        else
        {
            indexSegment = numberSegmentsForwardPropagation + ( currentTime - ( timeOfFlight / 2.0 ) ) / segmentDurationBackwardPropagation;
        }

        currentThrustVector = thrustMagnitudesPerSegment[ indexSegment ] * computeCurrentThrustNormTest( currentTime ).normalized( );


        return currentThrustVector;
    };

    // Define (constant) specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return specificImpulse;
    };

    std::string bodyToPropagate = "Vehicle";
    std::string centralBody = "Sun";


    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ centralBody ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ centralBody ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ centralBody ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ bodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( bodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );


//    // Ephemeris departure body.
//    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

//    // Ephemeris arrival body.
//    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Define state at departure and arrival.
    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    // Define integrator settings.
//    double stepSize = ( timeOfFlight ) / static_cast< double >( 50 );
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize / 10.0 );



    // Define optimisation algorithm.
    algorithm optimisationAlgorithm{ pagmo::de1220() };

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    SimsFlanagan simsFlanagan = SimsFlanagan( cartesianStateDepartureBody, cartesianStateArrivalBody, maximumThrust, specificImpulseFunction, numberSegments,
                                              timeOfFlight, bodyMap, bodyToPropagate, centralBody, optimisationAlgorithm, 20, 1000,
                                              1.0e-6, std::make_pair( functionInitialGuess, 0.1 ) );

    std::vector< double > fitnessVector = simsFlanagan.getBestIndividualFitness( );
    std::vector< double > bestIndividual = simsFlanagan.getBestIndividual( );

    std::vector< Eigen::Vector3d > bestThrottles;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        bestThrottles.push_back( ( Eigen::Vector3d( ) << bestIndividual[ i * 3 ], bestIndividual[ i * 3 + 1 ],
                bestIndividual[ i * 3 + 2 ] ).finished( ) );
    }

    ///! TEST SIMS FLANAGAN FULL PROPAGATION

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, std::vector< std::string >{ bodyToPropagate },
                std::vector< std::string >{ centralBody } );

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight, true );


    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Define pair of propagatorSettings for backward and forward propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings;

    propagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate }, stateAtDeparture,
                          terminationConditions.first, propagators::cowell, dependentVariablesToSave );

    propagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate }, stateAtArrival,
                          terminationConditions.second, propagators::cowell, dependentVariablesToSave );

    // Define empty maps to store the propagation results.
    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::Vector6d > simsFlanaganResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    // Compute full propagation.
    simsFlanagan.computeSimsFlanaganTrajectoryAndFullPropagation( integratorSettings, propagatorSettings, fullPropagationResults,
                                                                  simsFlanaganResults, dependentVariablesHistory );

    input_output::writeDataMapToTextFile( simsFlanaganResults,
                                          "simsFlanaganHFResults.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( fullPropagationResults,
                                          "fullPropagationSFResults.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( dependentVariablesHistory,
                                          "dependentVariablesHistorySimsFlanagan.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    std::cout << "DELTAV SIMS FLANAGAN: " << simsFlanagan.computeDeltaV( ) << "\n\n";
    std::cout << "DELTAV SHAPE BASED: " << VelocityShapingMethod.computeDeltaV( ) << "\n\n";


//    //! Test full propagation w.r.t. impulsive shots as thrust acceleration.

//    // Calculate number of segments for both the forward propagation (from departure to match point)
//    // and the backward propagation (from arrival to match point).
//    int numberSegmentsForwardPropagation = ( numberSegments + 1 ) / 2;
//    int numberSegmentsBackwardPropagation = numberSegments / 2;
//    int segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
//    int segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );

//    // Compute times at half of each segment.
//    std::vector< double > thrustMidTimes;
//    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
//    {
//        thrustMidTimes.push_back( segmentDurationForwardPropagation / 2.0 + i * segmentDurationForwardPropagation );
//    }
//    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
//    {
//        thrustMidTimes.push_back( segmentDurationBackwardPropagation / 2.0 + timeOfFlight / 2.0 + i * segmentDurationBackwardPropagation );
//    }


//    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );
//    SimsFlanaganLeg simsFlanaganLegTest = SimsFlanaganLeg( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction,
//                                                           timeOfFlight, bodyMap, bestThrottles, bodyToPropagate, centralBody );

//    // Compute state at half of the time of flight.
//    simsFlanaganLegTest.propagateForwardFromDepartureToMatchPoint( );
//    simsFlanaganLegTest.propagateBackwardFromArrivalToMatchPoint( );
//    Eigen::Vector6d stateAtHalfTimeOfFlight = simsFlanaganLegTest.getStateAtMatchPointForwardPropagation( );
//    Eigen::Vector6d stateAtHalfTimeOfFlightBackwardPropagation = simsFlanaganLegTest.getStateAtMatchPointBackwardPropagation( );
//    std::cout << "state half TOF test unit test: " << stateAtHalfTimeOfFlight.transpose() << "\n\n";
//    std::cout << "state half TOF test unit test: " << stateAtHalfTimeOfFlightBackwardPropagation.transpose() << "\n\n";

//    // Compute mass at half of the time of flight.
//    double massAtHalfTimeOfFlight = simsFlanaganLegTest.getMassAtMatchPointForwardPropagation( );

//    // Compute deltaVs.
//    double totalManeuverTime = 90.0;
//    double maneuverRiseTime = 15.0;
//    double currentMass = mass;
//    std::vector< Eigen::Vector3d > deltaVs;

//    // Compute deltaVs for the forward propagation half.
//    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
//    {
//        Eigen::Vector3d currentDeltaVvector = maximumThrust * bestThrottles[ i ] * segmentDurationForwardPropagation / currentMass;
//        deltaVs.push_back( currentDeltaVvector );

//        // Update mass.
//        currentMass *= std::exp( - currentDeltaVvector.norm() /
//                                 ( specificImpulseFunction( 0.0 ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );
//    }
//    // Compute deltaVs for the backward propagation half.
//    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
//    {
//        Eigen::Vector3d currentDeltaVvector = maximumThrust * bestThrottles[ i + numberSegmentsForwardPropagation ] *
//                segmentDurationBackwardPropagation / currentMass;
//        deltaVs.push_back( currentDeltaVvector );

//        // Update mass.
//        currentMass *= std::exp( - currentDeltaVvector.norm() /
//                                 ( specificImpulseFunction( 0.0 ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );
//    }
//    std::cout << "current mass: " << currentMass << "\n\n";


//    bodyToPropagateAccelerations.clear();
//    bodyToPropagateAccelerations[ centralBody ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
//                                                                basic_astrodynamics::central_gravity ) );
//    bodyToPropagateAccelerations[ bodyToPropagate ].push_back( std::make_shared< simulation_setup::MomentumWheelDesaturationAccelerationSettings >(
//                                                                   thrustMidTimes, deltaVs, totalManeuverTime, maneuverRiseTime ) );

//    accelerationMap.clear();
//    accelerationMap[ bodyToPropagate ] = bodyToPropagateAccelerations;

//    // Create the acceleration map.
//    accelerationModelMap.clear();
//    accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, std::vector< std::string >{ bodyToPropagate },
//                std::vector< std::string >{ centralBody } );

//    // BACKWARD PROPAGATION.

//    bodyMap[ bodyToPropagate ]->setConstantBodyMass( massAtHalfTimeOfFlight );

//    // Create termination conditions settings.
//    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings
//            = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );

//    //  Create propagator settings.
//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsImpulsiveDeltaV =
//            std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
//                std::vector< std::string >{ centralBody }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate },
//                stateAtHalfTimeOfFlight, terminationSettings, propagators::cowell );

//    integratorSettings->initialTimeStep_ = - std::fabs( integratorSettings->initialTimeStep_ /* / 10000.0 */ );
//    integratorSettings->initialTime_ = timeOfFlight / 2.0;

//    Eigen::Vector6d currentState = stateAtHalfTimeOfFlight;

//    // Backward propagation (corresponding to forward propagation in Sims-Flanagan method).
//    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
//    {
//        integratorSettings->initialTime_ = timeOfFlight / 2.0 - i * segmentDurationForwardPropagation;
//        terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >(
//                    ( timeOfFlight / 2.0 ) - ( i + 1 ) * segmentDurationForwardPropagation, true );
//        propagatorSettingsImpulsiveDeltaV->resetTerminationSettings( terminationSettings );
//        propagatorSettingsImpulsiveDeltaV->resetInitialStates( currentState );

//        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettingsImpulsiveDeltaV );

//        std::map< double, Eigen::VectorXd > impulsiveDeltaVsResults = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
//        currentState = impulsiveDeltaVsResults.begin( )->second;

//        std::cout << "impulsive deltaV segment " << std::to_string( numberSegmentsForwardPropagation - 1 - i ) << " : " <<
//                     impulsiveDeltaVsResults[ ( timeOfFlight / 2.0 ) - ( i + 1 ) * segmentDurationForwardPropagation ].transpose() << "\n\n"; //.rbegin()->second.transpose() << "\n\n";
//        std::cout << "Sims-Flanagan " << std::to_string( numberSegmentsForwardPropagation - 1 - i ) << " : " <<
//                     simsFlanaganResults[ ( timeOfFlight / 2.0 ) - ( i + 1 ) * segmentDurationForwardPropagation ].transpose() << "\n\n";
//    }


//    // FORWARD PROPAGATION.

//    bodyMap[ bodyToPropagate ]->setConstantBodyMass( massAtHalfTimeOfFlight );

//    // Create termination conditions settings.
//    terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight, true );

//    //  Create propagator settings.
//    propagatorSettingsImpulsiveDeltaV = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
//                std::vector< std::string >{ centralBody }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate },
//                stateAtHalfTimeOfFlightBackwardPropagation, terminationSettings, propagators::cowell );

//    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );
//    integratorSettings->initialTime_ = timeOfFlight / 2.0;

//    currentState = stateAtHalfTimeOfFlightBackwardPropagation;

//    // Forward propagation (corresponding to backward propagation in Sims-Flanagan method).
//    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
//    {
//        integratorSettings->initialTime_ = ( timeOfFlight / 2.0 ) +  i * segmentDurationBackwardPropagation;
//        terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >(
//                    ( timeOfFlight / 2.0 ) + ( i + 1 ) * segmentDurationBackwardPropagation, true );
//        propagatorSettingsImpulsiveDeltaV->resetTerminationSettings( terminationSettings );
//        propagatorSettingsImpulsiveDeltaV->resetInitialStates( currentState );

//        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettingsImpulsiveDeltaV );

//        currentState = dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin()->second;

//        std::cout << "impulsive deltaV segment " << std::to_string( i + numberSegmentsForwardPropagation ) << " : " <<
//                     currentState.transpose() << "\n\n";
//        std::cout << "Sims-Flanagan " << std::to_string( i + numberSegmentsForwardPropagation ) << " : " <<
//                     simsFlanaganResults[ ( timeOfFlight / 2.0 ) + ( i + 1 ) * segmentDurationBackwardPropagation ].transpose() << "\n\n";
//    }



}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
