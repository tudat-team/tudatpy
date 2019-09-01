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
#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethodLeg.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethod.h"
#include "pagmo/algorithms/de1220.hpp"

namespace tudat
{
namespace unit_tests
{

//! Test hybrid method implementation.
BOOST_AUTO_TEST_SUITE( test_hybrid_method )

BOOST_AUTO_TEST_CASE( test_hybrid_method_implementation )
{
    using namespace low_thrust_direct_methods;

    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.450;
    double specificImpulse = 3000.0;
    double mass = 1800.0;

    // Define (constant) specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return specificImpulse;
    };

    double julianDate = 1000.0 * physical_constants::JULIAN_DAY; //2458849.5;
    double timeOfFlight = 100.0 * physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate, julianDate - 300.0, julianDate + timeOfFlight + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    std::string bodyToPropagate = "Vehicle";
    std::string centralBody = "Earth";

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

    // Initial and final elements.
    Eigen::Vector6d initialKeplerianElements = ( Eigen::Vector6d( ) << 24505.9e3, 0.725, 7.0 * mathematical_constants::PI / 180.0,
                                                 0.0, 0.0, 0.0 ).finished( );
    Eigen::Vector6d finalKeplerianElements = ( Eigen::Vector6d( ) << 42164.65e3, 5.53e-4, 7.41e-5 * mathematical_constants::PI / 180.0,
                                               0.0, 0.0, 0.0 ).finished( );

    stateAtDeparture = orbital_element_conversions::convertKeplerianToCartesianElements(
                initialKeplerianElements, bodyMap[ "Earth" ]->getGravityFieldModel()->getGravitationalParameter() );
    stateAtArrival = orbital_element_conversions::convertKeplerianToCartesianElements(
                finalKeplerianElements, bodyMap[ "Earth" ]->getGravityFieldModel()->getGravitationalParameter() );

    // Define integrator settings.
    double stepSize = ( timeOfFlight ) / static_cast< double >( 50 );
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > >
            ( numerical_integrators::rungeKutta4, 0.0, stepSize / 1000.0 );


    Eigen::VectorXd initialCostates; initialCostates.resize( 5 );
    Eigen::VectorXd finalCostates; finalCostates.resize( 5 );

    for ( int i = 0 ; i < 5 ; i++ )
    {
        initialCostates[ i ] = 0.0;
        finalCostates[ i ] = 1.0;
    }

    HybridMethodLeg hybridMethodLeg = HybridMethodLeg( stateAtDeparture, stateAtArrival, initialCostates, finalCostates, /*0,*/
                                                       maximumThrust, specificImpulseFunction, timeOfFlight, bodyMap,
                                                       bodyToPropagate, centralBody );

    std::cout.precision( 20 );
    std::cout << "deltaV: " << hybridMethodLeg.getTotalDeltaV( ) << "\n\n";

    Eigen::Vector6d finalState = hybridMethodLeg.propagateTrajectory( integratorSettings );

    std::cout << "final state after propagation: " << finalState << "\n\n";

    Eigen::Vector6d finalStateTest = hybridMethodLeg.propagateTrajectory( 0.0, timeOfFlight, stateAtDeparture, mass, integratorSettings );
    std::cout << "final state after propagation test: " << finalStateTest << "\n\n";

    std::cout << "deltaV after propagation: " << hybridMethodLeg.getTotalDeltaV() << "\n\n";
    std::cout << "confirmation computation deltaV: " << hybridMethodLeg.computeTotalDeltaV() << "\n\n";
//    std::cout << "mass at time of flight: " << hybridMethodLeg.getMassAtTimeOfFlight() << "\n\n";
//    std::cout << "initial mass: " << mass << "\n\n";
//    std::cout << "delta m: " << mass - hybridMethodLeg.getMassAtTimeOfFlight() << "\n\n";
//    std::cout << "time of flight: " << timeOfFlight << "\n\n";

//    std::cout << "state at departure: " << stateAtDeparture << "\n\n";
//    std::cout << "state at arrival: " << stateAtArrival << "\n\n";


    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    HybridMethodProblem problem = HybridMethodProblem( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction,
                                                       /*1.0,*/ timeOfFlight, bodyMap, bodyToPropagate, centralBody, integratorSettings, 1.0e-6 );

//    std::vector< double > designVariables;
//    for ( int i = 0 ; i < 5 ; i++ )
//    {
//        designVariables.push_back( 0.0 );
//        designVariables.push_back( 1.0 );
//    }
//    std::vector< double > output = problem.fitness( designVariables );

//    std::cout << "size output vector: " << output.size() << "\n\n";
//    for ( int i = 0 ; i < output.size() ; i++ )
//    {
//        std::cout << "output: " << output[ i ] << "\n\n";
//    }


    // Define optimisation algorithm.
    algorithm optimisationAlgorithm{ pagmo::de1220() };

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    HybridMethod hybridMethod = HybridMethod( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction, /*1.0,*/
                                              timeOfFlight, bodyMap, bodyToPropagate, centralBody, optimisationAlgorithm,
                                              integratorSettings, 1.0e-3 );

    std::pair< std::vector< double >, std::vector< double > > champion = hybridMethod.performOptimisation();

    std::vector< double > fitnessVector = champion.first;
    std::vector< double > bestIndividual = champion.second;

    std::vector< double > bestDesignVariables = bestIndividual;
    std::vector< double > bestOutput = problem.fitness( bestDesignVariables );

    std::cout << "size output vector: " << bestOutput.size() << "\n\n";
    for ( int i = 0 ; i < bestOutput.size() ; i++ )
    {
        std::cout << "output: " << bestOutput[ i ] << "\n\n";
    }


    // Test full propagation.

    // Define pair of propagatorSettings for backward and forward propagation.

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Earth" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
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

    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings;

    propagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate }, stateAtDeparture,
                          terminationConditions.first, propagators::cowell, std::shared_ptr< propagators::DependentVariableSaveSettings >( ) );

    propagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate }, stateAtArrival,
                          terminationConditions.second, propagators::cowell, std::shared_ptr< propagators::DependentVariableSaveSettings >( ) );

    // Define empty maps to store the propagation results.
    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::Vector6d > hybridMethodResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    // Compute full propagation.
    hybridMethod.computeHybridMethodTrajectoryAndFullPropagation( propagatorSettings, fullPropagationResults,
                                                                  hybridMethodResults, dependentVariablesHistory );





    Eigen::VectorXd bestInitialMEEcostates; bestInitialMEEcostates.resize( 5 );
    Eigen::VectorXd bestFinalMEEcostates; bestFinalMEEcostates.resize( 5 );
    for ( int i = 0 ; i < 5 ; i++ )
    {
        bestInitialMEEcostates[ i ] = hybridMethod.getBestIndividual( )[ i ];
        bestFinalMEEcostates[ i ] = hybridMethod.getBestIndividual( )[ i + 5 ]; // championDesignVariables_[ i + 5 ];
    }

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    HybridMethodLeg hybridMethodLegTest = HybridMethodLeg( stateAtDeparture, stateAtArrival, bestInitialMEEcostates, bestFinalMEEcostates, /*0,*/
                                                       maximumThrust, specificImpulseFunction, timeOfFlight, bodyMap,
                                                       bodyToPropagate, centralBody );

    std::cout.precision( 20 );


    Eigen::Vector6d finalStateTestFullPropagation = hybridMethodLegTest.propagateTrajectory( 0.0, timeOfFlight,
                                                                                             stateAtDeparture, mass, integratorSettings );
    std::cout << "final state after propagation test: " << finalStateTestFullPropagation << "\n\n";

    std::cout << "state at departure hybrid method: " << hybridMethodResults.begin()->second << "\n\n";
    std::cout << "state at departure full propagation: " << fullPropagationResults.begin()->second << "\n\n";
    std::cout << "state at arrival hybrid method: " << hybridMethodResults.rbegin()->second << "\n\n";
    std::cout << "state at arrival full propagation: " << fullPropagationResults.rbegin()->second << "\n\n";

    std::cout << "state at departure: " << stateAtDeparture << "\n\n";


    input_output::writeDataMapToTextFile( hybridMethodResults,
                                          "hybridMethodResults.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( fullPropagationResults,
                                          "fullPropagationHybridMethodResults.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
