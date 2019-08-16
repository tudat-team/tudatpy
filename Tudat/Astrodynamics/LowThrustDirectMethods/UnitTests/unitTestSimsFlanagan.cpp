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
#include "pagmo/algorithms/de1220.hpp"

namespace tudat
{
namespace unit_tests
{

//! Test Sims Flanagan implementation.
BOOST_AUTO_TEST_SUITE( test_Sims_Flanagan )

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

//    std::vector< std::string > bodiesToPropagate;
//    bodiesToPropagate.push_back( bodyToPropagate );
//    std::vector< std::string > centralBodies;
//    centralBodies.push_back( centralBody );


    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );

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



//    SimsFlanaganProblem problem = SimsFlanaganProblem( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction,
//                                                       numberSegments, timeOfFlight, bodyMap, bodyToPropagate, centralBody, integratorSettings,
//                                                       propagators::cowell, false );

//    std::vector< double > designVariables;
//    for ( int i = 0 ; i < numberSegments ; i++ )
//    {
//        designVariables.push_back( 1.0 );
//        designVariables.push_back( 0.0 );
//        designVariables.push_back( 0.0 );
//    }
//    std::vector< double > output = problem.fitness( designVariables );

//    std::cout << "size output vector: " << output.size() << "\n\n";
//    for ( int i = 0 ; i < output.size() ; i++ )
//    {
//        std::cout << "output: " << output[ i ] << "\n\n";
//    }


    // Define optimisation algorithm.
    algorithm optimisationAlgorithm{ pagmo::de1220() };

    SimsFlanagan simsFlanagan = SimsFlanagan( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction, numberSegments,
                                              timeOfFlight, bodyMap, bodyToPropagate, centralBody, optimisationAlgorithm,
                                              integratorSettings, propagators::cowell, true );

    std::pair< std::vector< double >, std::vector< double > > champion = simsFlanagan.performOptimisation();

    std::vector< double > fitnessVector = champion.first;
    std::vector< double > bestIndividual = champion.second;

    std::vector< Eigen::Vector3d > bestThrottles;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        bestThrottles.push_back( ( Eigen::Vector3d( ) << bestIndividual[ i * 3 ], bestIndividual[ i * 3 + 1 ],
                bestIndividual[ i * 3 + 2 ] ).finished( ) );
    }

    /// TEST PROPAGATE SIMS FLANAGAN SOLUTION TO GIVEN TIME.
//    // Re-initialise mass of the spacecraft in the body map.
//    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

//    SimsFlanaganLeg simsFlanaganLeg = SimsFlanaganLeg( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulseFunction,
//                                                       timeOfFlight, bodyMap, bestThrottles, bodyToPropagate, centralBody );

//    std::cout << "propagate low order solution: " << simsFlanaganLeg.propagateTrajectory( timeOfFlight / 2.0 ) << "\n\n";
//    std::cout << "propagate high order solution: " << simsFlanaganLeg.propagateTrajectoryHighOrderSolution(
//                     timeOfFlight / 2.0, integratorSettings, propagators::cowell ) << "\n\n";


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
    simsFlanagan.computeSimsFlanaganTrajectoryAndFullPropagation( propagatorSettings, fullPropagationResults,
                                                                  simsFlanaganResults, dependentVariablesHistory );


    std::cout << "state at departure Sims Flanagan: " << simsFlanaganResults.begin()->second << "\n\n";
    std::cout << "state at departure full propagation: " << fullPropagationResults.begin()->second << "\n\n";
    std::cout << "state at arrival Sims Flanagan: " << simsFlanaganResults.rbegin()->second << "\n\n";
    std::cout << "state at arrival full propagation: " << fullPropagationResults.rbegin()->second << "\n\n";

    std::cout << "time of flight: " << timeOfFlight << "\n\n";


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
