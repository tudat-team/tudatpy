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
#include "Tudat/Astrodynamics/LowThrustDirectMethods/lowThrustLeg.h"

namespace tudat
{
namespace unit_tests
{

//! Test hybrid method implementation.
BOOST_AUTO_TEST_SUITE( test_hybrid_method )

BOOST_AUTO_TEST_CASE( test_hybrid_method_implementation )
{
    using namespace low_thrust_direct_methods;
    using namespace transfer_trajectories;

    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.450;
    double specificImpulse = 3000.0;
    double mass = 1800.0;

//    // Define (constant) specific impulse function.
//    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
//    {
//        return specificImpulse;
//    };

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


//    // Ephemeris departure body.
//    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

//    // Ephemeris arrival body.
//    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

//    // Define state at departure and arrival.
//    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
//    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    // Initial and final states in keplerian elements.
    Eigen::Vector6d initialKeplerianElements = ( Eigen::Vector6d( ) << 24505.9e3, 0.725, 7.0 * mathematical_constants::PI / 180.0,
                                                 0.0, 0.0, 0.0 ).finished( );
    Eigen::Vector6d finalKeplerianElements = ( Eigen::Vector6d( ) << 42164.65e3, 5.53e-4, 7.41e-5 * mathematical_constants::PI / 180.0,
                                               0.0, 0.0, 0.0 ).finished( );

    // Initial and final states in cartesian coordinates.
    Eigen::Vector6d stateAtDeparture = orbital_element_conversions::convertKeplerianToCartesianElements(
                initialKeplerianElements, bodyMap[ "Earth" ]->getGravityFieldModel()->getGravitationalParameter() );
    Eigen::Vector6d stateAtArrival = orbital_element_conversions::convertKeplerianToCartesianElements(
                finalKeplerianElements, bodyMap[ "Earth" ]->getGravityFieldModel()->getGravitationalParameter() );

    // Define integrator settings.
    double stepSize = ( timeOfFlight ) / static_cast< double >( 200 );
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > >
            ( numerical_integrators::rungeKutta4, 0.0, stepSize / 200.0 );



    Eigen::VectorXd initialCostates; initialCostates.resize( 5 );
    Eigen::VectorXd finalCostates; finalCostates.resize( 5 );

    for ( int i = 0 ; i < 5 ; i++ )
    {
        initialCostates[ i ] = 0.0;
        finalCostates[ i ] = 1.0;
    }

    HybridMethodLeg hybridMethodLeg = HybridMethodLeg( stateAtDeparture, stateAtArrival, initialCostates, finalCostates,
                                                       maximumThrust, specificImpulse, timeOfFlight, bodyMap,
                                                       bodyToPropagate, centralBody, integratorSettings );

    std::cout.precision( 20 );
    std::cout << "deltaV: " << hybridMethodLeg.getTotalDeltaV( ) << "\n\n";

    Eigen::Vector6d finalState = hybridMethodLeg.propagateTrajectory( );

    std::cout << "final state after propagation: " << finalState << "\n\n";

    Eigen::Vector6d finalStateTest = hybridMethodLeg.propagateTrajectory( 0.0, timeOfFlight, stateAtDeparture, mass/*, integratorSettings*/ );
    std::cout << "final state after propagation test: " << finalStateTest << "\n\n";

    std::cout << "deltaV after propagation: " << hybridMethodLeg.getTotalDeltaV() << "\n\n";
    std::cout << "confirmation computation deltaV: " << hybridMethodLeg.computeDeltaV() << "\n\n";


    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

//    Eigen::Matrix< double, 10, 1 > initialGuess = Eigen::MatrixXd::Zero( 10, 1 );

//    HybridMethodProblem problem = HybridMethodProblem( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulse,
//                                                       timeOfFlight, bodyMap, bodyToPropagate, centralBody, integratorSettings,
//                                                       std::make_pair( initialGuess, TUDAT_NAN ), 1.0e-6 );


    // Define optimisation algorithm.
    algorithm optimisationAlgorithm{ pagmo::de1220() };

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    std::shared_ptr< OptimisationSettings > optimisationSettings = std::make_shared< OptimisationSettings >( optimisationAlgorithm, 1, 10, 1.0e-3 );

    HybridMethod hybridMethod = HybridMethod( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulse,
                                              timeOfFlight, bodyMap, bodyToPropagate, centralBody, integratorSettings,
                                               /*optimisationAlgorithm, 1, 10,*/ optimisationSettings/*, 1.0e-3*/ );

    std::pair< std::vector< double >, std::vector< double > > champion = hybridMethod.performOptimisation();

    std::vector< double > fitnessVector =  hybridMethod.getBestIndividualFitness( );
    std::vector< double > bestIndividual = hybridMethod.getBestIndividual( );

    std::vector< double > bestDesignVariables = bestIndividual;
//    std::vector< double > bestOutput = problem.fitness( bestDesignVariables );

//    std::cout << "size output vector: " << bestOutput.size() << "\n\n";
//    for ( int i = 0 ; i < bestOutput.size() ; i++ )
//    {
//        std::cout << "output: " << bestOutput[ i ] << "\n\n";
//    }


    // Test full propagation.

    // Define pair of propagatorSettings for backward and forward propagation.

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    std::function< double( const double ) > specificImpulseFunction = [ = ] ( const double currentTime )
    {
        return specificImpulse;
    };


    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight, true );

    // Define empty maps to store the propagation results.
    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::Vector6d > hybridMethodResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    basic_astrodynamics::AccelerationMap perturbingAccelerationsMap;

    // Create complete propagation settings (backward and forward propagations).
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
            std::shared_ptr< propagators::PropagatorSettings< double > > > propagatorSettings = hybridMethod.createLowThrustPropagatorSettings(
                 specificImpulseFunction, perturbingAccelerationsMap, integratorSettings, dependentVariablesToSave );

    // Compute full propagation.
    hybridMethod.computeSemiAnalyticalAndFullPropagation( integratorSettings, propagatorSettings, fullPropagationResults,
                                                                  hybridMethodResults, dependentVariablesHistory );





    Eigen::VectorXd bestInitialMEEcostates; bestInitialMEEcostates.resize( 5 );
    Eigen::VectorXd bestFinalMEEcostates; bestFinalMEEcostates.resize( 5 );
    for ( int i = 0 ; i < 5 ; i++ )
    {
        bestInitialMEEcostates[ i ] = hybridMethod.getBestIndividual( )[ i ];
        bestFinalMEEcostates[ i ] = hybridMethod.getBestIndividual( )[ i + 5 ];
    }

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    HybridMethodLeg hybridMethodLegTest = HybridMethodLeg( stateAtDeparture, stateAtArrival, bestInitialMEEcostates, bestFinalMEEcostates,
                                                       maximumThrust, specificImpulse, timeOfFlight, bodyMap,
                                                       bodyToPropagate, centralBody, integratorSettings );

    std::cout.precision( 20 );


    Eigen::Vector6d finalStateTestFullPropagation = hybridMethodLegTest.propagateTrajectory( 0.0, timeOfFlight, stateAtDeparture, mass );
    std::cout << "final state after propagation test: " << finalStateTestFullPropagation << "\n\n";

    std::cout << "state at departure hybrid method: " << hybridMethodResults.begin()->second << "\n\n";
    std::cout << "state at departure full propagation: " << fullPropagationResults.begin()->second << "\n\n";
    std::cout << "state at arrival hybrid method: " << hybridMethodResults.rbegin()->second << "\n\n";
    std::cout << "state at arrival full propagation: " << fullPropagationResults.rbegin()->second << "\n\n";

    for ( int i = 0 ; i < 3 ; i++ )
    {
        BOOST_CHECK_SMALL( ( std::fabs(  hybridMethodResults.begin()->second[ i ] - fullPropagationResults.begin()->second[ i ] ) / stateAtDeparture.segment( 0, 3 ).norm( ) ), 2.0e-6 );
        BOOST_CHECK_SMALL( ( std::fabs(  hybridMethodResults.begin()->second[ i + 3 ] - fullPropagationResults.begin()->second[ i + 3 ] ) / stateAtDeparture.segment( 3, 3 ).norm( ) ), 2.0e-6 );

        BOOST_CHECK_SMALL( ( std::fabs(  hybridMethodResults.rbegin()->second[ i ] - fullPropagationResults.rbegin()->second[ i ] ) / stateAtArrival.segment( 0, 3 ).norm( ) ), 2.0e-6 );
        BOOST_CHECK_SMALL( ( std::fabs(  hybridMethodResults.rbegin()->second[ i + 3 ] - fullPropagationResults.rbegin()->second[ i + 3 ] ) / stateAtArrival.segment( 3, 3 ).norm( ) ), 2.0e-6 );

//        BOOST_CHECK_SMALL( ( std::fabs( trajectory[ itr->first ][ i + 3 ] - currentExpectedState[ i + 3 ] ) / currentExpectedState.segment( 3, 3 ).norm( ) ), 1.0e-6 );
//        BOOST_CHECK_SMALL( ( std::fabs( thrustProfile[ itr->first ][ i ] - currentExpectedThrust[ i ] ) ), 1.0e-8 );
//        BOOST_CHECK_SMALL( ( std::fabs( thrustAccelerationProfile[ itr->first ][ i ] - currentExpectedThrustAcceleration[ i ] ) ), 1.0e-10 );

//        BOOST_CHECK_SMALL( ( std::fabs( calculatedState[ i ] - currentExpectedState[ i ] ) / currentExpectedState.segment( 0, 3 ).norm( ) ), 1.0e-6 );
//        BOOST_CHECK_SMALL( ( std::fabs( calculatedState[ i + 3 ] - currentExpectedState[ i + 3 ] ) / currentExpectedState.segment( 3, 3 ).norm( ) ), 1.0e-6 );
//        BOOST_CHECK_SMALL( ( std::fabs( calculatedThrust[ i ] - currentExpectedThrust[ i ] ) ), 1.0e-8 );
//        BOOST_CHECK_SMALL( ( std::fabs( calculatedThrustAcceleration[ i ] - currentExpectedThrustAcceleration[ i ] ) ), 1.0e-10 );
    }

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


    //! TEST PROFILES

    int numberSteps = 10;
    std::vector< double > epochsVector;
//    for ( std::map< double, Eigen::Vector6d >::iterator itr = hybridMethodResults.begin( ) ; itr != hybridMethodResults.end( ) ; itr++ )
//    {
//        epochsVector.push_back( itr->first );
//    }
    for ( int i = 1 ; i <= numberSteps ; i++ )
    {
        epochsVector.push_back( timeOfFlight / numberSteps * i );
    }
//    epochsVector.push_back( 0.0 );
//    epochsVector.push_back( timeOfFlight / 4.0 );
//    epochsVector.push_back( timeOfFlight / 2.0 );
//    epochsVector.push_back( 3.0 * timeOfFlight / 4.0 );
//    epochsVector.push_back( timeOfFlight );

    std::map< double, Eigen::Vector6d > trajectory;
    std::map< double, Eigen::VectorXd > massProfile;
    std::map< double, Eigen::VectorXd > thrustProfile;
    std::map< double, Eigen::VectorXd > thrustAccelerationProfile;

    hybridMethod.getTrajectory( epochsVector, trajectory );
    hybridMethod.getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );
    hybridMethod.getThrustProfile( epochsVector, thrustProfile, specificImpulseFunction, integratorSettings );
    hybridMethod.getThrustAccelerationProfile( epochsVector, thrustAccelerationProfile, specificImpulseFunction, integratorSettings );


    std::function< Eigen::VectorXd( const double ) > costatesFunction = [ = ]( const double currentTime )
    {
        Eigen::VectorXd currentCostates;
        currentCostates.resize( 5 );

        for ( int i = 0 ; i < 5 ; i++ )
        {
            currentCostates[ i ] = bestInitialMEEcostates[ i ]
                    + ( currentTime / timeOfFlight ) * ( bestFinalMEEcostates[ i ] - bestInitialMEEcostates[ i ] );
        }
        return currentCostates;
    };

    // Define thrust direction settings from the MEE costates.
    std::shared_ptr< simulation_setup::MeeCostateBasedThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::MeeCostateBasedThrustDirectionSettings >( bodyToPropagate, centralBody, costatesFunction );

//    std::function< double ( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
//    {
//      return specificImpulse_;
//    };

    // Define bang-bang thrust magnitude settings based on MEE co-states.
    std::shared_ptr< simulation_setup::FromMeeCostatesBangBangThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromMeeCostatesBangBangThrustMagnitudeSettings >(
                maximumThrust, specificImpulseFunction, costatesFunction, bodyToPropagate, centralBody );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >( thrustDirectionSettings, thrustMagnitudeSettings );





    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ centralBody ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                    basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ bodyToPropagate ].push_back( thrustAccelerationSettings );

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, std::vector< std::string >{ bodyToPropagate }, std::vector< std::string >{ centralBody } );



//    // Retrieve low-thrust trajectory nominal accelerations (thrust + central gravity accelerations).
//    basic_astrodynamics::AccelerationMap lowThrustTrajectoryAccelerations = hybridMethod.retrieveLowThrustAccelerationMap( specificImpulseFunction, integratorSettings );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ bodyToPropagate ] = simulation_setup::createMassRateModel( bodyToPropagate, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap, accelerationModelMap );

    // Define list of dependent variables to save.
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_mass_rate_dependent_variables, bodyToPropagate ) );

    // Create object with list of dependent variables
    dependentVariablesToSave = std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );


    double initialTime = 0.0;
    Eigen::Vector6d currentState = stateAtDeparture;
    double currentMass = mass;
    for ( std::map< double, Eigen::Vector6d >::iterator itr = trajectory.begin( ) ; itr != trajectory.end( ) ; itr++ )
    {
        double currentEpoch = itr->first;

        // Create termination conditions settings.
        std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings =
                std::make_shared< propagators::PropagationTimeTerminationSettings >( currentEpoch, true );

        // Define translational state propagation settings
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
                std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                            ( std::vector< std::string >{ centralBody }, accelerationModelMap,
                              std::vector< std::string >{ bodyToPropagate }, currentState,
                              terminationSettings, propagators::gauss_modified_equinoctial, dependentVariablesToSave );

        // Create settings for propagating the mass of the vehicle.
        std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings =
                std::make_shared< propagators::MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate }, massRateModels, ( Eigen::Matrix< double, 1, 1 >( ) << currentMass ).finished( ),
                    terminationSettings );

        integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );
        integratorSettings->initialTime_ = initialTime;

        // Create list of propagation settings.
        std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > propagatorSettingsVector;

        // Backward propagator settings vector.
        propagatorSettingsVector.push_back( translationalStatePropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        // Define propagator settings.
        std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings = std::make_shared< propagators::MultiTypePropagatorSettings< double > >(
                    propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

        bodyMap[ bodyToPropagate ]->setConstantBodyMass( currentMass );

        // Perform the backward propagation.
        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
        std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

        Eigen::Vector6d currentExpectedState = stateHistory.rbegin( )->second.segment( 0, 6 );
        double currentExpectedMass = stateHistory.rbegin( )->second[ 6 ];
        Eigen::Vector3d currentExpectedThrustAcceleration =  dependentVariableHistory.rbegin( )->second.segment( 0, 3 );
        Eigen::Vector3d currentExpectedThrust = currentExpectedThrustAcceleration * currentExpectedMass;

//        Eigen::Vector6d calculatedState = hybridMethod.computeCurrentStateVector( currentEpoch );
//        double calculatedMass = hybridMethod.computeCurrentMass( currentEpoch, specificImpulseFunction, integratorSettings );
//        Eigen::Vector3d calculatedThrust = hybridMethod.computeCurrentThrust( currentEpoch, specificImpulseFunction, integratorSettings ).transpose( );
//        Eigen::Vector3d calculatedThrustAcceleration = hybridMethod.computeCurrentThrustAcceleration( currentEpoch, specificImpulseFunction, integratorSettings );

        for ( int i = 0 ; i < 3 ; i++ )
        {
            BOOST_CHECK_SMALL( ( std::fabs( trajectory[ itr->first ][ i ] - currentExpectedState[ i ] ) / currentExpectedState.segment( 0, 3 ).norm( ) ), 1.0e-12 );
            BOOST_CHECK_SMALL( ( std::fabs( trajectory[ itr->first ][ i + 3 ] - currentExpectedState[ i + 3 ] ) / currentExpectedState.segment( 3, 3 ).norm( ) ), 1.0e-12 );
            BOOST_CHECK_SMALL( ( std::fabs( thrustProfile[ itr->first ][ i ] - currentExpectedThrust[ i ] ) ), 1.0e-8 );
            BOOST_CHECK_SMALL( ( std::fabs( thrustAccelerationProfile[ itr->first ][ i ] - currentExpectedThrustAcceleration[ i ] ) ), 1.0e-10 );

//            BOOST_CHECK_SMALL( ( std::fabs( calculatedState[ i ] - currentExpectedState[ i ] ) / currentExpectedState.segment( 0, 3 ).norm( ) ), 1.0e-6 );
//            BOOST_CHECK_SMALL( ( std::fabs( calculatedState[ i + 3 ] - currentExpectedState[ i + 3 ] ) / currentExpectedState.segment( 3, 3 ).norm( ) ), 1.0e-6 );
//            BOOST_CHECK_SMALL( ( std::fabs( calculatedThrust[ i ] - currentExpectedThrust[ i ] ) ), 1.0e-7 );
//            BOOST_CHECK_SMALL( ( std::fabs( calculatedThrustAcceleration[ i ] - currentExpectedThrustAcceleration[ i ] ) ), 1.0e-10 );
        }
        BOOST_CHECK_SMALL( std::fabs( massProfile[ itr->first ][ 0 ] - currentExpectedMass ), 1.0e-15 );
//        BOOST_CHECK_SMALL( std::fabs( calculatedMass - currentExpectedMass ), 1.0e-11 );

//        double manuallyCalculatedMass = mass;
//        double previousEpoch;
//        for ( std::map< double, Eigen::VectorXd >::iterator itrDependentVariables = dependentVariableHistory.begin( ) ;
//              itrDependentVariables != dependentVariableHistory.end( ) ; itrDependentVariables++ )
//        {

//            if ( itrDependentVariables == dependentVariableHistory.begin( ) )
//            {
//                previousEpoch = itrDependentVariables->first;
//            }
//            else
//            {

//                double massRate;
//                Eigen::Vector3d currentAcceleration = itrDependentVariables->second.segment( 0, 3 );
//                Eigen::Vector3d currentThrust;
//                if ( currentAcceleration.norm( ) > 0.0 )
//                {
//                    currentThrust = maximumThrust * currentAcceleration.normalized( );
//                    massRate = - maximumThrust / ( specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );
//                }
//                else
//                {
//                    currentThrust = Eigen::Vector3d::Zero( );
//                    massRate = 0.0;
//                }

//                double currentEpoch = itrDependentVariables->first;
//                manuallyCalculatedMass += massRate * ( currentEpoch - previousEpoch  );

//                BOOST_CHECK_SMALL( std::fabs( manuallyCalculatedMass - stateHistory[ itrDependentVariables->first ][ 6 ] ), 1.0e-12 );

//                previousEpoch = currentEpoch;

//            }

//        }

        initialTime = currentEpoch;
        currentMass = currentExpectedMass;
        currentState = currentExpectedState;

    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
