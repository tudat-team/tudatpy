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
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethodModel.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethod.h"
#include "pagmo/algorithms/de1220.hpp"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLeg.h"

namespace tudat
{
namespace unit_tests
{

//! Test hybrid method implementation.
BOOST_AUTO_TEST_SUITE( test_hybrid_method )

BOOST_AUTO_TEST_CASE( test_hybrid_method_implementation )
{
    using namespace low_thrust_trajectories;

    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.450;
    double specificImpulse = 3000.0;
    double mass = 1800.0;

    std::function< double( const double ) > specificImpulseFunction = [ = ] ( const double currentTime )
    {
        return specificImpulse;
    };

    double julianDate = 1000.0 * physical_constants::JULIAN_DAY;
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


    // Create spacecraft object.
    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    std::string bodyToPropagate = "Vehicle";
    std::string centralBody = "Earth";

    // Set vehicle mass.
    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

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
    double stepSize = ( timeOfFlight ) / static_cast< double >( 40000 );
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > >
            ( numerical_integrators::rungeKutta4, 0.0, stepSize );


    // Define optimisation algorithm.
    algorithm optimisationAlgorithm{ pagmo::de1220() };

    std::shared_ptr< simulation_setup::OptimisationSettings > optimisationSettings =
            std::make_shared< simulation_setup::OptimisationSettings >( optimisationAlgorithm, 1, 10, 1.0e-3 );

    // Create hybrid method trajectory.
    HybridMethod hybridMethod = HybridMethod( stateAtDeparture, stateAtArrival, maximumThrust, specificImpulse,
                                              timeOfFlight, bodyMap, bodyToPropagate, centralBody, integratorSettings,
                                              optimisationSettings );


    // Retrieve optimisation output.
    std::vector< double > fitnessVector =  hybridMethod.getBestIndividualFitness( );
    std::vector< double > bestDesignVariables = hybridMethod.getBestIndividual( );


    // Define MEE costates function.
    Eigen::VectorXd bestInitialMEEcostates; bestInitialMEEcostates.resize( 5 );
    Eigen::VectorXd bestFinalMEEcostates; bestFinalMEEcostates.resize( 5 );
    for ( int i = 0 ; i < 5 ; i++ )
    {
        bestInitialMEEcostates[ i ] = hybridMethod.getBestIndividual( )[ i ];
        bestFinalMEEcostates[ i ] = hybridMethod.getBestIndividual( )[ i + 5 ];
    }

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



    // Retrieve trajectory, mass, thrust, and thrust acceleration for a given set of epochs.
    int numberSteps = 10;
    std::vector< double > epochsVector;
    for ( int i = 1 ; i <= numberSteps ; i++ )
    {
        epochsVector.push_back( timeOfFlight / numberSteps * i );
    }

    std::map< double, Eigen::Vector6d > trajectory;
    std::map< double, Eigen::VectorXd > massProfile;
    std::map< double, Eigen::VectorXd > thrustProfile;
    std::map< double, Eigen::VectorXd > thrustAccelerationProfile;

    hybridMethod.getTrajectory( epochsVector, trajectory );
    hybridMethod.getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );
    hybridMethod.getThrustForceProfile( epochsVector, thrustProfile, specificImpulseFunction, integratorSettings );
    hybridMethod.getThrustAccelerationProfile( epochsVector, thrustAccelerationProfile, specificImpulseFunction, integratorSettings );


    /// PROPAGATE THE TRAJECTORY NUMERICALLY (SHOULD BE EQUIVALENT TO HYBRID METHOD)

    // Define thrust direction settings from the MEE costates.
    std::shared_ptr< simulation_setup::MeeCostateBasedThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::MeeCostateBasedThrustDirectionSettings >( bodyToPropagate, centralBody, costatesFunction );

    // Define bang-bang thrust magnitude settings based on MEE co-states.
    std::shared_ptr< simulation_setup::FromMeeCostatesBangBangThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromMeeCostatesBangBangThrustMagnitudeSettings >(
                maximumThrust, specificImpulseFunction, costatesFunction, bodyToPropagate, centralBody );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >( thrustDirectionSettings, thrustMagnitudeSettings );


    // Define acceleration map.
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

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ bodyToPropagate ] = simulation_setup::createMassRateModel( bodyToPropagate, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap, accelerationModelMap );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

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

        // Perform propagation.
        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
        std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

        Eigen::Vector6d currentExpectedState = stateHistory.rbegin( )->second.segment( 0, 6 );
        double currentExpectedMass = stateHistory.rbegin( )->second[ 6 ];
        Eigen::Vector3d currentExpectedThrustAcceleration =  dependentVariableHistory.rbegin( )->second.segment( 0, 3 );
        Eigen::Vector3d currentExpectedThrust = currentExpectedThrustAcceleration * currentExpectedMass;

        // Check consistency between hybrid method and numerical propagation results.
        for ( int i = 0 ; i < 3 ; i++ )
        {
            BOOST_CHECK_SMALL( ( std::fabs( trajectory[ itr->first ][ i ] - currentExpectedState[ i ] ) / currentExpectedState.segment( 0, 3 ).norm( ) ), 1.0e-12 );
            BOOST_CHECK_SMALL( ( std::fabs( trajectory[ itr->first ][ i + 3 ] - currentExpectedState[ i + 3 ] ) / currentExpectedState.segment( 3, 3 ).norm( ) ), 1.0e-12 );
            BOOST_CHECK_SMALL( ( std::fabs( thrustProfile[ itr->first ][ i ] - currentExpectedThrust[ i ] ) ), 1.0e-8 );
            BOOST_CHECK_SMALL( ( std::fabs( thrustAccelerationProfile[ itr->first ][ i ] - currentExpectedThrustAcceleration[ i ] ) ), 1.0e-10 );
        }
        BOOST_CHECK_SMALL( std::fabs( massProfile[ itr->first ][ 0 ] - currentExpectedMass ), 1.0e-15 );

        initialTime = currentEpoch;
        currentMass = currentExpectedMass;
        currentState = currentExpectedState;

    }



    /// Test the computeSemiAnalyticalAndFullPropagation function
    /// (the difference between hybrid method and full propagation results should be integration errors only)

    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

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



    for ( int i = 0 ; i < 3 ; i++ )
    {
        BOOST_CHECK_SMALL( ( std::fabs(  hybridMethodResults.begin()->second[ i ] - fullPropagationResults.begin()->second[ i ] ) / stateAtDeparture.segment( 0, 3 ).norm( ) ), 2.0e-6 );
        BOOST_CHECK_SMALL( ( std::fabs(  hybridMethodResults.begin()->second[ i + 3 ] - fullPropagationResults.begin()->second[ i + 3 ] ) / stateAtDeparture.segment( 3, 3 ).norm( ) ), 2.0e-6 );

        BOOST_CHECK_SMALL( ( std::fabs(  hybridMethodResults.rbegin()->second[ i ] - fullPropagationResults.rbegin()->second[ i ] ) / stateAtArrival.segment( 0, 3 ).norm( ) ), 2.0e-6 );
        BOOST_CHECK_SMALL( ( std::fabs(  hybridMethodResults.rbegin()->second[ i + 3 ] - fullPropagationResults.rbegin()->second[ i + 3 ] ) / stateAtArrival.segment( 3, 3 ).norm( ) ), 2.0e-6 );

    }


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
