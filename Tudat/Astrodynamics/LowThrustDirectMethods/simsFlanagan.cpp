/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanagan.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganLeg.h"

#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/applicationOutput.h"
#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/getAlgorithm.h"
#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/saveOptimizationResults.h"
#include "pagmo/problems/unconstrain.hpp"
#include "pagmo/algorithms/compass_search.hpp"


namespace tudat
{
namespace low_thrust_direct_methods
{


//! Perform optimisation.
std::pair< std::vector< double >, std::vector< double > > SimsFlanagan::performOptimisation( )
{
    using namespace tudat_pagmo_applications;

    //Set seed for reproducible results
    pagmo::random_device::set_seed( 456 );

    // Create object to compute the problem fitness
    problem prob{ SimsFlanaganProblem( stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulseFunction_, numberSegments_,
                                       timeOfFlight_, bodyMap_, bodyToPropagate_, centralBody_, integratorSettings_, propagatorType_,
                                       useHighOrderSolution_, optimiseTimeOfFlight_, timeOfFlightBounds_ )};



    std::vector< double > constraintsTolerance;
    for ( unsigned int i = 0 ; i < ( prob.get_nec() + prob.get_nic() ) ; i++ )
    {
        constraintsTolerance.push_back( 1.0e-3 );
    }
    prob.set_c_tol( constraintsTolerance );


//    algorithm algo{ pagmo::ihs() };
//    algo.set_verbosity( 10 );

        unconstrain unconstrainedProb{ prob, "ignore_o" };
        population pop{ unconstrainedProb, 10 };

//        algorithm algoUnconstrained{ pagmo::de1220() };
        algorithm algoUnconstrained = optimisationAlgorithm_;
        algoUnconstrained.set_verbosity( 10 );
        algoUnconstrained.evolve( pop );

        island islUnconstrained{ algoUnconstrained, unconstrainedProb, 50 /*5000*/ };

        // Evolve for 10 generations
        for( int i = 0 ; i < 1 /*600*/ ; i++ )
        {
            islUnconstrained.evolve( );
            while( islUnconstrained.status( ) != pagmo::evolve_status::idle &&
                   islUnconstrained.status( ) != pagmo::evolve_status::idle_error )
            {
                islUnconstrained.wait( );
            }
            islUnconstrained.wait_check( ); // Raises errors

            // Write current iteration results to file
            printPopulationToFile( islUnconstrained.get_population( ).get_x( ), "testSimsFlanagan_generation_" + std::to_string( i ) , false );
            printPopulationToFile( islUnconstrained.get_population( ).get_f( ), "testSimsFlanagan_generation_" + std::to_string( i ) , true );

            std::vector< double > championFitness = islUnconstrained.get_population().champion_f();
            for ( int i = 0 ; i < championFitness.size() ; i++ )
            {
                std::cout << "champion fitness: " << championFitness[ i ] << "\n\n";
            }
            std::cout << "TEST" << "\n\n";
            std::cout<< "current generation: " << i << std::endl;
        }
        std::cout << "TEST" << "\n\n";

        std::vector< double > championFitness = islUnconstrained.get_population().champion_f();
        std::vector< double > championDesignVariables = islUnconstrained.get_population().champion_x();
        for ( int i = 0 ; i < championFitness.size() ; i++ )
        {
            std::cout << "champion fitness: " << championFitness[ i ] << "\n\n";
        }
        for ( int i = 0 ; i < championDesignVariables.size() ; i++ )
        {
            std::cout << "champion design variables: " << championDesignVariables[ i ] << "\n\n";
        }

        std::vector< double > championFitnessConstrainedPb = prob.fitness( championDesignVariables );
        for ( int i = 0 ; i < championFitnessConstrainedPb.size() ; i++ )
        {
            std::cout << "champion fitness constrained problem: " << championFitnessConstrainedPb[ i ] << "\n\n";
        }

        std::cout << "champion fitness: " << championFitness[ 0 ] << "\n\n";



//        // Create an island with 1024 individuals
//        island isl{ algo, prob, 100}; //1024};

//        // Evolve for 100 generations
//        for( int i = 0 ; i < 1 ; i++ ) //300 ; i++) //100; i++ )
//        {
//            isl.evolve( );
//            while( isl.status( ) != pagmo::evolve_status::idle &&
//                   isl.status( ) != pagmo::evolve_status::idle_error )
//            {
//                isl.wait( );
//            }
//            isl.wait_check( ); // Raises errors

//            // Write current iteration results to file
//            printPopulationToFile( isl.get_population( ).get_x( ), "testSimsFlanagan_generation_" + std::to_string( i ) , false );
//            printPopulationToFile( isl.get_population( ).get_f( ), "testSimsFlanagan_generation_" + std::to_string( i ) , true );

//            std::cout<< "current generation: " << i << std::endl;
//        }

        std::pair< std::vector< double >, std::vector< double > > output;
        output.first = championFitness;
        output.second = championDesignVariables;

        championFitness_ = championFitness;
        championDesignVariables_ = championDesignVariables;

        return output;

}

//! Function to compute the Sims Flanagan trajectory and the propagation fo the full problem.
void SimsFlanagan::computeSimsFlanaganTrajectoryAndFullPropagation(
     std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
     std::map< double, Eigen::VectorXd >& fullPropagationResults,
     std::map< double, Eigen::Vector6d >& SimsFlanaganResults,
     std::map< double, Eigen::VectorXd>& dependentVariablesHistory )
{

    fullPropagationResults.clear();
    SimsFlanaganResults.clear();
    dependentVariablesHistory.clear();

//    integratorSettings_->initialTimeStep_ = integratorSettings_->initialTimeStep_ / 10.0;

    std::vector< Eigen::Vector3d > bestThrottles;
    for ( int i = 0 ; i < numberSegments_ ; i++ )
    {
        bestThrottles.push_back( ( Eigen::Vector3d( ) << championDesignVariables_[ i * 3 ],
                                 championDesignVariables_[ i * 3 + 1 ], championDesignVariables_[ i * 3 + 2 ] ).finished( ) );
    }

    // Re-initialise spacecraft mass in body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    // Create SimsFlanaganLeg object.
    SimsFlanaganLeg simsFlanaganLeg = SimsFlanaganLeg( stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulseFunction_,
                                                       timeOfFlight_, bodyMap_, bestThrottles, bodyToPropagate_,
                                                       centralBody_ );

    // Retrieve Sims Flanagan acceleration map.
    basic_astrodynamics::AccelerationMap accelerationMapSimsFlanagan = simsFlanaganLeg.getAccelerationModelFullLeg( );

    // Compute state at half of the time of flight.
    Eigen::Vector6d stateHalvedTimeOfFlight;
    if ( useHighOrderSolution_ )
    {
        stateHalvedTimeOfFlight = simsFlanaganLeg.propagateTrajectoryHighOrderSolution(
                    timeOfFlight_ / 2.0, integratorSettings_, propagatorType_ );
    }
    else
    {
//        std::cout << "state halved time of flight high fidelity solution: " << simsFlanaganLeg.propagateTrajectoryHighOrderSolution(
//                         timeOfFlight_ / 2.0, integratorSettings_, propagatorType_ ) << "\n\n";
        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );
        stateHalvedTimeOfFlight = simsFlanaganLeg.propagateTrajectory( 0.0,  timeOfFlight_ / 2.0, stateAtDeparture_ );
        std::cout << "state halved time of flight low fidelity solution: " << stateHalvedTimeOfFlight << "\n\n";
//        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );
//        simsFlanaganLeg.propagateForwardFromDepartureToMatchPoint( );
    }

    // Re-initialise spacecraft mass in body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    // Retrieve acceleration map for full propagation.
    basic_astrodynamics::AccelerationMap accelerationMapFullPropagation = propagators::getAccelerationMapFromPropagatorSettings(
                std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< double > >( propagatorSettings.first ) );

    // Add thrust acceleration model to the full propagation acceleration map.
    accelerationMapFullPropagation[ bodyToPropagate_ ][ bodyToPropagate_ ] =
            accelerationMapSimsFlanagan[ bodyToPropagate_ ][ bodyToPropagate_ ];

    // Create complete propagation settings (backward and forward propagations).
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
            std::shared_ptr< propagators::PropagatorSettings< double > > > completePropagatorSettings;


    // Define translational state propagation settings
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody_ }, accelerationMapFullPropagation,
                          std::vector< std::string >{ bodyToPropagate_ }, stateHalvedTimeOfFlight,
                          propagatorSettings.first->getTerminationSettings(),
                          propagatorType_, propagatorSettings.first->getDependentVariablesToSave() );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody_ }, accelerationMapFullPropagation,
                          std::vector< std::string >{ bodyToPropagate_ }, stateHalvedTimeOfFlight,
                          propagatorSettings.second->getTerminationSettings(),
                          propagatorType_, propagatorSettings.second->getDependentVariablesToSave() );


    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ bodyToPropagate_ ] = createMassRateModel(
                bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                bodyMap_, accelerationMapFullPropagation );

    // Propagate mass until half of the time of flight.
    std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettingsToHalvedTimeOfFlight =
            std::make_shared< propagators::MassPropagatorSettings< double > >( std::vector< std::string >{ "Vehicle" }, massRateModels,
                ( Eigen::Vector1d() << bodyMap_[ bodyToPropagate_ ]->getBodyMass() ).finished(),
                std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_ / 2.0, true ) );

    integratorSettings_->initialTime_ = 0.0;

    // Create dynamics simulation object.
    propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodyMap_, integratorSettings_, massPropagatorSettingsToHalvedTimeOfFlight, true, false, false );

    // Propagate spacecraft mass until half of the time of flight.
    std::map< double, Eigen::VectorXd > propagatedMass = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    double massAtHalvedTimeOfFlight = propagatedMass.rbegin()->second[ 0 ];

    // Create settings for propagating the mass of the vehicle.
    std::pair< std::shared_ptr< propagators::MassPropagatorSettings< double > >,
            std::shared_ptr< propagators::MassPropagatorSettings< double > > > massPropagatorSettings;

    // Define backward mass propagation settings.
    massPropagatorSettings.first = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate_ }, massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << massAtHalvedTimeOfFlight ).finished( ),
                propagatorSettings.first->getTerminationSettings() );

    // Define forward mass propagation settings.
    massPropagatorSettings.second = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate_ },
                massRateModels, ( Eigen::Matrix< double, 1, 1 >( ) << massAtHalvedTimeOfFlight ).finished( ),
                propagatorSettings.second->getTerminationSettings() );


    // Create list of propagation settings.
    std::pair< std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > >,
            std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > > propagatorSettingsVector;

    // Backward propagator settings vector.
    propagatorSettingsVector.first.push_back( translationalStatePropagatorSettings.first );
    propagatorSettingsVector.first.push_back( massPropagatorSettings.first );

    // Forward propagator settings vector.
    propagatorSettingsVector.second.push_back( translationalStatePropagatorSettings.second );
    propagatorSettingsVector.second.push_back( massPropagatorSettings.second );


    // Backward hybrid propagation settings.
    completePropagatorSettings.first = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.first,
                propagatorSettings.first->getTerminationSettings(), propagatorSettings.first->getDependentVariablesToSave() );

    // Forward hybrid propagation settings.
    completePropagatorSettings.second = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.second,
                propagatorSettings.second->getTerminationSettings(), propagatorSettings.second->getDependentVariablesToSave() );


    // Define backward propagator settings variables.
    integratorSettings_->initialTimeStep_ = - integratorSettings_->initialTimeStep_;
    integratorSettings_->initialTime_ = timeOfFlight_ / 2.0;

    // Initialise spacecraft mass to mass at halved time of flight in body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( massAtHalvedTimeOfFlight );

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards(
                bodyMap_, integratorSettings_, completePropagatorSettings.first );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation =
            dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryBackwardPropagation =
            dynamicsSimulatorIntegrationBackwards.getDependentVariableHistory( );

    // Declare vector of epochs at which the trajectory has to be calculated.
//    std::vector< double > epochsVectorBackwardPropagation;

    // Declare vector of epochs at which the trajectory has to be calculated.
    std::vector< double > epochsVector;

    // Compute and save full propagation and shaping method results along the backward propagation direction
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
    {
        if ( useHighOrderSolution_ )
        {
//            SimsFlanaganResults[ itr->first ] = simsFlanaganLeg.propagateTrajectoryHighOrderSolution(
//                        itr->first, integratorSettings_, propagatorType_ );
            epochsVector.push_back( itr->first );
        }
        else
        {
//            SimsFlanaganResults[ itr->first ] = simsFlanaganLeg.propagateTrajectory( itr->first );
            epochsVector.push_back( itr->first );
        }
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryBackwardPropagation[ itr->first ];
    }
//    if ( !useHighOrderSolution_ )
//    {
//        simsFlanaganLeg.propagateTrajectory( epochsVector, SimsFlanaganResults );
//        std::cout << "size sims flanagan results map: " << SimsFlanaganResults.size() << "\n\n";
//    }

    // Reset initial integrator settings.
    integratorSettings_->initialTimeStep_ = - integratorSettings_->initialTimeStep_;

    // Define forward propagator settings variables.
    integratorSettings_->initialTime_ = timeOfFlight_ / 2.0;

    // Initialise spacecraft mass to mass at halved time of flight in body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( massAtHalvedTimeOfFlight );

//    std::cout << "mass halved time of flight: " << massAtHalvedTimeOfFlight << "\n\n";

    // Perform forward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards( bodyMap_, integratorSettings_, completePropagatorSettings.second );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation = dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );

    // Compute and save full propagation and shaping method results along the forward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {
        if ( useHighOrderSolution_ )
        {
//            SimsFlanaganResults[ itr->first ] = simsFlanaganLeg.propagateTrajectoryHighOrderSolution(
//                        itr->first, integratorSettings_, propagatorType_ );
            epochsVector.push_back( itr->first );
        }
        else
        {
            epochsVector.push_back( itr->first );
//            SimsFlanaganResults[ itr->first ] = simsFlanaganLeg.propagateTrajectory( itr->first );
//            std::cout << "current time: " << itr->first << "\n\n";
        }
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
    }
//    if ( !useHighOrderSolution_ )
//    {
//        simsFlanaganLeg.propagateTrajectory( epochsVector, SimsFlanaganResults );
//        std::cout << "size sims flanagan results map: " << SimsFlanaganResults.size() << "\n\n";
//    }


    if ( useHighOrderSolution_ )
    {
        simsFlanaganLeg.propagateTrajectoryHighOrderSolution( epochsVector, SimsFlanaganResults, integratorSettings_, propagatorType_ );
    }
    else
    {
        simsFlanaganLeg.propagateTrajectory( epochsVector, SimsFlanaganResults );
        std::cout << "size sims flanagan results map: " << SimsFlanaganResults.size() << "\n\n";
    }


}



} // namespace low_thrust_direct_methods
} // namespace tudat
