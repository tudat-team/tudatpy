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

//! Transform thrust model as a function of time into Sims Flanagan thrust model.
std::vector< Eigen::Vector3d > convertToSimsFlanaganThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime,
                                                                 const double maximumThrust,
                                                                 const double timeOfFlight, const int numberSegmentsForwardPropagation,
                                                                 const int numberSegmentsBackwardPropagation )
{
    double segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    double segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );

    std::vector< Eigen::Vector3d > SimsFlanaganThrustModel;
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        SimsFlanaganThrustModel.push_back( ( thrustModelWrtTime( segmentDurationForwardPropagation * ( i + 1.0 / 2.0 ) ) ) / maximumThrust );
    }
    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        SimsFlanaganThrustModel.push_back( ( thrustModelWrtTime( timeOfFlight / 2.0 + segmentDurationBackwardPropagation * ( i + 1.0 / 2.0 ) ) )
                                           / maximumThrust );
    }

    return SimsFlanaganThrustModel;
}


//! Perform optimisation.
std::pair< std::vector< double >, std::vector< double > > SimsFlanagan::performOptimisation( )
{
    using namespace tudat_pagmo_applications;

    //Set seed for reproducible results
    pagmo::random_device::set_seed( 456 );

    // Create object to compute the problem fitness
    problem prob{ SimsFlanaganProblem( stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulseFunction_, numberSegments_,
                                       timeOfFlight_, bodyMap_, bodyToPropagate_, centralBody_, initialGuessThrustModel_ )};



    std::vector< double > constraintsTolerance;
    for ( unsigned int i = 0 ; i < ( prob.get_nec() + prob.get_nic() ) ; i++ )
    {
        constraintsTolerance.push_back( 1.0e-3 );
    }
    prob.set_c_tol( constraintsTolerance );


    algorithm algo{ optimisationAlgorithm_ };

    unsigned long long populationSize = numberOfIndividualsPerPopulation_;

    island islUnconstrained{ algo, prob, populationSize };

    // Evolve for a given number of generations.
    for( int i = 0 ; i < numberOfGenerations_ ; i++ )
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
    }

    std::vector< double > championFitness = islUnconstrained.get_population().champion_f();
    std::vector< double > championDesignVariables = islUnconstrained.get_population().champion_x();

    std::vector< double > championFitnessConstrainedPb = prob.fitness( championDesignVariables );

    std::pair< std::vector< double >, std::vector< double > > output;
    output.first = championFitness;
    output.second = championDesignVariables;

    championFitness_ = championFitness;
    championDesignVariables_ = championDesignVariables;

    return output;

}

//! Function to compute the Sims Flanagan trajectory and the propagation fo the full problem.
void SimsFlanagan::computeSimsFlanaganTrajectoryAndFullPropagation(
     std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
     std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
     std::map< double, Eigen::VectorXd >& fullPropagationResults,
     std::map< double, Eigen::Vector6d >& SimsFlanaganResults,
     std::map< double, Eigen::VectorXd>& dependentVariablesHistory )
{

    fullPropagationResults.clear();
    SimsFlanaganResults.clear();
    dependentVariablesHistory.clear();

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

    // Compute state at half of the time of flight after forward propagation.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );
    Eigen::Vector6d stateHalfOfTimeOfFlightForwardPropagation = simsFlanaganLeg.propagateTrajectoryForward(
                0.0, timeOfFlight_ / 2.0, stateAtDeparture_, timeOfFlight_ / ( 2.0 * numberSegmentsForwardPropagation_ ) );

    // Compute state at half of the time of flight after backward propagation.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );
    Eigen::Vector6d stateHalfOfTimeOfFlightBackwardPropagation = simsFlanaganLeg.propagateTrajectoryBackward(
                timeOfFlight_, timeOfFlight_ / 2.0, stateAtArrival_, timeOfFlight_ / ( 2.0 * numberSegmentsBackwardPropagation_ ) );

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
                          std::vector< std::string >{ bodyToPropagate_ }, stateHalfOfTimeOfFlightForwardPropagation,
                          propagatorSettings.first->getTerminationSettings(),
                          propagatorSettings.first->propagator_, propagatorSettings.first->getDependentVariablesToSave() );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody_ }, accelerationMapFullPropagation,
                          std::vector< std::string >{ bodyToPropagate_ }, stateHalfOfTimeOfFlightBackwardPropagation,
                          propagatorSettings.second->getTerminationSettings(),
                          propagatorSettings.first->propagator_, propagatorSettings.second->getDependentVariablesToSave() );


    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ bodyToPropagate_ ] = createMassRateModel(
                bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                bodyMap_, accelerationMapFullPropagation );

    // Propagate mass until half of the time of flight.
    std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettingsToHalvedTimeOfFlight =
            std::make_shared< propagators::MassPropagatorSettings< double > >( std::vector< std::string >{ bodyToPropagate_ }, massRateModels,
                ( Eigen::Vector1d() << bodyMap_[ bodyToPropagate_ ]->getBodyMass() ).finished(),
                std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_ / 2.0, true ) );

    integratorSettings->initialTime_ = 0.0;

    // Create dynamics simulation object.
    propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodyMap_, integratorSettings, massPropagatorSettingsToHalvedTimeOfFlight, true, false, false );

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
    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;
    integratorSettings->initialTime_ = timeOfFlight_ / 2.0;

    // Initialise spacecraft mass to mass at halved time of flight in body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( massAtHalvedTimeOfFlight );

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards(
                bodyMap_, integratorSettings, completePropagatorSettings.first );
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
        epochsVector.push_back( itr->first );
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryBackwardPropagation[ itr->first ];
    }

    std::vector< double > epochsBackwardPropagation = epochsVector;
//    for ( int i = epochsVector.size() - 1 ; i >= 0 ; i-- )
//    {
//        epochsBackwardPropagation.push_back( epochsVector[ i ] );
//    }

    // Reset initial integrator settings.
    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;

    // Define forward propagator settings variables.
    integratorSettings->initialTime_ = timeOfFlight_ / 2.0;

    // Initialise spacecraft mass to mass at halved time of flight in body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( massAtHalvedTimeOfFlight );

    // Perform forward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards( bodyMap_, integratorSettings, completePropagatorSettings.second );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation = dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );


    std::vector< double > epochsForwardPropagation;

    // Compute and save full propagation and shaping method results along the forward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {
        epochsVector.push_back( itr->first );
        epochsForwardPropagation.push_back( itr->first );
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
    }

    std::vector< double > epochsForwardPropagationTest;
    for ( int i = epochsForwardPropagation.size() - 1 ; i >= 0 ; i-- )
    {
        epochsForwardPropagationTest.push_back( epochsForwardPropagation[ i ] );
    }

    double massAtTimeOfFlight = simsFlanaganLeg.propagateMassToSegment( numberSegments_ );

    int numberSegmentsForwardPropagation = ( numberSegments_ ) / 2;
    simsFlanaganLeg.propagateTrajectoryForward( epochsBackwardPropagation, SimsFlanaganResults, stateAtDeparture_,
                                                 initialSpacecraftMass_, 0.0, ( timeOfFlight_ / 2.0 ) / numberSegmentsForwardPropagation );

    int numberSegmentsBackwardPropagation = ( numberSegments_ + 1 ) / 2;
    simsFlanaganLeg.propagateTrajectoryBackward( epochsForwardPropagationTest, SimsFlanaganResults, stateAtArrival_,
                                                massAtTimeOfFlight, timeOfFlight_, ( timeOfFlight_ / 2.0 ) / numberSegmentsBackwardPropagation );

    deltaV_ = simsFlanaganLeg.getTotalDeltaV();

}



} // namespace low_thrust_direct_methods
} // namespace tudat
