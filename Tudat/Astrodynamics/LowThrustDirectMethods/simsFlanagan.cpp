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
//#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganLeg.h"

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


//! Compute direction thrust acceleration in cartesian coordinates.
Eigen::Vector3d SimsFlanagan::computeCurrentThrustAccelerationDirection( double currentTime )
{
    Eigen::Vector3d thrustAcceleration;

    basic_astrodynamics::SingleBodyAccelerationMap vehicleAccelerationMap =
            simsFlanaganLeg_->getLowThrustTrajectoryAccelerationMap( )[ bodyToPropagate_ ];

    if ( vehicleAccelerationMap[ bodyToPropagate_ ].size( ) != 1 )
    {
        throw std::runtime_error( "Error when retrieving the thrust acceleration for hybrid method, more than one acceleration exerted by"
                                  "the body Vehicle on Vehicle" );
    }
    vehicleAccelerationMap[ bodyToPropagate_ ][ 0 ]->updateMembers( currentTime );
    thrustAcceleration = vehicleAccelerationMap[ bodyToPropagate_ ][ 0 ]->getAcceleration( );

    return thrustAcceleration;
}


//! Compute magnitude thrust acceleration.
double SimsFlanagan::computeCurrentThrustAccelerationMagnitude( double currentTime )
{
    basic_astrodynamics::SingleBodyAccelerationMap vehicleAccelerationMap =
            simsFlanaganLeg_->getLowThrustTrajectoryAccelerationMap( )[ bodyToPropagate_ ];

    if ( vehicleAccelerationMap[ bodyToPropagate_ ].size( ) != 1 )
    {
        throw std::runtime_error( "Error when retrieving the thrust acceleration for hybrid method, more than one acceleration exerted by"
                                  "the body Vehicle on Vehicle" );
    }
    vehicleAccelerationMap[ bodyToPropagate_ ][ 0 ]->updateMembers( currentTime );
    Eigen::Vector3d thrustAcceleration = vehicleAccelerationMap[ bodyToPropagate_ ][ 0 ]->getAcceleration( );
    double thrustAccelerationMagnitude = thrustAcceleration.norm( );

    return thrustAccelerationMagnitude;
}


//! Compute current cartesian state.
Eigen::Vector6d SimsFlanagan::computeCurrentStateVector( const double currentTime )
{
    std::vector< double > epochsVector;
    epochsVector.push_back( 0.0 );
    epochsVector.push_back( currentTime );

    std::map< double, Eigen::Vector6d > propagatedTrajectory;
    getTrajectory( epochsVector, propagatedTrajectory );

    return propagatedTrajectory.rbegin( )->second;
}


//! Retrieve acceleration map (thrust and central gravity accelerations).
basic_astrodynamics::AccelerationMap SimsFlanagan::retrieveLowThrustAccelerationMap(
        std::function< double ( const double ) > specificImpulseFunction )
{
    basic_astrodynamics::AccelerationMap SimsFlanaganAccelerationMap = simsFlanaganLeg_->getLowThrustTrajectoryAccelerationMap( );
    return SimsFlanaganAccelerationMap;
}


//! Function to compute the Sims Flanagan trajectory and the propagation fo the full problem.
void SimsFlanagan::computeSimsFlanaganTrajectoryAndFullPropagation(
     std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
     std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
        std::shared_ptr< propagators::PropagatorSettings< double > > >& propagatorSettings,
     std::map< double, Eigen::VectorXd >& fullPropagationResults,
     std::map< double, Eigen::Vector6d >& SimsFlanaganResults,
     std::map< double, Eigen::VectorXd>& dependentVariablesHistory )
{

    fullPropagationResults.clear();
    SimsFlanaganResults.clear();
    dependentVariablesHistory.clear();

    double massAtHalvedTimeOfFlight = computeCurrentMass( timeOfFlight_ / 2.0, specificImpulseFunction_, integratorSettings );

    // Define backward propagator settings variables.
    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;
    integratorSettings->initialTime_ = timeOfFlight_ / 2.0;

    // Initialise spacecraft mass to mass at halved time of flight in body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( massAtHalvedTimeOfFlight );

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards( bodyMap_, integratorSettings, propagatorSettings.first );
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

//    std::vector< double > epochsBackwardPropagation = epochsVector;
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
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards( bodyMap_, integratorSettings, propagatorSettings.second );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation = dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );


//    std::vector< double > epochsForwardPropagation;

    // Compute and save full propagation and shaping method results along the forward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {
        epochsVector.push_back( itr->first );
//        epochsForwardPropagation.push_back( itr->first );
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
    }

    simsFlanaganLeg_->propagateTrajectory( epochsVector, SimsFlanaganResults );

    deltaV_ = simsFlanaganLeg_->getTotalDeltaV();

}



} // namespace low_thrust_direct_methods
} // namespace tudat
