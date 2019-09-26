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


//! Define appropriate translational state propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > SimsFlanagan::createLowThrustTranslationalStatePropagatorSettings(
        basic_astrodynamics::AccelerationMap accelerationModelMap,
        std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave )
{

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true );


    // Compute state at half of the time of flight after forward propagation.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );
    Eigen::Vector6d stateHalfOfTimeOfFlightForwardPropagation = simsFlanaganLeg_->propagateTrajectoryForward(
                0.0, timeOfFlight_ / 2.0, stateAtDeparture_, timeOfFlight_ / ( 2.0 * numberSegmentsForwardPropagation_ ) );

    // Compute state at half of the time of flight after backward propagation.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );
    Eigen::Vector6d stateHalfOfTimeOfFlightBackwardPropagation =  simsFlanaganLeg_->propagateTrajectoryBackward(
                timeOfFlight_, timeOfFlight_ / 2.0, stateAtArrival_, timeOfFlight_ / ( 2.0 * numberSegmentsBackwardPropagation_ ) );


    // Compute state at half of the time of flight (averaged between forward and backward propagation results).
    Eigen::Vector6d stateHalfOfTimeOfFlight = 1.0 / 2.0 * ( stateHalfOfTimeOfFlightForwardPropagation + stateHalfOfTimeOfFlightBackwardPropagation );

    // Re-initialise spacecraft mass in body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );


    // Define translational state propagator settings.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody_ }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate_ }, /*stateHalfOfTimeOfFlight*/ stateHalfOfTimeOfFlightForwardPropagation,
                          terminationConditions.first,
                          propagators::cowell, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody_ }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate_ }, /*stateHalfOfTimeOfFlight*/ stateHalfOfTimeOfFlightBackwardPropagation,
                          terminationConditions.second,
                          propagators::cowell, dependentVariablesToSave );

    return translationalStatePropagatorSettings;
}



} // namespace low_thrust_direct_methods
} // namespace tudat
