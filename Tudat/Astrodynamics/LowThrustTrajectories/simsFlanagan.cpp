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


#include "Tudat/Astrodynamics/LowThrustTrajectories/simsFlanagan.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/simsFlanaganOptimisationSetup.h"
#include "pagmo/problems/unconstrain.hpp"
#include "pagmo/algorithms/compass_search.hpp"


namespace tudat
{
namespace low_thrust_trajectories
{

//! Transform thrust model as a function of time into Sims Flanagan thrust model.
std::vector< double > convertToSimsFlanaganThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime,
                                                        const double maximumThrust,
                                                        const double timeOfFlight, const int numberSegmentsForwardPropagation,
                                                        const int numberSegmentsBackwardPropagation )
{
    double segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    double segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );

    std::vector< double > simsFlanaganThrustModel;
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        Eigen::Vector3d currentThrustVector = thrustModelWrtTime( segmentDurationForwardPropagation * ( i + 1.0 / 2.0 ) ) / maximumThrust;
        for ( int j = 0 ; j < 3 ; j++ )
        {
            simsFlanaganThrustModel.push_back( currentThrustVector[ j ] );
        }
    }
    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        Eigen::Vector3d currentThrustVector = thrustModelWrtTime( timeOfFlight / 2.0
                                                                  + segmentDurationBackwardPropagation * ( i + 1.0 / 2.0 ) ) / maximumThrust;
        for ( int j = 0 ; j < 3 ; j++ )
        {
            simsFlanaganThrustModel.push_back( currentThrustVector[ j ] );
        }
    }

    return simsFlanaganThrustModel;
}


std::function< Eigen::Vector3d( const double ) > getInitialGuessFunctionFromShaping(
        std::shared_ptr< shape_based_methods::ShapeBasedMethod > shapeBasedLeg,
        const int numberSegmentsSimsFlanagan,
        const double timeOfFlight,
        std::function< double( const double ) > specificImpulseFunction,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    // Calculate number of segments for both the forward propagation (from departure to match point)
    // and the backward propagation (from arrival to match point).
    int numberSegmentsForwardPropagation = ( numberSegmentsSimsFlanagan + 1 ) / 2;
    int numberSegmentsBackwardPropagation = numberSegmentsSimsFlanagan / 2;
    int segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    int segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );

    // Compute time at each node of the Sims-Flanagan trajectory.
    std::vector< double > timesAtNodes;
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++)
    {
        timesAtNodes.push_back( i * segmentDurationForwardPropagation );
    }
    for ( int i = 0 ; i <= numberSegmentsBackwardPropagation ; i++ )
    {
        timesAtNodes.push_back( timeOfFlight / 2.0 + i * segmentDurationBackwardPropagation );
    }

    // Retrieve thrust profile from shaping method.
    std::map< double, Eigen::VectorXd > thrustProfileShapedTrajectory;
    shapeBasedLeg->getThrustForceProfile( timesAtNodes, thrustProfileShapedTrajectory, specificImpulseFunction, integratorSettings );

    // Compute averaged thrust vector per segment.
    std::vector< Eigen::Vector3d > thrustVectorPerSegment;
    for ( int i = 0 ; i < numberSegmentsSimsFlanagan ; i++ )
    {
        thrustVectorPerSegment.push_back( ( thrustProfileShapedTrajectory[ timesAtNodes[ i ] ]
                                          + thrustProfileShapedTrajectory[ timesAtNodes[ i + 1 ] ] ) / 2.0 );
    }

    // Define function returning Sims-Flanagan thrust vector as a function of time.
    std::function< Eigen::Vector3d( const double ) > initialGuessThrustFromShaping = [ = ]( const double currentTime )
    {

        // Compute index of the current segment, dependent on currentTime.
        int indexSegment;

        if ( currentTime <= timeOfFlight / 2.0 )
        {
            indexSegment = currentTime / segmentDurationForwardPropagation;
        }
        else if ( currentTime == timeOfFlight )
        {
            indexSegment = numberSegmentsSimsFlanagan - 1;
        }
        else
        {
            indexSegment = numberSegmentsForwardPropagation + ( currentTime - ( timeOfFlight / 2.0 ) ) / segmentDurationBackwardPropagation;
        }

        // Retrieve the constant thrust vector corresponding to this segment.
        Eigen::Vector3d thrustVectorFromInitialGuess = thrustVectorPerSegment[ indexSegment ];
        return thrustVectorFromInitialGuess;

    };

    return initialGuessThrustFromShaping;
}


//! Perform optimisation.
std::pair< std::vector< double >, std::vector< double > > SimsFlanagan::performOptimisation( )
{
    //Set seed for reproducible results
    pagmo::random_device::set_seed( 456 );

    // Create object to compute the problem fitness
    problem prob{ SimsFlanaganProblem(
                    stateAtDeparture_, stateAtArrival_, centralBodyGravitationalParameter_, initialMass_, maximumThrust_,
                    specificImpulseFunction_, numberSegments_, timeOfFlight_, initialGuessThrustModel_,
                    optimisationSettings_->relativeToleranceConstraints_ ) };



    std::vector< double > constraintsTolerance;
    for ( unsigned int i = 0 ; i < ( prob.get_nec() + prob.get_nic() ) ; i++ )
    {
        constraintsTolerance.push_back( 1.0e-3 );
    }
    prob.set_c_tol( constraintsTolerance );


    algorithm algo{ optimisationSettings_->optimisationAlgorithm_ };

    unsigned long long populationSize = optimisationSettings_->numberOfIndividualsPerPopulation_;
    island islUnconstrained{ algo, prob, populationSize };

    // Evolve for a given number of generations.
    for( int i = 0 ; i < optimisationSettings_->numberOfGenerations_ ; i++ )
    {
        islUnconstrained.evolve( );
        while( islUnconstrained.status( ) != pagmo::evolve_status::idle &&
               islUnconstrained.status( ) != pagmo::evolve_status::idle_error )
        {
            islUnconstrained.wait( );
        }
        islUnconstrained.wait_check( ); // Raises errors
    }

    std::vector< double > championFitness = islUnconstrained.get_population().champion_f();
    std::vector< double > championDesignVariables = islUnconstrained.get_population().champion_x();

    std::pair< std::vector< double >, std::vector< double > > output;
    output.first = championFitness;
    output.second = championDesignVariables;

    championFitness_ = championFitness;
    championDesignVariables_ = championDesignVariables;

    return output;

}


//! Compute direction thrust acceleration in cartesian coordinates.
Eigen::Vector3d SimsFlanagan::computeCurrentThrustAccelerationDirection(
        double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    Eigen::Vector3d currentThrustVector = computeCurrentThrustForce( currentTime, specificImpulseFunction, integratorSettings );

    Eigen::Vector3d thrustAcceleration = currentThrustVector.normalized( );

    return thrustAcceleration.normalized( );
}


//! Compute magnitude thrust acceleration.
double SimsFlanagan::computeCurrentThrustAccelerationMagnitude(
        double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    double currentMass = computeCurrentMass( currentTime, specificImpulseFunction, integratorSettings );
    Eigen::Vector3d currentThrustVector = computeCurrentThrustForce( currentTime, specificImpulseFunction, integratorSettings );

    return currentThrustVector.norm( ) / currentMass;
}


//! Compute current thrust vector.
Eigen::Vector3d SimsFlanagan::computeCurrentThrustForce(
        double currentTime,
                                                    std::function< double ( const double ) > specificImpulseFunction,
                                                    std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    int indexSegment = convertTimeToLegSegment( currentTime );
    return maximumThrust_ * simsFlanaganModel_->getThrottles( ).at( indexSegment );
}



//! Return thrust acceleration profile.
void SimsFlanagan::getThrustAccelerationProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    thrustAccelerationProfile.clear();

    std::map< double, Eigen::VectorXd > thrustProfile;
    getThrustForceProfile( epochsVector, thrustProfile, specificImpulseFunction, integratorSettings );

    std::map< double, Eigen::VectorXd > massProfile;
    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );

    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a low-thrust trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        Eigen::Vector3d currentThrustVector = thrustProfile[ epochsVector[ i ] ];
        double currentMass = massProfile[ epochsVector[ i ] ][ 0 ];

        Eigen::Vector3d currentThrustAccelerationVector = currentThrustVector / currentMass;

        thrustAccelerationProfile[ epochsVector[ i ] ] = currentThrustAccelerationVector;

    }
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
        const simulation_setup::NamedBodyMap& bodyMapTest,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::function< double ( const double ) > specificImpulseFunction,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    basic_astrodynamics::AccelerationMap SimsFlanaganAccelerationMap = simsFlanaganModel_->getLowThrustTrajectoryAccelerationMap(
                bodyMapTest, bodyToPropagate, centralBody );
    return SimsFlanaganAccelerationMap;
}


//! Define appropriate translational state propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >
SimsFlanagan::createLowThrustTranslationalStatePropagatorSettings(
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave )
{

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true );


    // Compute state at half of the time of flight after forward propagation.
    Eigen::Vector6d stateHalfOfTimeOfFlightForwardPropagation = simsFlanaganModel_->propagateTrajectoryForward(
                0.0, timeOfFlight_ / 2.0, stateAtDeparture_, timeOfFlight_ / ( 2.0 * numberSegmentsForwardPropagation_ ) );

    // Compute state at half of the time of flight after backward propagation.
    Eigen::Vector6d stateHalfOfTimeOfFlightBackwardPropagation =  simsFlanaganModel_->propagateTrajectoryBackward(
                timeOfFlight_, timeOfFlight_ / 2.0, stateAtArrival_, timeOfFlight_ / ( 2.0 * numberSegmentsBackwardPropagation_ ) );


    // Compute state at half of the time of flight (averaged between forward and backward propagation results).
    Eigen::Vector6d stateHalfOfTimeOfFlight = 1.0 / 2.0 * ( stateHalfOfTimeOfFlightForwardPropagation + stateHalfOfTimeOfFlightBackwardPropagation );

    // Define translational state propagator settings.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody }, accelerationModelMap,
              std::vector< std::string >{ bodyToPropagate }, stateHalfOfTimeOfFlightForwardPropagation,
              terminationConditions.first,
              propagators::cowell, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody }, accelerationModelMap,
              std::vector< std::string >{ bodyToPropagate }, stateHalfOfTimeOfFlightBackwardPropagation,
              terminationConditions.second,
              propagators::cowell, dependentVariablesToSave );

    return translationalStatePropagatorSettings;
}



} // namespace low_thrust_trajectories
} // namespace tudat
