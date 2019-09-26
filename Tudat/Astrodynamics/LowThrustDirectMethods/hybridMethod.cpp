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


#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethod.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridOptimisationSetup.h"
//#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethodLeg.h"

#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/applicationOutput.h"
#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/getAlgorithm.h"
#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/saveOptimizationResults.h"
#include "pagmo/problems/unconstrain.hpp"
#include "pagmo/algorithms/compass_search.hpp"


namespace tudat
{
namespace low_thrust_direct_methods
{


//! Transform thrust model as a function of time into hybrid method thrust model.
Eigen::Matrix< double, 10, 1 > convertToHybridMethodThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime )
{
    Eigen::Matrix< double, 10, 1 > meeInitialAndFinalCostates;

    return meeInitialAndFinalCostates;
}


//! Perform optimisation.
std::pair< std::vector< double >, std::vector< double > > HybridMethod::performOptimisation( )
{
    using namespace tudat_pagmo_applications;

    //Set seed for reproducible results
    pagmo::random_device::set_seed( 456 );

    // Create object to compute the problem fitness
    problem prob{ HybridMethodProblem( stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulse_, timeOfFlight_, bodyMap_,
                                       bodyToPropagate_, centralBody_, integratorSettings_, initialGuessThrustModel_, relativeToleranceConstraints_ )};

    std::vector< double > constraintsTolerance;
    for ( unsigned int i = 0 ; i < ( prob.get_nec() + prob.get_nic() ) ; i++ )
    {
        constraintsTolerance.push_back( 1.0e-3 );
    }
    prob.set_c_tol( constraintsTolerance );


//    unconstrain unconstrainedProb{ prob, "ignore_o" };
//    population pop{ unconstrainedProb, 10 };

    algorithm algo = optimisationAlgorithm_;
//    algoUnconstrained.set_verbosity( 10 );
//    algoUnconstrained.evolve( pop );

    unsigned long long populationSize = numberOfIndividualsPerPopulation_;

    island island{ algo, prob, populationSize };

    // Evolve for 10 generations
    for( int i = 0 ; i < numberOfGenerations_ ; i++ )
    {
        island.evolve( );
        while( island.status( ) != pagmo::evolve_status::idle &&
               island.status( ) != pagmo::evolve_status::idle_error )
        {
            island.wait( );
        }
        island.wait_check( ); // Raises errors

        // Write current iteration results to file
        printPopulationToFile( island.get_population( ).get_x( ), "testHybridMethod_generation_" + std::to_string( i ) , false );
        printPopulationToFile( island.get_population( ).get_f( ), "testHybridMethod_generation_" + std::to_string( i ) , true );

        std::vector< double > championFitness = island.get_population().champion_f();
        for ( int i = 0 ; i < championFitness.size() ; i++ )
        {
            std::cout << "champion fitness: " << championFitness[ i ] << "\n\n";
        }
        std::cout << "TEST" << "\n\n";
        std::cout<< "current generation: " << i << std::endl;
    }
    std::cout << "TEST" << "\n\n";

    std::vector< double > championFitness = island.get_population().champion_f();
    std::vector< double > championDesignVariables = island.get_population().champion_x();
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


////! Compute current mass of the spacecraft between two epochs.
//double HybridMethod::computeCurrentMass(
//        const double timeInitialEpoch,
//        const double timeFinalEpoch,
//        const double massInitialEpoch,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    // Retrieve acceleration map.
//    basic_astrodynamics::AccelerationMap accelerationMap = retrieveLowThrustAccelerationMap( specificImpulseFunction ); // hybridMethodLeg_->getLowThrustTrajectoryAccelerationMap( );

//    // Create mass rate models
//    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
//    massRateModels[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
//                                                       bodyMap_, accelerationMap );

//    // Define mass propagator settings.
//    std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettings =
//            std::make_shared< propagators::MassPropagatorSettings< double > >( std::vector< std::string >{ bodyToPropagate_ }, massRateModels,
//                ( Eigen::Vector1d() << massInitialEpoch ).finished(),
//                std::make_shared< propagators::PropagationTimeTerminationSettings >( timeFinalEpoch, true ) );

//    // Re-initialise integrator settings.
//    integratorSettings->initialTime_ = timeInitialEpoch;
//    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

//    // Create dynamics simulation object.
//    propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
//                bodyMap_, integratorSettings, massPropagatorSettings, true, false, false );

//    // Propagate spacecraft mass.
//    std::map< double, Eigen::VectorXd > propagatedMass = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
//    double currentMass = propagatedMass.rbegin()->second[ 0 ];

//    return currentMass;

//}


////! Return mass profile.
//void HybridMethod::getMassProfile(
//        std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& massProfile,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    massProfile.clear( );

//    double currentMass = initialMass_;

//    for ( int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
//        {
//            throw std::runtime_error( "Error when retrieving the mass profile of a hybrid trajectory, "
//                                      "epochs are not provided in increasing order." );
//        }

//        if ( i == 0 )
//        {
//            currentMass = computeCurrentMass( 0.0, epochsVector[ i ], currentMass, specificImpulseFunction, integratorSettings );
//            massProfile[ epochsVector[ i ] ] = ( Eigen::Vector1d( ) << currentMass ).finished( );
//        }
//        else
//        {
//            currentMass = computeCurrentMass( epochsVector[ i - 1 ], epochsVector[ i ], currentMass, specificImpulseFunction, integratorSettings );
//            massProfile[ epochsVector[ i ] ] = ( Eigen::Vector1d( ) << currentMass ).finished();
//        }
//    }

//}


//! Compute magnitude thrust acceleration.
double HybridMethod::computeCurrentThrustAccelerationMagnitude( double currentTime )
{
    basic_astrodynamics::SingleBodyAccelerationMap vehicleAccelerationMap =
            hybridMethodLeg_->getLowThrustTrajectoryAccelerationMap( )[ bodyToPropagate_ ];

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


//! Compute direction thrust acceleration in cartesian coordinates.
Eigen::Vector3d HybridMethod::computeCurrentThrustAccelerationDirection( double currentTime )
{
    Eigen::Vector3d thrustAcceleration;

    basic_astrodynamics::SingleBodyAccelerationMap vehicleAccelerationMap =
            hybridMethodLeg_->getLowThrustTrajectoryAccelerationMap( )[ bodyToPropagate_ ];

    if ( vehicleAccelerationMap[ bodyToPropagate_ ].size( ) != 1 )
    {
        throw std::runtime_error( "Error when retrieving the thrust acceleration for hybrid method, more than one acceleration exerted by"
                                  "the body Vehicle on Vehicle" );
    }
    vehicleAccelerationMap[ bodyToPropagate_ ][ 0 ]->updateMembers( currentTime );
    thrustAcceleration = vehicleAccelerationMap[ bodyToPropagate_ ][ 0 ]->getAcceleration( );

    return thrustAcceleration;
}

//! Compute current cartesian state.
Eigen::Vector6d HybridMethod::computeCurrentStateVector( const double currentTime )
{

    return hybridMethodLeg_->propagateTrajectory( 0.0, currentTime, stateAtDeparture_, initialSpacecraftMass_/*, integratorSettings_*/ );

//    std::vector< double > epochsVector;
////    epochsVector.push_back( 0.0 );
//    epochsVector.push_back( currentTime );

//    std::map< double, Eigen::Vector6d > propagatedTrajectory;
//    getTrajectory( epochsVector, propagatedTrajectory );

//    std::cout << "size propagated trajectory: " << propagatedTrajectory.size(  ) << "\n\n";

//    return propagatedTrajectory.rbegin( )->second;
}


////! Compute current thrust vector.
//Eigen::Vector3d HybridMethod::computeCurrentThrustAcceleration( double time )
//{

//    double independentVariable = convertTimeToIndependentVariable( time );
//    return computeCurrentThrustAccelerationMagnitude( independentVariable ) * computeCurrentThrustAccelerationDirection( independentVariable );
//}


////! Return thrust acceleration profile.
//void HybridMethod::getThrustAccelerationProfile(
//        std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& thrustAccelerationProfile )
//{
//    thrustAccelerationProfile.clear();

//    for ( int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
//        {
//            throw std::runtime_error( "Error when retrieving the thrust profile of a shape-based trajectories, "
//                                      "epochs are not provided in increasing order." );
//        }

//        Eigen::Vector3d currentThrustAccelerationVector = computeCurrentThrustAcceleration( epochsVector[ i ] );
//        thrustAccelerationProfile[ epochsVector[ i ] ] = currentThrustAccelerationVector;

//    }

//}


////! Compute current thrust vector.
//Eigen::Vector3d HybridMethod::computeCurrentThrust( double time,
//                                                    std::function< double ( const double ) > specificImpulseFunction,
//                                                    std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    double independentVariable = convertTimeToIndependentVariable( time );
//    return computeCurrentMass( 0.0, time, initialSpacecraftMass_, specificImpulseFunction, integratorSettings )
//            * computeCurrentThrustAccelerationMagnitude( independentVariable ) * computeCurrentThrustAccelerationDirection( independentVariable );
//}


////! Return thrust profile.
//void HybridMethod::getThrustProfile(
//        std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& thrustProfile,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings)
//{
//    thrustProfile.clear( );

//    // Retrieve corresponding mass profile.
//    std::map< double, Eigen::VectorXd > massProfile;
//    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );
//    std::vector< double > massProfileVector;
//    for ( std::map< double, Eigen::VectorXd >::iterator itr = massProfile.begin( ) ; itr != massProfile.end( ) ; itr++ )
//    {
//        massProfileVector.push_back( itr->second[ 0 ] );
//    }

//    for ( int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
//        {
//            throw std::runtime_error( "Error when retrieving the thrust profile of a shape-based trajectories, "
//                                      "epochs are not provided in increasing order." );
//        }

//        Eigen::Vector3d currentThrustVector = computeCurrentThrustAcceleration( epochsVector[ i ] ) * massProfileVector[ i ];
//        thrustProfile[ epochsVector[ i ] ] = currentThrustVector;

//    }
//}


//! Retrieve acceleration map (thrust and central gravity accelerations).
basic_astrodynamics::AccelerationMap HybridMethod::retrieveLowThrustAccelerationMap(
        std::function< double ( const double ) > specificImpulseFunction )
{
    basic_astrodynamics::AccelerationMap hybridMethodAccelerationMap = hybridMethodLeg_->getLowThrustTrajectoryAccelerationMap( );
    return hybridMethodAccelerationMap;
}


//! Define appropriate translational state propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > HybridMethod::createLowThrustTranslationalStatePropagatorSettings(
        basic_astrodynamics::AccelerationMap accelerationModelMap,
        std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave )
{

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true );

    // Compute state vector at half of the time of flight.
    double independentVariableAtHalfTimeOfFlight = convertTimeToIndependentVariable( timeOfFlight_ / 2.0 );
    Eigen::Vector6d stateAtHalfOfTimeOfFlight = computeCurrentStateVector( timeOfFlight_ / 2.0 );

    // Define translational state propagator settings.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody_ }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate_ }, stateAtHalfOfTimeOfFlight,
                          terminationConditions.first, propagators::gauss_modified_equinoctial, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( std::vector< std::string >{ centralBody_ }, accelerationModelMap,
                          std::vector< std::string >{ bodyToPropagate_ }, stateAtHalfOfTimeOfFlight,
                          terminationConditions.second, propagators::gauss_modified_equinoctial, dependentVariablesToSave );

    return translationalStatePropagatorSettings;
}


} // namespace low_thrust_direct_methods
} // namespace tudat
