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

//#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/applicationOutput.h"
//#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/getAlgorithm.h"
//#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/saveOptimizationResults.h"
#include "pagmo/problems/unconstrain.hpp"
#include "pagmo/algorithms/compass_search.hpp"
#include "Tudat/Astrodynamics/Propulsion/thrustMagnitudeWrapper.h"
#include "Tudat/Astrodynamics/Propulsion/costateBasedThrustGuidance.h"


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
//    using namespace tudat_pagmo_applications;

    //Set seed for reproducible results
    pagmo::random_device::set_seed( 456 );

    // Create object to compute the problem fitness
    problem prob{ HybridMethodProblem( stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulse_, timeOfFlight_, bodyMap_,
                                       bodyToPropagate_, centralBody_, integratorSettings_, initialGuessThrustModel_, optimisationSettings_->relativeToleranceConstraints_
                                       /*relativeToleranceConstraints_*/ )};

    std::vector< double > constraintsTolerance;
    for ( unsigned int i = 0 ; i < ( prob.get_nec() + prob.get_nic() ) ; i++ )
    {
        constraintsTolerance.push_back( 1.0e-3 );
    }
    prob.set_c_tol( constraintsTolerance );


//    unconstrain unconstrainedProb{ prob, "ignore_o" };
//    population pop{ unconstrainedProb, 10 };

    algorithm algo = optimisationSettings_->optimisationAlgorithm_; // optimisationAlgorithm_;
//    algoUnconstrained.set_verbosity( 10 );
//    algoUnconstrained.evolve( pop );

    unsigned long long populationSize = optimisationSettings_->numberOfIndividualsPerPopulation_; // numberOfIndividualsPerPopulation_;

    island island{ algo, prob, populationSize };

    // Evolve for 10 generations
    for( int i = 0 ; i < optimisationSettings_->numberOfGenerations_ /*numberOfGenerations_*/ ; i++ )
    {
        island.evolve( );
        while( island.status( ) != pagmo::evolve_status::idle &&
               island.status( ) != pagmo::evolve_status::idle_error )
        {
            island.wait( );
        }
        island.wait_check( ); // Raises errors

//        // Write current iteration results to file
//        printPopulationToFile( island.get_population( ).get_x( ), "testHybridMethod_generation_" + std::to_string( i ) , false );
//        printPopulationToFile( island.get_population( ).get_f( ), "testHybridMethod_generation_" + std::to_string( i ) , true );

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


//! Compute current thrust vector.
Eigen::Vector3d HybridMethod::computeCurrentThrust( double time,
                                                    std::function< double ( const double ) > specificImpulseFunction,
                                                    std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{

//    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

   std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings = hybridMethodLeg_->getMEEcostatesBasedThrustMagnitudeSettings( );

   std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection = simulation_setup::getBodyFixedThrustDirection(
           thrustMagnitudeSettings, bodyMap_, bodyToPropagate_ );

   std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction = [ = ] ( )
   {
     return computeCurrentStateVector( time );
   };

   std::function< Eigen::Vector6d( ) > centralBodyStateFunction = [ = ] ( )
   {
     return Eigen::Vector6d::Zero( );
   };

   std::function< double( ) > centralBodyGravitationalParameterFunction = [ = ]( )
   {
        return bodyMap_[ centralBody_ ]->getGravityFieldModel( )->getGravitationalParameter( );
   };

   std::function< double( ) > thrustingBodyMassFunction = std::bind( &simulation_setup::Body::getBodyMass, bodyMap_.at( bodyToPropagate_ ) );


   propulsion::MeeCostatesBangBangThrustMagnitudeWrapper thrustMagnitudeWrapper = propulsion::MeeCostatesBangBangThrustMagnitudeWrapper(
           thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
           hybridMethodLeg_->getCostatesFunction_( ), maximumThrust_, specificImpulseFunction, thrustingBodyMassFunction);

   propulsion::MeeCostateBasedThrustGuidance thrustGuidance = propulsion::MeeCostateBasedThrustGuidance(
               thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
               hybridMethodLeg_->getCostatesFunction_( ), bodyFixedThrustDirection );

   thrustGuidance.updateCalculator( time );
   thrustMagnitudeWrapper.update( time );

   return thrustMagnitudeWrapper.getCurrentThrustMagnitude( ) * thrustGuidance.getCurrentForceDirectionInPropagationFrame( );
}



//! Return thrust profile.
void HybridMethod::getThrustProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    thrustProfile.clear( );
//    std::map< double, Eigen::VectorXd > massProfile;

//    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );

    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings = hybridMethodLeg_->getMEEcostatesBasedThrustMagnitudeSettings( );

    std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection = simulation_setup::getBodyFixedThrustDirection(
            thrustMagnitudeSettings, bodyMap_, bodyToPropagate_ );

    std::map< double, Eigen::Vector6d > trajectory;
    getTrajectory( epochsVector, trajectory );

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a shape-based trajectories, "
                                      "epochs are not provided in increasing order." );
        }

//        double independentVariable = convertTimeToIndependentVariable( epochsVector[ i ] );

//        double currentMass = massProfile[ epochsVector[ i ] ][ 0 ];
//        thrustProfile[ epochsVector[ i ] ] = currentMass * computeCurrentThrustAccelerationMagnitude( independentVariable, specificImpulseFunction, integratorSettings )
//                * computeCurrentThrustAccelerationDirection( independentVariable, specificImpulseFunction, integratorSettings );


        Eigen::Vector6d currentStateVector = trajectory[ epochsVector[ i ] ];


        std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction = [ = ] ( )
        {
            return currentStateVector;
        };

        std::function< Eigen::Vector6d( ) > centralBodyStateFunction = [ = ] ( )
        {
          return Eigen::Vector6d::Zero( );
        };

        std::function< double( ) > centralBodyGravitationalParameterFunction = [ = ]( )
        {
             return bodyMap_[ centralBody_ ]->getGravityFieldModel( )->getGravitationalParameter( );
        };

        std::function< double( ) > thrustingBodyMassFunction = std::bind( &simulation_setup::Body::getBodyMass, bodyMap_.at( bodyToPropagate_ ) );


        propulsion::MeeCostatesBangBangThrustMagnitudeWrapper thrustMagnitudeWrapper = propulsion::MeeCostatesBangBangThrustMagnitudeWrapper(
                thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                hybridMethodLeg_->getCostatesFunction_( ), maximumThrust_, specificImpulseFunction, thrustingBodyMassFunction);

        propulsion::MeeCostateBasedThrustGuidance thrustGuidance = propulsion::MeeCostateBasedThrustGuidance(
                    thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                    hybridMethodLeg_->getCostatesFunction_( ), bodyFixedThrustDirection );

        thrustGuidance.updateCalculator( epochsVector[ i ] );
        thrustMagnitudeWrapper.update( epochsVector[ i ] );

        thrustProfile[ epochsVector[ i ] ] = thrustMagnitudeWrapper.getCurrentThrustMagnitude( ) * thrustGuidance.getCurrentForceDirectionInPropagationFrame( );

    }
}

//! Compute magnitude thrust acceleration.
double HybridMethod::computeCurrentThrustAccelerationMagnitude(
        double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{

//    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
//    {
//        return specificImpulse_;
//    };

    double currentMass = computeCurrentMass( currentTime, specificImpulseFunction, integratorSettings );
    Eigen::Vector3d currentThrustVector = computeCurrentThrust( currentTime, specificImpulseFunction, integratorSettings_ );

    return currentThrustVector.norm( ) / currentMass;

}


//! Compute direction thrust acceleration in cartesian coordinates.
Eigen::Vector3d HybridMethod::computeCurrentThrustAccelerationDirection(
        double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{

//    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
//    {
//        return specificImpulse_;
//    };

    Eigen::Vector3d currentThrustVector = computeCurrentThrust( currentTime, specificImpulseFunction, integratorSettings );

    Eigen::Vector3d thrustAcceleration = currentThrustVector.normalized( );

    return thrustAcceleration.normalized( );
}



//! Return thrust acceleration profile.
void HybridMethod::getThrustAccelerationProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    thrustAccelerationProfile.clear();

    std::map< double, Eigen::VectorXd > thrustProfile;
    getThrustProfile( epochsVector, thrustProfile, specificImpulseFunction, integratorSettings );

    std::map< double, Eigen::VectorXd > massProfile;
    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
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
Eigen::Vector6d HybridMethod::computeCurrentStateVector( const double currentTime )
{
    Eigen::Vector6d stateVector;
    if ( currentTime == 0.0 )
    {
        stateVector = stateAtDeparture_;
    }
    else
    {
        stateVector = hybridMethodLeg_->propagateTrajectory( 0.0, currentTime, stateAtDeparture_, initialSpacecraftMass_/*, integratorSettings_*/ );
    }

    return stateVector;

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
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
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
